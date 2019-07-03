#include "denovo_discovery/candidate_region.h"


size_t std::hash<CandidateRegionIdentifier>::operator()(const CandidateRegionIdentifier &id) const {
    return std::hash<int>()(std::get<0>(id).start) ^ std::hash<int>()(std::get<0>(id).get_end()) ^
           std::hash<std::string>()(std::get<1>(id));
}


CandidateRegion::CandidateRegion(const Interval &interval, std::string name)
        : interval { interval }, name { std::move(name) }, interval_padding { 0 } {
    initialise_filename();
    omp_init_lock(&add_pileup_entry_lock);
}


CandidateRegion::CandidateRegion(const Interval &interval, std::string name, const uint_least16_t &interval_padding)
        : interval { interval }, name { std::move(name) }, interval_padding { interval_padding } {
    initialise_filename();
    omp_init_lock(&add_pileup_entry_lock);
}

CandidateRegion::~CandidateRegion() {
    omp_destroy_lock(&add_pileup_entry_lock);
}

void CandidateRegion::initialise_filename() {
    const auto i { get_interval() };
    this->filename = fs::path(std::string(get_name()).append(".").append(std::to_string(i.start)).append("-")
                                                     .append(std::to_string(i.get_end())).append("_denovo_discovery")
                                                     .append(".fa"));
}


Interval CandidateRegion::get_interval() const {
    const auto interval_start { (interval.start <= interval_padding) ? 0 : interval.start - interval_padding };
    const auto interval_end { interval.get_end() + interval_padding };
    return { interval_start, interval_end };
}


const std::string &CandidateRegion::get_name() const {
    return name;
}


CandidateRegionIdentifier CandidateRegion::get_id() const {
    return std::make_tuple(this->get_interval(), this->get_name());
}


std::string CandidateRegion::get_max_likelihood_sequence_with_flanks() const {
    return std::string(left_flanking_sequence).append(max_likelihood_sequence).append(right_flanking_sequence);
}


bool CandidateRegion::operator==(const CandidateRegion &rhs) const {
    return interval == rhs.interval && name == rhs.name;
}


bool CandidateRegion::operator!=(const CandidateRegion &rhs) const {
    return !(rhs == *this);
}


CandidateRegions find_candidate_regions_for_pan_node(const TmpPanNode &pangraph_node_components,
                                                     const uint_least16_t &candidate_region_interval_padding) {
    const auto &pangraph_node { pangraph_node_components.pangraph_node };
    const auto &local_prg { pangraph_node_components.local_prg };
    const auto &local_node_max_likelihood_path { pangraph_node_components.local_node_max_likelihood_path };
    const auto &kmer_node_max_likelihood_path { pangraph_node_components.kmer_node_max_likelihood_path };
    const uint32_t sample_id { 0 };

    const auto covgs_along_localnode_path {
            get_covgs_along_localnode_path(pangraph_node, local_node_max_likelihood_path, kmer_node_max_likelihood_path,
                                           sample_id) };
    const auto candidate_intervals { identify_low_coverage_intervals(covgs_along_localnode_path) };

    CandidateRegions candidate_regions;

    BOOST_LOG_TRIVIAL(debug) << "there are " << candidate_intervals.size() << " intervals";

    for (const auto &current_interval: candidate_intervals) {
        CandidateRegion candidate_region { current_interval, pangraph_node->get_name(),
                                           candidate_region_interval_padding };
        BOOST_LOG_TRIVIAL(debug) << "Looking at interval: " << candidate_region.get_interval();
        BOOST_LOG_TRIVIAL(debug) << "For gene: " << candidate_region.get_name();

        const auto interval_path_components { find_interval_and_flanks_in_localpath(candidate_region.get_interval(),
                                                                                    local_node_max_likelihood_path) };

        candidate_region.read_coordinates = pangraph_node->get_read_overlap_coordinates(interval_path_components.slice);

        BOOST_LOG_TRIVIAL(debug) << "there are " << candidate_region.read_coordinates.size()
                                 << " read_overlap_coordinates";

        candidate_region.max_likelihood_sequence = local_prg->string_along_path(interval_path_components.slice);
        candidate_region.left_flanking_sequence = local_prg->string_along_path(interval_path_components.flank_left);
        candidate_region.right_flanking_sequence = local_prg->string_along_path(interval_path_components.flank_right);
        candidate_regions.insert(std::make_pair(candidate_region.get_id(), candidate_region));
    }
    return candidate_regions;
}


std::vector<Interval>
identify_low_coverage_intervals(const std::vector<uint32_t> &covg_at_each_position, const uint32_t &min_required_covg,
                                const uint32_t &min_length) {
    std::vector<Interval> identified_regions;
    const auto predicate { [min_required_covg](uint_least32_t x) { return x <= min_required_covg; }};
    auto first { covg_at_each_position.begin() };
    auto previous { covg_at_each_position.begin() };
    auto current { first };

    //TODO: merging candidate regions that are close enough together
    while (current != covg_at_each_position.end()) {
        previous = current;
        current = std::find_if_not(current, covg_at_each_position.end(), predicate);
        if (current - previous >= min_length) {
            identified_regions.emplace_back(previous - first, current - first);
        }
        if (current == covg_at_each_position.end()) {
            break;
        }
        ++current;
    }
    return identified_regions;
}

void CandidateRegion::add_pileup_entry(const std::string &read, const ReadCoordinate &read_coordinate) {
    const bool read_coord_start_is_past_read_end { read_coordinate.start >= read.length() };
    if (not read_coord_start_is_past_read_end) {
        const auto end_pos_of_region_in_read{std::min(read_coordinate.end, (uint32_t) read.length())};
        std::string sequence_in_read_overlapping_region{
                read.substr(read_coordinate.start, end_pos_of_region_in_read - read_coordinate.start)};

        if (!read_coordinate.is_forward) {
            sequence_in_read_overlapping_region = reverse_complement(sequence_in_read_overlapping_region);
        }

        omp_set_lock(&add_pileup_entry_lock);
            this->pileup.push_back(sequence_in_read_overlapping_region);
        omp_unset_lock(&add_pileup_entry_lock);
    }
}


void CandidateRegion::write_denovo_paths_to_file(const fs::path &output_directory) {
    fs::create_directories(output_directory);

    if (denovo_paths.empty()) {
        BOOST_LOG_TRIVIAL(debug) << "No denovo paths for " << filename;
        return;
    }

    auto fasta { generate_fasta_for_denovo_paths() };

    const auto discovered_paths_filepath { output_directory / filename };

    fasta.save(discovered_paths_filepath.string());

    BOOST_LOG_TRIVIAL(debug) << "Denovo path for " << name << " written to: " << discovered_paths_filepath;

}


Fastaq CandidateRegion::generate_fasta_for_denovo_paths() {
    const bool gzip { false };
    const bool fastq { false };
    Fastaq fasta { gzip, fastq };

    for (size_t i = 0; i < denovo_paths.size(); ++i) {
        const auto &path { denovo_paths.at(i) };
        const std::string path_name { std::string(name).append(".").append(std::to_string(i)) };
        fasta.add_entry(path_name, path, "");
    }

    return fasta;
}


TmpPanNode::TmpPanNode(const PanNodePtr &pangraph_node, const shared_ptr<LocalPRG> &local_prg,
                       const vector<KmerNodePtr> &kmer_node_max_likelihood_path,
                       const vector<LocalNodePtr> &local_node_max_likelihood_path)
        : pangraph_node { pangraph_node }, local_prg { local_prg },
          kmer_node_max_likelihood_path { kmer_node_max_likelihood_path },
          local_node_max_likelihood_path { local_node_max_likelihood_path } {}




PileupConstructionMap construct_pileup_construction_map(CandidateRegions &candidate_regions) {
    PileupConstructionMap pileup_construction_map;
    for (auto &element : candidate_regions) {
        auto &candidate_region{element.second};
        for (const auto &read_coordinate : candidate_region.read_coordinates) {
            pileup_construction_map[read_coordinate.id].emplace_back(&candidate_region, &read_coordinate);
        }
    }
    return pileup_construction_map;
}



void
load_all_candidate_regions_pileups_from_fastq(const fs::path &reads_filepath, const CandidateRegions &candidate_regions,
                                              const PileupConstructionMap &pileup_construction_map, const uint32_t threads) {
    if (candidate_regions.empty() or pileup_construction_map.empty())
        return;

    const uint32_t nb_reads_to_map_in_a_batch = 1000; //nb of reads to map in a batch

    //shared variables - controlled by critical(ReadFileMutex)
    FastaqHandler fh(reads_filepath.string());
    uint32_t id{0};

    //parallel region
    //TODO: this is duplicated code with pangraph_from_read_file(), refactor
    #pragma omp parallel num_threads(threads)
    {
        //will hold the reads batch
        std::vector<std::pair<uint32_t, std::string>> sequencesBuffer(nb_reads_to_map_in_a_batch, make_pair(0, ""));
        while (true) {
            //read the next batch of reads
            uint32_t nbOfReads = 0;

            //read the reads in batch
            #pragma omp critical(ReadFileMutex)
            {
                //TODO: we need to read only until the max read id
                for (auto &id_and_sequence : sequencesBuffer) {
                    //did we reach the end already?
                    if (!fh.eof()) { //no
                        //print some logging
                        if (id && id % 100000 == 0)
                            BOOST_LOG_TRIVIAL(info) << id << " reads processed...";

                        //read the read
                        fh.get_next();
                        id_and_sequence = make_pair(id, fh.read);
                        ++nbOfReads;
                        ++id;
                    } else { //yes
                        break; //we read everything already, exit this loop
                    }
                }
            }


            if (nbOfReads == 0) break; //we reached the end of the file, nothing else to map

            //process nbOfReads reads
            for (uint32_t i=0; i<nbOfReads; i++) {
                uint32_t id;
                std::string sequence;
                std::tie(id, sequence)  = sequencesBuffer[i];

                auto pileup_construction_map_iterator = pileup_construction_map.find(id);
                const bool is_read_required_for_pileup = pileup_construction_map_iterator != pileup_construction_map.end();

                if (not is_read_required_for_pileup)
                    continue;

                //create all pileups for this read
                for (const auto &pair : pileup_construction_map_iterator->second) {
                    CandidateRegion * candidate_region;
                    const ReadCoordinate * read_coordinate;
                    std::tie(candidate_region, read_coordinate) = pair;

                    candidate_region->add_pileup_entry(sequence, *read_coordinate);
                }
            }

        }
    }
}