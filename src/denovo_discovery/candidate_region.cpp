#include "denovo_discovery/candidate_region.h"
#include "utils.h"

size_t std::hash<CandidateRegionIdentifier>::operator()(
    const CandidateRegionIdentifier& id) const
{
    return std::hash<int>()(std::get<0>(id).start)
        ^ std::hash<int>()(std::get<0>(id).get_end())
        ^ std::hash<std::string>()(std::get<1>(id));
}

CandidateRegion::CandidateRegion(const Interval& interval, std::string name)
    : interval { interval }
    , name { std::move(name) }
    , interval_padding { 0 }
{
    initialise_filename();
#ifndef NO_OPENMP
    omp_init_lock(&add_pileup_entry_lock);
#endif
}

CandidateRegion::CandidateRegion(
    const Interval& interval, std::string name, const uint_least16_t& interval_padding,
    const std::vector<LocalNodePtr> &local_node_max_likelihood_path,
    const std::string &local_node_max_likelihood_sequence)
    : interval { interval }
    , name { std::move(name) }
    , interval_padding { interval_padding }
    , local_node_max_likelihood_path(local_node_max_likelihood_path)
    , local_node_max_likelihood_sequence { local_node_max_likelihood_sequence }
{
    initialise_filename();
#ifndef NO_OPENMP
    omp_init_lock(&add_pileup_entry_lock);
#endif
}

CandidateRegion::~CandidateRegion()
{
#ifndef NO_OPENMP
    omp_destroy_lock(&add_pileup_entry_lock);
#endif
}

void CandidateRegion::initialise_filename()
{
    const auto i { get_interval() };
    this->filename = fs::path(std::string(get_name())
                                  .append(".")
                                  .append(std::to_string(i.start))
                                  .append("-")
                                  .append(std::to_string(i.get_end()))
                                  .append("_denovo_discovery")
                                  .append(".fa"));
}

Interval CandidateRegion::get_interval() const
{
    const auto interval_start {
        (interval.start <= interval_padding) ? 0 : interval.start - interval_padding
    };
    const auto interval_end { interval.get_end() + interval_padding };
    return { interval_start, interval_end };
}

const std::string& CandidateRegion::get_name() const { return name; }

CandidateRegionIdentifier CandidateRegion::get_id() const
{
    return std::make_tuple(this->get_interval(), this->get_name());
}

std::string CandidateRegion::get_max_likelihood_sequence_with_flanks() const
{
    return std::string(left_flanking_sequence)
        .append(max_likelihood_sequence)
        .append(right_flanking_sequence);
}

bool CandidateRegion::operator==(const CandidateRegion& rhs) const
{
    return interval == rhs.interval && name == rhs.name;
}

bool CandidateRegion::operator!=(const CandidateRegion& rhs) const
{
    return !(rhs == *this);
}

CandidateRegions Discover::find_candidate_regions_for_pan_node(
    const TmpPanNode& pangraph_node_components)
{
    const auto& pangraph_node { pangraph_node_components.pangraph_node };
    const auto& local_prg { pangraph_node_components.local_prg };
    const auto& local_node_max_likelihood_path {
        pangraph_node_components.local_node_max_likelihood_path
    };
    const auto& local_node_max_likelihood_seq {
        pangraph_node_components.local_node_max_likelihood_seq
    };
    const auto& kmer_node_max_likelihood_path {
        pangraph_node_components.kmer_node_max_likelihood_path
    };
    const uint32_t sample_id { 0 };

    const auto covgs_along_localnode_path { get_covgs_along_localnode_path(
        pangraph_node, local_node_max_likelihood_path, kmer_node_max_likelihood_path,
        sample_id) };

    uint i { 1 };
    BOOST_LOG_TRIVIAL(trace) << "Coverage along localnode path";
    for (const auto& cov : covgs_along_localnode_path) {
        BOOST_LOG_TRIVIAL(trace) << i << ": " << cov;
        i++;
    }

    auto candidate_intervals { this->identify_low_coverage_intervals(
        covgs_along_localnode_path) };
    BOOST_LOG_TRIVIAL(trace) << candidate_intervals.size()
                             << " candidate intervals before merging";
    merge_intervals_within(candidate_intervals, this->merge_dist);
    BOOST_LOG_TRIVIAL(trace) << candidate_intervals.size()
                             << " candidate intervals after merging";

    CandidateRegions candidate_regions;

    for (const auto& current_interval : candidate_intervals) {
        CandidateRegion candidate_region { current_interval, pangraph_node->get_name(),
            this->candidate_padding , local_node_max_likelihood_path, local_node_max_likelihood_seq};

        auto interval_path_components { find_interval_and_flanks_in_localpath(
            candidate_region.get_interval(), local_node_max_likelihood_path) };

        candidate_region.read_coordinates = pangraph_node->get_read_overlap_coordinates(
            interval_path_components.slice);

        BOOST_LOG_TRIVIAL(trace)
            << "Candidate region with interval " << candidate_region.get_interval()
            << " for gene " << candidate_region.get_name() << " has "
            << candidate_region.read_coordinates.size() << " read overlaps";

        candidate_region.max_likelihood_sequence
            = local_prg->string_along_path(interval_path_components.slice);
        candidate_region.left_flanking_sequence
            = local_prg->string_along_path(interval_path_components.flank_left);
        candidate_region.right_flanking_sequence
            = local_prg->string_along_path(interval_path_components.flank_right);
        candidate_regions.insert(
            std::make_pair(candidate_region.get_id(), candidate_region));
    }
    return candidate_regions;
}

std::vector<Interval> Discover::identify_low_coverage_intervals(
    const std::vector<uint32_t>& covg_at_each_position)
{
    std::vector<Interval> identified_regions;
    const auto predicate { [this](uint_least32_t x) {
        return x < this->min_required_covg;
    } };
    auto first { covg_at_each_position.begin() };
    auto previous { covg_at_each_position.begin() };
    auto current { first };

    while (current != covg_at_each_position.end()) {
        previous = current;
        current = std::find_if_not(current, covg_at_each_position.end(), predicate);
        const auto length_of_interval { current - previous };
        const bool interval_between_min_and_max_len
            = (length_of_interval >= this->min_candidate_len)
            and (length_of_interval <= this->max_candidate_len);

        if (interval_between_min_and_max_len) {
            identified_regions.emplace_back(previous - first, current - first);
        }
        if (current == covg_at_each_position.end()) {
            break;
        }
        ++current;
    }
    return identified_regions;
}

void CandidateRegion::add_pileup_entry(
    const std::string& read, const ReadCoordinate& read_coordinate)
{
    const bool read_coord_start_is_past_read_end { read_coordinate.start
        >= read.length() };
    if (not read_coord_start_is_past_read_end) {
        const auto end_pos_of_region_in_read { std::min(
            read_coordinate.end, (uint32_t)read.length()) };
        std::string sequence_in_read_overlapping_region { read.substr(
            read_coordinate.start, end_pos_of_region_in_read - read_coordinate.start) };

        if (!read_coordinate.is_forward) {
            sequence_in_read_overlapping_region
                = reverse_complement(sequence_in_read_overlapping_region);
        }
#ifndef NO_OPENMP
        omp_set_lock(&add_pileup_entry_lock);
        this->pileup.push_back(sequence_in_read_overlapping_region);
        omp_unset_lock(&add_pileup_entry_lock);
#else
        this->pileup.push_back(sequence_in_read_overlapping_region);
#endif
    }
}

void CandidateRegion::write_denovo_paths_to_buffer(CandidateRegionWriteBuffer &buffer,
                                                   const fs::path& temp_dir)
{
    if (denovo_paths.empty()) {
        return;
    }

    const fs::path temp_dir_for_alignment = temp_dir / (filename.string() + "_aln_tmp");
    fs::create_directories(temp_dir_for_alignment);

    {
        const std::string ref_filepath = (temp_dir_for_alignment / "ref.fa").string();
        ofstream ref_filehandler;
        open_file_for_writing(ref_filepath, ref_filehandler);
        ref_filehandler << ">" << buffer.get_sample_name() << "." << name << std::endl;
        ref_filehandler << local_node_max_likelihood_sequence << std::endl;
        ref_filehandler.close();
    }

    {
        const std::string denovo_filepath = (temp_dir_for_alignment / "denovo_paths.fa").string();
        uint32_t denovo_path_index = 0;
        ofstream denovo_filehandler;
        open_file_for_writing(denovo_filepath, denovo_filehandler);
        for (const std::string &denovo_path : denovo_paths) {
            std::string denovo_path_header = string("denovo_path_") + to_string(denovo_path_index);
            denovo_filehandler << ">" << denovo_path_header << std::endl;
            denovo_filehandler << denovo_path << std::endl;
            ++denovo_path_index;
        }
        denovo_filehandler.close();
    }

    std::vector<std::string> variants = get_variants(temp_dir_for_alignment);
    const std::string local_node_max_likelihood_path_as_str =
        LocalNode::to_string_vector(local_node_max_likelihood_path);
    for (const std::string &variant : variants) {
        buffer.add_new_variant(name, local_node_max_likelihood_path_as_str, variant);
    }
}

std::vector<std::string> CandidateRegion::get_variants(const fs::path &temp_dir_for_alignment) const {
    const std::string temp_dir_for_alignment_as_str = temp_dir_for_alignment.string();
    execute_command("cd " + temp_dir_for_alignment_as_str + " && " +
                    "bowtie2-build -q ref.fa ref_index",
                    true, "bowtie2-build failed");
    execute_command("cd " + temp_dir_for_alignment_as_str + " && " +
                    "bowtie2 -x ref_index -S denovo_paths.sam --quiet --end-to-end -f -U denovo_paths.fa",
                    true, "bowtie2 failed");
    execute_command("cd " + temp_dir_for_alignment_as_str + " && " +
                    "samtools view -b denovo_paths.sam -o denovo_paths.bam",
                    true, "samtools view failed");
    execute_command("cd " + temp_dir_for_alignment_as_str + " && " +
                    "bcftools mpileup -f ref.fa denovo_paths.bam -o denovo_paths.mpileup",
                    true, "bcftools mpileup failed");
    execute_command("cd " + temp_dir_for_alignment_as_str + " && " +
                    "bcftools call --ploidy 1 -v -m -f GQ,GP -O v -o denovo_paths.vcf denovo_paths.mpileup",
                    true, "bcftools call failed");

    std::vector<std::string> variants;
    ifstream vcf_filehandler;
    open_file_for_reading((temp_dir_for_alignment/"denovo_paths.vcf").string(), vcf_filehandler);
    std::string line;
    while (getline(vcf_filehandler, line).good()) {
        bool is_comment = line[0]=='#';
        if (is_comment) {
            continue;
        }
        variants.push_back(line);
    }
    vcf_filehandler.close();

    return variants;
}

Fastaq CandidateRegion::generate_fasta_for_denovo_paths()
{
    const bool gzip { false };
    const bool fastq { false };
    Fastaq fasta { gzip, fastq };

    for (size_t i = 0; i < denovo_paths.size(); ++i) {
        const auto& path { denovo_paths.at(i) };
        const std::string path_name { std::string(name).append(".").append(
            std::to_string(i)) };
        fasta.add_entry(path_name, path, "");
    }

    return fasta;
}

TmpPanNode::TmpPanNode(const PanNodePtr& pangraph_node,
    const shared_ptr<LocalPRG>& local_prg,
    const vector<KmerNodePtr>& kmer_node_max_likelihood_path,
    const vector<LocalNodePtr>& local_node_max_likelihood_path,
    const std::string &local_node_max_likelihood_seq)
    : pangraph_node { pangraph_node }
    , local_prg { local_prg }
    , kmer_node_max_likelihood_path { kmer_node_max_likelihood_path }
    , local_node_max_likelihood_path { local_node_max_likelihood_path }
    , local_node_max_likelihood_seq { local_node_max_likelihood_seq }
{
}

PileupConstructionMap Discover::pileup_construction_map(
    CandidateRegions& candidate_regions)
{
    PileupConstructionMap pileup_construction_map;
    for (auto& element : candidate_regions) {
        auto& candidate_region { element.second };
        for (const auto& read_coordinate : candidate_region.read_coordinates) {
            pileup_construction_map[read_coordinate.id].emplace_back(
                &candidate_region, &read_coordinate);
        }
    }
    return pileup_construction_map;
}

void Discover::load_candidate_region_pileups(
    const fs::path& reads_filepath, const CandidateRegions& candidate_regions,
    const PileupConstructionMap& pileup_construction_map, uint32_t threads)
{
    if (candidate_regions.empty() or pileup_construction_map.empty())
        return;

    const uint32_t nb_reads_to_map_in_a_batch = 1000; // nb of reads to map in a batch

    // shared variables - controlled by critical(ReadFileMutex)
    FastaqHandler fh(reads_filepath.string());
    uint32_t id { 0 };

// parallel region
// TODO: this is duplicated code with pangraph_from_read_file(), refactor
#pragma omp parallel num_threads(threads)
    {
        // will hold the reads batch
        std::vector<std::pair<uint32_t, std::string>> sequencesBuffer(
            nb_reads_to_map_in_a_batch, make_pair(0, ""));
        while (true) {
            // read the next batch of reads
            uint32_t nbOfReads = 0;

// read the reads in batch
#pragma omp critical(ReadFileMutex)
            {
                // TODO: we need to read only until the max read id
                for (auto& id_and_sequence : sequencesBuffer) {
                    try {
                        fh.get_next();
                    } catch (std::out_of_range& err) {
                        break;
                    }
                    id_and_sequence = make_pair(id, fh.read);
                    ++nbOfReads;
                    ++id;
                }
            }

            if (nbOfReads == 0)
                break; // we reached the end of the file, nothing else to map

            // process nbOfReads reads
            for (uint32_t i = 0; i < nbOfReads; i++) {
                uint32_t id;
                std::string sequence;
                std::tie(id, sequence) = sequencesBuffer[i];

                auto pileup_construction_map_iterator
                    = pileup_construction_map.find(id);
                const bool is_read_required_for_pileup
                    = pileup_construction_map_iterator != pileup_construction_map.end();

                if (not is_read_required_for_pileup)
                    continue;

                // create all pileups for this read
                for (const auto& pair : pileup_construction_map_iterator->second) {
                    CandidateRegion* candidate_region;
                    const ReadCoordinate* read_coordinate;
                    std::tie(candidate_region, read_coordinate) = pair;

                    candidate_region->add_pileup_entry(sequence, *read_coordinate);
                }
            }
        }
    }
    BOOST_LOG_TRIVIAL(trace) << "Loaded all candidate regions pileups from "
                             << reads_filepath.string();
}

Discover::Discover(uint32_t min_required_covg, uint32_t min_candidate_len,
    uint32_t max_candidate_len, uint16_t candidate_padding, uint32_t merge_dist)
    : min_required_covg { min_required_covg }
    , min_candidate_len { min_candidate_len }
    , max_candidate_len { max_candidate_len }
    , candidate_padding { candidate_padding }
    , merge_dist { merge_dist }
{
}

void CandidateRegionWriteBuffer::add_new_variant(
    const std::string &locus_name,
    const std::string &ML_path,
    const std::string &variant) {
    bool locus_was_not_found = locus_name_to_ML_path.find(locus_name) == locus_name_to_ML_path.end();
    if (locus_was_not_found) {
        locus_name_to_ML_path[locus_name] = ML_path;
    }
    locus_name_to_variants[locus_name].push_back(variant);
}


void CandidateRegionWriteBuffer::write_to_file(const fs::path& output_file) {
    ofstream output_filehandler;
    open_file_for_writing(output_file.string(), output_filehandler);

    output_filehandler << "Sample " << sample_name << std::endl;
    output_filehandler << locus_name_to_ML_path.size() << " loci with denovo variants" << std::endl;
    for (const auto &locus_name_and_ML_path : locus_name_to_ML_path) {
        const auto &locus_name = locus_name_and_ML_path.first;
        output_filehandler << locus_name << std::endl;
        output_filehandler << locus_name_and_ML_path.second << std::endl;
        output_filehandler << locus_name_to_variants[locus_name].size() << " denovo variants for this locus" << std::endl;
        for (const auto &variant : locus_name_to_variants[locus_name]) {
            output_filehandler << variant << std::endl;
        }
    }

    output_filehandler.close();
}