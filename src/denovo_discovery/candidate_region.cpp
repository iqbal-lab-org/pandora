#include "denovo_discovery/candidate_region.h"


size_t std::hash<CandidateRegionIdentifier>::operator()(const CandidateRegionIdentifier &id) const {
    return std::hash<int>()(std::get<0>(id).start) ^ std::hash<int>()(std::get<0>(id).get_end()) ^
           std::hash<std::string>()(std::get<1>(id));
}


CandidateRegion::CandidateRegion(const Interval &interval, std::string name)
        : interval { interval }, name { std::move(name) }, interval_padding { 0 } {
    initialise_filename();
}


CandidateRegion::CandidateRegion(const Interval &interval, std::string name, const uint_least16_t &interval_padding)
        : interval { interval }, name { std::move(name) }, interval_padding { interval_padding } {
    initialise_filename();
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
    //    const auto candidate_intervals { identify_low_coverage_intervals(covgs_along_localnode_path) };
    const auto candidate_intervals { find_coverage_anomalies(covgs_along_localnode_path, 22) };

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


std::vector<Interval> find_coverage_anomalies(const std::vector<uint32_t> &per_base_coverage,
                                              const uint_least32_t min_dist_between_candidates) {
    const auto log_change_in_covg { transform_to_log_change(per_base_coverage) };
    const double covg_threshold { 1.0 };
    std::vector<uint_least32_t> candidate_idxs;

    for (uint_least32_t i = 0; i < log_change_in_covg.size(); i++) {
        if (std::abs(log_change_in_covg[i]) >= covg_threshold) {
            candidate_idxs.push_back(i + 1); // add one as log_change_in_covg has length 1 less than per_base_coverage
        }
    }

    return collapse_candidate_indicies_into_intervals(candidate_idxs, min_dist_between_candidates);
}


std::vector<Interval> collapse_candidate_indicies_into_intervals(const std::vector<uint_least32_t> &candidate_idxs,
                                                                 const uint_least32_t min_dist_between_candidates) {
    std::vector<Interval> identified_regions;

    if (candidate_idxs.empty()) {
        return identified_regions;
    }
    auto start { candidate_idxs.at(0) };
    for (size_t i = 0; i < candidate_idxs.size() - 1; i++) {
        const bool close_to_next { (candidate_idxs[i + 1] - candidate_idxs[i]) <= min_dist_between_candidates };
        if (not close_to_next) {
            identified_regions.emplace_back(start, candidate_idxs[i] + 1);
            start = candidate_idxs[i + 1];
        }
    }
    identified_regions.emplace_back(start, candidate_idxs.back() + 1);

    return identified_regions;
}


std::vector<double> transform_to_log_change(const std::vector<uint32_t> &per_base_coverage) {
    if (per_base_coverage.size() < 2) {
        return std::vector<double>();
    }

    const auto log_division_of { [](uint32_t val, uint32_t prev) { return std::log((val + 1.0) / (prev + 1.0)); }};
    std::vector<double> result(per_base_coverage.size());

    std::adjacent_difference(per_base_coverage.begin(), per_base_coverage.end(), result.begin(), log_division_of);
    result.erase(result.begin());
    return result;
}


void CandidateRegion::generate_read_pileup(const fs::path &reads_filepath) {
    FastaqHandler readfile(reads_filepath.string());
    if (readfile.eof()) return;

    uint32_t last_id { 0 };

    for (auto &read_coordinate : this->read_coordinates) {
        assert(last_id <= read_coordinate.id);  // stops from iterating through fastq multiple times
        readfile.get_id(read_coordinate.id);

        const bool read_coord_start_is_past_read_end { read_coordinate.start >= readfile.read.length() };
        if (read_coord_start_is_past_read_end) continue;

        const auto end_pos_of_region_in_read { std::min(read_coordinate.end, (uint32_t) readfile.read.length()) };
        std::string sequence_in_read_overlapping_region {
                readfile.read.substr(read_coordinate.start, end_pos_of_region_in_read - read_coordinate.start) };

        if (!read_coordinate.is_forward) {
            sequence_in_read_overlapping_region = reverse_complement(sequence_in_read_overlapping_region);
        }
        this->pileup.push_back(sequence_in_read_overlapping_region);
        last_id = read_coordinate.id;
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
