#include "denovo_discovery/denovo_discovery.h"


fs::path get_discovered_paths_fname(const GeneIntervalInfo &info,
                                    const uint32_t &local_assembly_kmer_size) {
    return fs::path(
            info.pnode->get_name() + "."
            + std::to_string(info.interval.start) + "-" + std::to_string(info.interval.get_end())
            + "_local_assembly_K" + std::to_string(local_assembly_kmer_size)
            + ".fa");
}

void denovo_discovery::find_candidates(
        const std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
        const std::string &readfilepath, const fs::path &output_directory, const double &error_rate,
        const uint32_t &local_assembly_kmer_size) {

    if (not fs::exists(output_directory)) {
        fs::create_directories(output_directory);
    }

    auto pileups = denovo_discovery::collect_read_pileups(pangraph_coordinate_pairs, readfilepath);
    for (const auto &pileup: pileups) {
        auto &sequences = pileup.second;
        auto &info = pileup.first;
        auto &interval_sequence = info.seq;

        auto fname = get_discovered_paths_fname(info, local_assembly_kmer_size);
        const auto denovo_paths_dir = output_directory / "denovo_paths";
        fs::create_directories(denovo_paths_dir);
        auto discovered_paths_fpath = denovo_paths_dir / fname;

        const uint32_t read_covg = sequences.size();
        const uint32_t ref_length = interval_sequence.length();
        const double expected_kmer_covg = denovo_discovery::calculate_kmer_coverage(read_covg, ref_length,
                                                                                    local_assembly_kmer_size,
                                                                                    error_rate);

        BOOST_LOG_TRIVIAL(debug) << "Running local assembly for: " << info.pnode->get_name() << " - interval ["
                                 << info.interval.start << ", " << info.interval.get_end() << "]";

        local_assembly(sequences, interval_sequence, info.flank_seq_left, info.flank_seq_right, discovered_paths_fpath, local_assembly_kmer_size,
                       expected_kmer_covg, interval_sequence.length() * 2);
    }
}

double
denovo_discovery::calculate_kmer_coverage(const uint32_t &read_covg, const uint32_t &ref_length, const uint32_t &k,
                                          const double &error_rate) {
    if (ref_length == 0) {
        throw std::invalid_argument("ref_length should be greater than 0.");
    } else if (k == 0) {
        throw std::invalid_argument("K should be greater than 0.");
    } else if (error_rate < 0) {
        throw std::invalid_argument("error_rate should not be a negative value.");
    }

    const auto numerator = read_covg * (ref_length - k + 1) * std::pow(1 - error_rate, k);
    const auto denominator = ref_length;
    return numerator / denominator;
}