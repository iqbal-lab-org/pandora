#include "denovo_discovery/denovo_discovery.h"


fs::path get_discovered_paths_fname(const GeneIntervalInfo &info,
                                    const uint32_t &local_assembly_kmer_size) {
    return fs::path(
            info.pnode->get_name() + "."
            + std::to_string(info.interval.start) + "-" + to_string(info.interval.get_end())
            + "_local_assembly_K" + std::to_string(local_assembly_kmer_size)
            + ".fa");
}

void denovo_discovery::find_candidates(
        const std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
        const std::string &readfilepath, const fs::path &output_directory, const double &error_rate,
        const uint32_t &local_assembly_kmer_size, const uint32_t &kmer_attempts_count) {
    logging::core::get()->set_filter(logging::trivial::severity >= g_log_level);

    if (not fs::exists(output_directory)) {
        fs::create_directories(output_directory);
    }

    auto pileups = denovo_discovery::collect_read_pileups(pangraph_coordinate_pairs, readfilepath);
    for (const auto &pileup: pileups) {
        auto &sequences = pileup.second;
        auto &info = pileup.first;
        auto &interval_sequence = info.seq;

        auto start_kmers{
                generate_start_kmers(interval_sequence,
                                     local_assembly_kmer_size,
                                     kmer_attempts_count)
        };
        auto end_kmers{
                generate_end_kmers(interval_sequence,
                                   local_assembly_kmer_size,
                                   kmer_attempts_count)
        };
        // todo: create directory to put paths into i.e denovo/
        auto fname = get_discovered_paths_fname(info, local_assembly_kmer_size);
        auto discovered_paths_fpath = output_directory / fname;

        const uint32_t read_covg = sequences.size();
        const uint32_t ref_length = interval_sequence.length();
        const double expected_kmer_covg = denovo_discovery::calculate_kmer_coverage(read_covg, ref_length,
                                                                                    local_assembly_kmer_size,
                                                                                    error_rate);

        BOOST_LOG_TRIVIAL(debug) << "Running local assembly for: " << info.pnode->get_name() << " - interval ["
                                 << info.interval.start << ", " << info.interval.get_end() << "]";

        local_assembly(sequences, start_kmers, end_kmers, discovered_paths_fpath, local_assembly_kmer_size,
                       expected_kmer_covg);
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