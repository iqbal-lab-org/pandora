#include <boost/filesystem/operations.hpp>

#include "denovo_discovery/denovo_discovery.h"
#include "denovo_discovery/local_assembly.h"


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
        const std::string &readfilepath, const fs::path &output_directory,
        const double &error_rate, const uint32_t &local_assembly_kmer_size,
        const uint32_t &kmer_attempts_count, const uint32_t &padding_size) {

    if (not fs::exists(output_directory)) {
        fs::create_directories(output_directory);
    }

    auto pileups = collect_read_pileups(pangraph_coordinate_pairs, readfilepath);
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

        auto fname = get_discovered_paths_fname(info, local_assembly_kmer_size);
        auto discovered_paths_fpath = output_directory / fname;

        //todo: revisit expected maximum path length with zam
        const uint32_t expected_max_path_len = interval_sequence.length() * 2;
        const uint32_t max_path_length = (expected_max_path_len > g_max_length) ? g_max_length : expected_max_path_len;

        //todo: add expected coverage calculation
        // calculate coverage for the slice
        //const auto slice_coverage{fa.calculate_kmer_coverage(sub_lmp_as_string.length(), g_local_assembly_kmer_size)};
        /*This is how I was doing kmer coverage before
         * // returns coverage as just the number of reads in the fastaq
double Fastaq::calculate_coverage() const {
    return sequences.size();
}

// calculates coverage as number of bases / length of a given reference
double
Fastaq::calculate_kmer_coverage(const unsigned long &ref_length, const unsigned int k, const double &error_rate) const {
    const auto D{this->calculate_coverage()};

    return (D * (ref_length - k + 1)) / (ref_length * pow(1-error_rate, k));
}


         *
         */
        const auto read_covg = sequences.size();
        const auto ref_length = interval_sequence.length();
        const auto expected_kmer_covg = denovo_discovery::calculate_kmer_coverage(read_covg, ref_length,
                                                                                  local_assembly_kmer_size,
                                                                                  error_rate);

        local_assembly(sequences,
                       start_kmers,
                       end_kmers,
                       discovered_paths_fpath,
                       local_assembly_kmer_size,
                       max_path_length,
                       expected_kmer_covg);
    }
}

double
denovo_discovery::calculate_kmer_coverage(const uint32_t &read_covg, const uint32_t &ref_length, const uint32_t k,
                                          const double &error_rate) {
    const auto numerator = read_covg * (ref_length - k + 1);
    const auto denominator = ref_length * std::pow(1 - error_rate, k);
    return numerator / denominator;
}