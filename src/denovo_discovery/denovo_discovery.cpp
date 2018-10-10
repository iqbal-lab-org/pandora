#include <boost/filesystem/operations.hpp>
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
        const std::string &readfilepath,
        const fs::path &output_directory,
        const uint32_t &local_assembly_kmer_size,
        const uint32_t &kmer_attempts_count,
        const uint32_t &max_path_length,
        const uint32_t &padding_size) {

    if (not fs::exists(output_directory))
        fs::create_directories(output_directory);

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

        /*
        local_assembly(sequences,
                       start_kmers,
                       end_kmers,
                       discovered_paths_fpath,
                       local_assembly_kmer_size,
                       max_path_length);
        */
    }
}