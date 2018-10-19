#ifndef PANDORA_DENOVO_DISCOVERY_H
#define PANDORA_DENOVO_DISCOVERY_H

#include "denovo_discovery/gene_interval_info.h"
#include "denovo_discovery/extract_reads.h"
#include "denovo_discovery/local_assembly.h"

namespace fs = boost::filesystem;

namespace denovo_discovery {
    void find_candidates(
            const std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
            const std::string &readfilepath,
            const fs::path &output_directory,
            const double &error_rate,
            const uint32_t &local_assembly_kmer_size = g_local_assembly_kmer_size,
            const uint32_t &kmer_attempts_count = g_kmer_attempts_count,
            const uint32_t &padding_size = 0);

    double calculate_kmer_coverage(const uint32_t &read_covg, const uint32_t &ref_length, const uint32_t k,
                                   const double &error_rate);
}

#endif //PANDORA_DENOVO_DISCOVERY_H
