#ifndef PANDORA_DENOVO_DISCOVERY_H
#define PANDORA_DENOVO_DISCOVERY_H

#include "denovo_discovery/gene_interval_info.h"
#include "denovo_discovery/extract_reads.h"
#include "denovo_discovery/local_assembly.h"

namespace fs = boost::filesystem;

namespace denovo_discovery {
    void find_candidates(const std::set <std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
                         const std::string &readfilepath,
                         const fs::path &output_directory,
                         const uint32_t &local_assembly_kmer_size = 0,
                         const uint32_t &kmer_attempts_count = 0,
                         const uint32_t &max_path_length = 0,
                         const uint32_t &padding_size = 0);
}

#endif //PANDORA_DENOVO_DISCOVERY_H
