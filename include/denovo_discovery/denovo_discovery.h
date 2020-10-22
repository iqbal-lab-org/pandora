#ifndef PANDORA_DENOVO_DISCOVERY_H
#define PANDORA_DENOVO_DISCOVERY_H

#include "denovo_discovery/candidate_region.h"
#include "denovo_discovery/denovo_utils.h"
#include "denovo_discovery/local_assembly.h"
#include "fastaq_handler.h"
#include <boost/filesystem.hpp>
#include <gatb/debruijn/impl/Simplifications.hpp>
#include <gatb/gatb_core.hpp>
#include <stdexcept>

namespace fs = boost::filesystem;

class DenovoDiscovery {
public:
    const uint16_t min_covg_for_node_in_assembly_graph;
    bool clean_assembly_graph { false };
    const uint16_t max_insertion_size;

    DenovoDiscovery(const uint16_t& kmer_size, const double& read_error_rate,
        uint16_t max_insertion_size = 15, uint16_t min_covg_for_node_in_assembly_graph = 1);

    void find_paths_through_candidate_region(
        CandidateRegion& candidate_region, const fs::path& denovo_output_directory);

    double calculate_kmer_coverage(
        const uint32_t& read_covg, const uint32_t& ref_length) const;

private:
    const uint16_t kmer_size;
    const double read_error_rate;
};

#endif // PANDORA_DENOVO_DISCOVERY_H
