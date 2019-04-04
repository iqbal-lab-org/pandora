#ifndef PANDORA_DENOVO_DISCOVERY_H
#define PANDORA_DENOVO_DISCOVERY_H

#include <stdexcept>
#include <boost/filesystem.hpp>
#include <gatb/gatb_core.hpp>
#include <gatb/debruijn/impl/Simplifications.hpp>
#include "fastaq_handler.h"
#include "denovo_discovery/denovo_utils.h"
#include "denovo_discovery/local_assembly.h"
#include "denovo_discovery/candidate_region.h"


namespace fs = boost::filesystem;


class DenovoDiscovery {
public:
    const uint_least8_t min_covg_for_node_in_assembly_graph { 2 };
    bool clean_assembly_graph { false };

    DenovoDiscovery(const uint_least8_t &kmer_size, const double &read_error_rate);

    void find_paths_through_candidate_region(CandidateRegion &candidate_region);

    double calculate_kmer_coverage(const uint32_t &read_covg, const uint32_t &ref_length) const;

private:
    const uint_least8_t kmer_size;
    const double read_error_rate;

};


#endif //PANDORA_DENOVO_DISCOVERY_H
