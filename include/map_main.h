#ifndef PANDORA_MAP_MAIN_H
#define PANDORA_MAP_MAIN_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <cassert>

#include "utils.h"
#include "localPRG.h"
#include "localgraph.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "index.h"
#include "estimate_parameters.h"
#include "noise_filtering.h"

#include "denovo_discovery/denovo_utils.h"
#include "denovo_discovery/denovo_discovery.h"
#include "CLI11.hpp"

using std::set;
using std::vector;

namespace fs = boost::filesystem;

/// Collection of all options of map subcommand.
struct MapOptions {
    std::string prgfile;
    std::string readsfile;
    std::string outdir { "pandora" };
    uint32_t window_size { 14 };
    uint32_t kmer_size { 15 };
    uint32_t threads { 1 };
    std::string vcf_refs_file;
    uint8_t verbosity { 0 };
    float error_rate { 0.11 };
    uint32_t genome_size { 5000000 };
    uint32_t max_diff { 250 };
    bool output_kg { false };
    bool output_vcf { false };
    bool output_comparison_paths { false };
    bool output_covgs { false };
    bool output_mapped_read_fa { false };
    bool illumina { false };
    bool clean { false };
    bool binomial { false };
    uint32_t max_covg { 300 };
    std::string genotype;
    bool snps_only { false };
    bool discover { false };
    uint32_t denovo_kmer_size { 11 };
    uint32_t min_cluster_size { 10 };
    uint32_t max_num_kmers_to_avg { 100 };
    uint32_t min_allele_covg_gt { 0 };
    uint32_t min_total_covg_gt { 0 };
    uint32_t min_diff_covg_gt { 0 };
    float min_allele_fraction_covg_gt { 0 };
    float genotyping_error_rate { 0.01 };
    uint16_t confidence_threshold { 1 };
};

void setup_map_subcommand(CLI::App& app);
int pandora_map(MapOptions& opt);

#endif // PANDORA_MAP_MAIN_H
