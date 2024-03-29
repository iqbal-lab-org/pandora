#ifndef PANDORA_MAP_MAIN_H
#define PANDORA_MAP_MAIN_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include "utils.h"
#include "localPRG.h"
#include "localgraph.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "index.h"
#include "estimate_parameters.h"
#include "CLI11.hpp"

using std::set;
using std::vector;

namespace fs = boost::filesystem;

/// Collection of all options of map subcommand.
struct MapOptions {
    fs::path index_file;
    fs::path readsfile;
    fs::path outdir { "pandora" };
    uint32_t threads { 1 };
    fs::path vcf_refs_file;
    uint8_t verbosity { 0 };
    float error_rate { 0.11 };
    uint32_t rng_seed { 0 };
    uint32_t genome_size { 5000000 };
    uint32_t max_diff { 250 };
    float conflicting_clusters_overlap_threshold { 0.8 };
    float conflicting_clusters_minimiser_tolerance { 0.05 };
    bool output_kg { false };
    bool output_vcf { false };
    bool illumina { false };
    bool binomial { false };
    bool do_not_auto_update_params { false };
    uint32_t max_covg { 300 };
    float min_absolute_gene_coverage { 3.0 };
    float min_relative_gene_coverage { 0.05 };
    float max_relative_gene_coverage { 100 };
    float min_gene_coverage_proportion { 0.8 };
    bool no_gene_coverage_filtering { false };
    bool genotype { false };
    bool local_genotype { false };
    bool snps_only { false };
    uint32_t min_cluster_size { 10 };
    uint32_t max_num_kmers_to_avg { 100 };
    uint32_t min_allele_covg_gt { 0 };
    uint32_t min_total_covg_gt { 0 };
    uint32_t min_diff_covg_gt { 0 };
    float min_allele_fraction_covg_gt { 0 };
    float genotyping_error_rate { 0.01 };
    uint16_t confidence_threshold { 1 };
    float partial_matching_lower_bound { 0.5 };
    bool keep_extra_debugging_files { false };
};

void setup_map_subcommand(CLI::App& app);
int pandora_map(MapOptions& opt);

#endif // PANDORA_MAP_MAIN_H
