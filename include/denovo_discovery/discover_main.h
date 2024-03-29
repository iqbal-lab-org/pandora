#ifndef PANDORA_DISCOVER_MAIN_H
#define PANDORA_DISCOVER_MAIN_H
#include <boost/filesystem.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include "CLI11.hpp"
#include "denovo_discovery/denovo_utils.h"
#include "utils.h"
#include "index.h"
#include "pangenome/pangraph.h"
#include "estimate_parameters.h"

namespace fs = boost::filesystem;

/// Collection of all options of discover subcommand.
struct DiscoverOptions {
    fs::path index_file;
    fs::path reads_idx_file;
    fs::path outdir { "pandora_discover" };
    uint32_t threads { 1 };
    uint8_t verbosity { 0 };
    float error_rate { 0.11 };
    uint32_t rng_seed { 0 };
    uint32_t genome_size { 5000000 };
    uint32_t max_diff { 250 };
    float conflicting_clusters_overlap_threshold { 0.8 };
    float conflicting_clusters_minimiser_tolerance { 0.05 };
    bool output_kg { false };
    bool illumina { false };
    bool binomial { false };
    bool do_not_auto_update_params { false };
    uint32_t max_covg { 600 };
    float min_absolute_gene_coverage { 3.0 };
    float min_relative_gene_coverage { 0.05 };
    float max_relative_gene_coverage { 100 };
    float min_gene_coverage_proportion { 0.8 };
    bool no_gene_coverage_filtering { false };
    uint32_t min_cluster_size { 10 };
    uint32_t max_num_kmers_to_avg { 100 };
    float partial_matching_lower_bound { 0.5 };
    bool keep_extra_debugging_files { false };
};

void setup_discover_subcommand(CLI::App& app);
int pandora_discover(DiscoverOptions& opt);

#endif // PANDORA_DISCOVER_MAIN_H
