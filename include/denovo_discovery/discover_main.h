#ifndef PANDORA_DISCOVER_MAIN_H
#define PANDORA_DISCOVER_MAIN_H
#include <boost/filesystem.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include "CLI11.hpp"
#include "utils.h"
#include "index.h"
#include "pangenome/pangraph.h"
#include "noise_filtering.h"
#include "estimate_parameters.h"
#include "denovo_discovery/candidate_region.h"
#include "denovo_discovery/denovo_discovery.h"

constexpr auto MAX_DENOVO_K { 32 };
namespace fs = boost::filesystem;

/// Collection of all options of discover subcommand.
struct DiscoverOptions {
    fs::path prgfile;
    fs::path reads_idx_file;
    fs::path outdir { "pandora_discover" };
    uint32_t window_size { 14 };
    uint32_t kmer_size { 15 };
    uint32_t threads { 1 };
    uint8_t verbosity { 0 };
    float error_rate { 0.11 };
    uint32_t genome_size { 5000000 };
    uint32_t max_diff { 250 };
    bool output_kg { false };
    bool output_mapped_read_fa { false };
    bool illumina { false };
    bool clean { false };
    bool binomial { false };
    uint32_t max_covg { 600 };
    uint16_t denovo_kmer_size { 15 };
    uint16_t max_insertion_size { 15 };
    uint16_t min_covg_for_node_in_assembly_graph { 2 };
    uint32_t min_candidate_covg { 3 };
    uint32_t min_candidate_len { 1 };
    uint32_t max_candidate_len { 30 };
    int max_num_candidate_paths { 25 };
    uint32_t merge_dist { 15 };
    uint32_t min_cluster_size { 10 };
    uint32_t max_num_kmers_to_avg { 100 };
    bool clean_dbg { false };
};

void setup_discover_subcommand(CLI::App& app);
int pandora_discover(DiscoverOptions& opt);

#endif // PANDORA_DISCOVER_MAIN_H
