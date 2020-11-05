#ifndef PANDORA_COMPARE_MAIN_H
#define PANDORA_COMPARE_MAIN_H
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <tuple>
#include <functional>
#include <cctype>
#include <fstream>
#include <algorithm>
#include <map>
#include <cassert>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include "utils.h"
#include "localPRG.h"
#include "localgraph.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "index.h"
#include "noise_filtering.h"
#include "estimate_parameters.h"
#include "OptionsAggregator.h"
#include "CLI11.hpp"

using std::set;
using std::vector;
using SampleIdText = std::string;
using SampleFpath = std::string;

struct CompareOptions {
    fs::path prgfile;
    fs::path reads_idx_file;
    fs::path outdir { "pandora" };
    uint32_t window_size { 14 };
    uint32_t kmer_size { 15 };
    uint32_t threads { 1 };
    fs::path vcf_refs_file;
    uint8_t verbosity { 0 };
    float error_rate { 0.11 };
    uint32_t genome_size { 5000000 };
    uint32_t max_diff { 250 };
    bool output_vcf { false };
    bool illumina { false };
    bool clean { false };
    bool binomial { false };
    uint32_t max_covg { 300 };
    bool genotype { false };
    bool local_genotype { false };
    uint32_t min_cluster_size { 10 };
    uint32_t max_num_kmers_to_avg { 100 };
    uint32_t min_allele_covg_gt { 0 };
    uint32_t min_total_covg_gt { 0 };
    uint32_t min_diff_covg_gt { 0 };
    float min_allele_fraction_covg_gt { 0 };
    float genotyping_error_rate { 0.01 };
    uint16_t confidence_threshold { 1 };
};

std::vector<std::pair<SampleIdText, SampleFpath>> load_read_index(
    const fs::path& read_index_fpath);
void setup_compare_subcommand(CLI::App& app);
int pandora_compare(CompareOptions& opt);

#endif // PANDORA_COMPARE_MAIN_H
