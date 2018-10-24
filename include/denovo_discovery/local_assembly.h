#ifndef PANDORA_LOCAL_ASSEMBLY_H
#define PANDORA_LOCAL_ASSEMBLY_H

#include <iostream>
#include <fstream>
#include <stack>
#include <unordered_map>
#include <unordered_set>

#include <gatb/gatb_core.hpp>
#include <gatb/debruijn/impl/Simplifications.hpp>
#include <sys/stat.h>


#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem/path.hpp>


using DfsTree = std::unordered_map<std::string, GraphVector<Node>>;
using Paths = std::vector<std::string>;

namespace logging = boost::log;
namespace fs = boost::filesystem;

const uint32_t g_max_length{300};
const uint32_t g_local_assembly_kmer_size{11};
const auto g_log_level{logging::trivial::debug};
const float COVG_SCALING_FACTOR{0.2};
const uint32_t g_kmer_attempts_count{10};

std::pair<Node, bool> get_node(const std::string &kmer, const Graph &graph);


bool has_ending(std::string const &fullString, std::string const &ending);


DfsTree DFS(const Node &start_node, const Graph &graph);


Paths get_paths_between(const std::string &start_kmer, const std::string &end_kmer,
                        std::unordered_map<string, GraphVector<Node>> &tree, const Graph &graph,
                        const unsigned long max_path_length, const double &expected_coverage = 1);

void get_paths_between_util(const std::string &node, const std::string &end_kmer, std::string path_accumulator,
                            const Graph &graph, std::unordered_map<string, GraphVector<Node>> &tree, Paths &full_paths,
                            const unsigned long &max_path_length = g_max_length, const double &expected_kmer_covg = 1,
                            unsigned int kmers_below_threshold = 0);


void write_paths_to_fasta(const boost::filesystem::path &filepath,
                          const Paths &paths,
                          unsigned long line_width = 80);

void local_assembly(const std::vector<std::string> &sequences,
                    std::unordered_set<std::string> &start_kmers,
                    std::unordered_set<std::string> &end_kmers,
                    const fs::path &out_path,
                    const unsigned int kmer_size,
                    const unsigned long max_path_length,
                    const double &expected_coverage = 1,
                    const bool clean_graph = false,
                    const unsigned int min_coverage = 2);

void do_graph_clean(Graph &graph, const int num_cores = 1);

std::string reverse_complement(const std::string &forward);

bool file_exists(const std::string &name);

void remove_graph_file(const std::string &filepath = "");

std::unordered_set<std::string> generate_start_kmers(const std::string &sequence, const unsigned int k, unsigned int n);

std::unordered_set<std::string> generate_end_kmers(const std::string &sequence, const unsigned int k, unsigned int n);

#endif //PANDORA_LOCAL_ASSEMBLY_H
