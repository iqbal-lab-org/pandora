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
#include <boost/filesystem.hpp>

using DfsTree = std::unordered_map<std::string, GraphVector<Node>>;
using Paths = std::vector<std::string>;

namespace logging = boost::log;
namespace fs = boost::filesystem;

constexpr float g_covg_scaling_factor{0.1};
constexpr auto g_max_num_paths{50};

std::pair<Node, bool> get_node(const std::string &kmer, const Graph &graph);

bool has_ending(std::string const &fullString, std::string const &ending);

DfsTree DFS(const Node &start_node, const Graph &graph);

Paths get_paths_between(const std::string &start_kmer,
                        const std::string &end_kmer,
                        std::unordered_map<string, GraphVector<Node> > &tree,
                        const Graph &graph,
                        const uint32_t &max_path_length,
                        const double &expected_coverage = 1);

void get_paths_between_util(const std::string &node,
                            const std::string &end_kmer,
                            std::string path_accumulator,
                            const Graph &graph,
                            std::unordered_map<string, GraphVector<Node>> &tree,
                            Paths &full_paths,
                            const uint32_t &max_path_length,
                            const double &expected_kmer_covg,
                            const float &covg_scaling_factor = g_covg_scaling_factor,
                            uint32_t kmers_below_threshold = 0);

void write_paths_to_fasta(const boost::filesystem::path &filepath,
                          const Paths &paths,
                          const uint32_t &line_width = 80);

void local_assembly(const std::vector<std::string> &sequences, const std::string &slice_sequence,
                    const std::string &flank_left, const std::string &flank_right, const fs::path &out_path,
                    const uint32_t &kmer_size, const uint32_t &max_path_length,
                    const double &expected_coverage = 1, const bool &clean_graph = false,
                    const uint32_t &min_coverage = 2);

void do_graph_clean(Graph &graph, const uint16_t &num_cores = 1);

std::string reverse_complement(const std::string &forward);

void remove_graph_file();

std::vector<std::string> generate_start_kmers(const std::string &sequence,
                                              const uint32_t &k,
                                              uint32_t num_to_generate);

std::vector<std::string> generate_end_kmers(const std::string &sequence,
                                            const uint32_t &k,
                                            uint32_t num_to_generate);

#endif //PANDORA_LOCAL_ASSEMBLY_H
