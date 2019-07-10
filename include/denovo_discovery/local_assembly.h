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


namespace logging = boost::log;
namespace fs = boost::filesystem;
using DfsTree = std::unordered_map<std::string, GraphVector<Node>>;
using DenovoPaths = std::vector<std::string>;

constexpr float COVG_SCALING_FACTOR { 0.1 };
constexpr auto MAX_NUMBER_CANDIDATE_PATHS { 25 };


class LocalAssemblyGraph : public Graph {
public:
    LocalAssemblyGraph& operator=(const Graph &graph);

    std::pair<Node, bool> get_node(const std::string &query_kmer);

    DenovoPaths get_paths_between(const Node &start_node, const Node &end_node,
                                  const uint32_t &max_path_length,
                                  const double &expected_coverage = 1);

private:
    DfsTree depth_first_search_from(const Node &start_node, bool reverse=false);

    void build_paths_between(const std::string &start_kmer, const std::string &end_kmer, std::string path_accumulator,
                             DfsTree &tree, DfsTree &reverse_tree, DenovoPaths &paths_between_queries,
                             const uint32_t &max_path_length, const double &expected_kmer_covg,
                             const float &required_percent_of_expected_covg = COVG_SCALING_FACTOR,
                             uint32_t num_kmers_below_threshold = 0);
};


void clean(Graph &graph, const uint16_t &num_cores = 1);

bool string_ends_with(std::string const &query, std::string const &ending);

std::string reverse_complement(const std::string &forward);

void remove_graph_file();

std::vector<std::string> generate_start_kmers(const std::string &sequence, const uint16_t &k, uint32_t num_to_generate);

std::vector<std::string> generate_end_kmers(const std::string &sequence, const uint32_t &k, uint32_t num_to_generate);

std::vector<std::string> all_kmers_in(const std::string &query, const uint_least8_t k_size);

#endif //PANDORA_LOCAL_ASSEMBLY_H
