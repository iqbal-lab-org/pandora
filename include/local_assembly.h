#ifndef PANDORA_LOCAL_ASSEMBLY_H
#define PANDORA_LOCAL_ASSEMBLY_H

#include <iostream>
#include <stack>
#include <unordered_map>
#include <unordered_set>

#include <gatb/gatb_core.hpp>


const long g_max_length{50};
const int g_kmer_size = 5;


using DfsTree = std::unordered_map<std::string, GraphVector<Node>>;
using Paths = std::unordered_set<std::string>;


std::pair<Node, bool> get_node(const std::string &kmer, const Graph &graph);


bool has_ending(std::string const &fullString, std::string const &ending);


DfsTree DFS(const Node &start_node, const Graph &graph);


void get_paths_between(const std::string &start_kmer,
                       const std::string &end_kmer,
                       DfsTree &tree,
                       const Graph &graph,
                       Paths &result);

void get_paths_between_util(const std::string &node,
                            const std::string &end_kmer,
                            std::string acc,
                            const Graph &graph,
                            DfsTree &tree,
                            Paths &full_paths);

#endif //PANDORA_LOCAL_ASSEMBLY_H
