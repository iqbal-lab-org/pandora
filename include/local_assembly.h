#ifndef PANDORA_LOCAL_ASSEMBLY_H
#define PANDORA_LOCAL_ASSEMBLY_H

#include <iostream>
#include <stack>
#include <unordered_map>
#include <set>

#include <gatb/gatb_core.hpp>


std::pair<Node, bool> get_node(const std::string &kmer, Graph &graph);

std::unordered_map<std::string, GraphVector<Node>> &DFS(Node &start_node, Graph &graph);

void print_path(std::unordered_map<std::string, GraphVector<Node>> &tree,
                const std::string start_node,
                Graph &graph,
                std::vector<std::string> &result);

void helper(std::string node,
            std::string acc,
            Graph &graph,
            std::unordered_map<std::string, GraphVector<Node>> &tree,
            std::vector<std::string> &result);

#endif //PANDORA_LOCAL_ASSEMBLY_H
