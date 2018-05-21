//
// Created by Michael Benjamin Hall on 11/05/2018.
//

#ifndef PANDORA_LOCAL_ASSEMBLY_H
#define PANDORA_LOCAL_ASSEMBLY_H
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <stack>
#include <unordered_map>
#include <set>


bool get_node(Node &node, Graph &graph);
std::unordered_map<std::string, GraphVector<Node>>& DFS(Node &start_node, Graph &graph);
void print_path(std::unordered_map<std::string, GraphVector<Node>> &tree, const std::string start_node,
                Graph &graph, std::vector<std::string> &result);
void helper(std::string node, std::string acc, Graph &graph, std::unordered_map<std::string, GraphVector<Node>> &tree,
            std::vector<std::string> &result);

#endif //PANDORA_LOCAL_ASSEMBLY_H
