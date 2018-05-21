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
class DfsNode {
    DfsNode *parent;
    DfsNode *left_child;
    DfsNode *right_sibling;

public:
    DfsNode(std::string kmer);
    std::string sequence;
    DfsNode *const get_parent();
    DfsNode *const get_left_child();
    DfsNode *const get_right_sibling();
    void set_parent(DfsNode *node);
    void set_left_child(DfsNode *node);
    void set_right_sibling(DfsNode *node);
    bool is_leaf();
    bool is_rightmost_child();
    bool is_root();

};


class DfsTree {
    DfsNode *m_root;
public:
    DfsTree();
    DfsTree(DfsNode *root);
    void set_root(DfsNode *root);
    DfsNode* get_root();
    // add function to return a list of paths through tree
};
#endif //PANDORA_LOCAL_ASSEMBLY_H
