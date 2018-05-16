//
// Created by Michael Benjamin Hall on 11/05/2018.
//

#ifndef PANDORA_LOCAL_ASSEMBLY_H
#define PANDORA_LOCAL_ASSEMBLY_H
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <stack>


bool kmer_in_graph(const char *kmer, Graph &graph);

class DfsNode {
    DfsNode *parent;
    DfsNode *left_child;
    DfsNode *right_sibling;

public:
    DfsNode();
    DfsNode *const get_parent();
    DfsNode *const get_left_child();
    DfsNode *const get_right_sibling();
    void set_parent(DfsNode *node);
    void set_left_child(DfsNode *node);
    void set_right_sibling(DfsNode *node);
    bool has_children();
    bool is_rightmost_child();
    bool is_root();

};


class DfsTree {
    DfsNode *root;
public:
    DfsTree();
    // add function to return a list of paths through tree
};
#endif //PANDORA_LOCAL_ASSEMBLY_H
