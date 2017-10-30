#ifndef __PANREAD_H_INCLUDED__   // if panread.h hasn't been included yet...
#define __PANREAD_H_INCLUDED__

#include <vector>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include "minihits.h"
#include "pangenome/ns.cpp"

class pangenome::Read {
    const uint32_t id; // corresponding the the read id
    vector<NodePtr> nodes;
    vector<bool> node_orientations;
    std::unordered_map<uint32_t,std::set<MinimizerHitPtr, pComp_path>> hits; // from node id to cluster of hits against that node in this read
public:
    Read(const uint32_t);
    void add_hits(const uint32_t, const std::set<MinimizerHitPtr, pComp>&);
    void remove_node(NodePtr);
    // remove nodes
    //void replace_node(NodePtr, NodePtr);
    // replace nodes

    bool operator == (const Read& y) const;
    bool operator != (const Read& y) const;
    bool operator < (const Read& y) const;

    friend std::ostream& operator<< (std::ostream& out, const Read& r);
    friend class pangenome::Graph;
    friend class pangenome::Node;
};

#endif
