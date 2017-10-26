#ifndef __DBGRAPH_H_INCLUDED__   // if de_bruijn/graph.h hasn't been included yet...
#define __DBGRAPH_H_INCLUDED__

#include <vector>
#include <memory>
#include <unordered_map>
#include "de_bruijn/ns.cpp"
#include "de_bruijn/node.h"


class debruijn::Graph {
    uint16_t next_id;
    unordered_map<uint16_t, NodePtr> nodes;
public:
    Graph();
    ~Graph();

    NodePtr add_node(const vector<uint16_t>&, uint32_t);
    void add_edge (NodePtr, NodePtr);

    bool operator == (const Graph& y) const;
    bool operator != (const Graph& y) const;

};

#endif