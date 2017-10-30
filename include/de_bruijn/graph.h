#ifndef __DBGRAPH_H_INCLUDED__   // if de_bruijn/graph.h hasn't been included yet...
#define __DBGRAPH_H_INCLUDED__

#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <deque>
#include "de_bruijn/ns.cpp"
#include "de_bruijn/node.h"

class debruijn::Graph {
    uint16_t next_id;
    uint8_t size;
    unordered_map<uint16_t, NodePtr> nodes;
public:
    Graph(uint8_t);
    ~Graph();

    NodePtr add_node(const vector<uint16_t>&, uint32_t);
    void add_edge (NodePtr, NodePtr);
    void remove_node(const uint16_t);

    unordered_set<uint16_t> get_leaves();
    set<deque<uint16_t>> get_unitigs();
    void extend_unitig(deque<uint16_t>&);

    bool operator == (const Graph& y) const;
    bool operator != (const Graph& y) const;

};

#endif