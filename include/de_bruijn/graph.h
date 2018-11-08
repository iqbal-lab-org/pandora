#ifndef __DBGRAPH_H_INCLUDED__   // if de_bruijn/graph.h hasn't been included yet...
#define __DBGRAPH_H_INCLUDED__

#include <cstdint>
#include <deque>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>
#include "de_bruijn/ns.cpp"
#include "de_bruijn/node.h"


class debruijn::Graph {
protected:
    uint32_t next_id;
public:
    uint8_t size;
    sequence_of_genes_to_data_map <std::deque<uint_least32_t>, uint32_t> node_hash;
    std::unordered_map<uint32_t, NodePtr> nodes;

    Graph(uint8_t);

    ~Graph();

    OrientedNodePtr add_node(const std::deque<uint_least32_t> &, uint32_t);

    void add_edge(OrientedNodePtr, OrientedNodePtr);

    void remove_node(const uint32_t);

    void remove_read_from_node(const uint32_t, const uint32_t);

    std::unordered_set<uint32_t> get_leaves(uint_least32_t covg_thresh = 1);

    std::unordered_set<uint32_t> get_leaf_tips();

    std::set<std::deque<uint32_t>> get_unitigs();

    void extend_unitig(std::deque<uint32_t> &);

    bool found_in_out_nodes(const NodePtr, const NodePtr) const;

    bool found_in_in_nodes(const NodePtr, const NodePtr) const;

    bool operator==(const Graph &y) const;

    bool operator!=(const Graph &y) const;

    friend std::ostream &operator<<(std::ostream &, const Graph &);
};

#endif
