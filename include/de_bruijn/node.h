#ifndef __DBNODE_H_INCLUDED__   // if de_bruijn/node.h hasn't been included yet...
#define __DBNODE_H_INCLUDED__

#include <vector>
#include <unordered_set>
#include <memory>
#include "de_bruijn/ns.cpp"

typedef std::shared_ptr<debruijn::Node> NodePtr;

class debruijn::Node {
    uint16_t id;
    std::vector<uint16_t> hashed_node_ids;
    std::unordered_multiset<uint32_t> read_ids;
public:
    std::unordered_set<uint16_t> out_nodes;

    Node(const uint16_t, const vector<uint16_t>&, const uint32_t);

    bool operator == (const Node& y) const;
    bool operator != (const Node& y) const;

    friend class Graph;

};

#endif