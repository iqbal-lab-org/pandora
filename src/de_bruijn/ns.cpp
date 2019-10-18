#ifndef __DBGNS_CPP_INCLUDED__ // if de_bruijn/ns.cpp hasn't been included yet...
#define __DBGNS_CPP_INCLUDED__

#include <boost/functional/hash.hpp>
#include <iostream>
#include <memory>
#include <unordered_map>

namespace debruijn {
class Node;

class Graph;

typedef std::shared_ptr<debruijn::Node> NodePtr;
typedef std::pair<debruijn::NodePtr, bool> OrientedNodePtr;

template <typename SEQUENCE_OF_GENES> struct seq_hash {
    std::size_t operator()(const SEQUENCE_OF_GENES& seq) const
    {
        std::size_t hash = 0;
        boost::hash_range(hash, seq.begin(), seq.end());
        return hash;
    }
};

template <typename SEQUENCE_OF_GENES, typename T>
using sequence_of_genes_to_data_map
    = std::unordered_map<SEQUENCE_OF_GENES, T, seq_hash<SEQUENCE_OF_GENES>>;

}

#endif
