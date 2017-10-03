#include <iostream>
#include "kmernode.h"
#include "utils.h" // for pointer_values_equal

using namespace std;

KmerNode::KmerNode(uint32_t i, const Path &p) : id(i), path(p), covg({0, 0}),
                                                khash(std::numeric_limits<uint64_t>::max()), num_AT(0) {}

// copy constructor
KmerNode::KmerNode(const KmerNode &other) {
    id = other.id;
    path = other.path;
    khash = other.khash;
    num_AT = other.num_AT;
    covg = {0, 0};
    // NB we don't do edges
}

// Assignment operator
KmerNode &KmerNode::operator=(const KmerNode &other) {
    // check for self-assignment
    if (this == &other)
        return *this;

    id = other.id;
    path = other.path;
    khash = other.khash;
    num_AT = other.num_AT;
    covg = {0, 0};
    // NB we don't do edges

    return *this;
}

std::ostream &operator<<(std::ostream &out, KmerNode const &n) {
    out << n.id << " " << n.path << " " << (unsigned) n.covg[0] << ", " << (unsigned) n.covg[1] << endl;
    for (uint32_t i = 0; i != n.outNodes.size(); ++i) {
        out << n.id << " -> " << n.outNodes[i]->id << endl;
    }
    return out;
}

bool KmerNode::operator==(const KmerNode &y) const {
    return path == y.path;
}
