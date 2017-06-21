#include <iostream>
#include <string>
#include <algorithm>
#include <limits>
#include <cassert>
#include "kmernode.h"
#include "path.h"
#include "utils.h" // for pointer_values_equal

using namespace std;

KmerNode::KmerNode (uint32_t i, const Path& p): id(i), path(p), covg({0,0}), khash(std::numeric_limits<uint64_t>::max()), num_AT(0) {}

std::ostream& operator<< (std::ostream & out, KmerNode const& n) {
    out << n.id << " " << n.path << " " << (unsigned)n.num_AT << endl;
    for (uint32_t i=0; i!=n.outNodes.size(); ++i)
    {
        out << n.id << " -> " << n.outNodes[i]->id << endl;
    }
    return out ;
}

bool KmerNode::operator == (const KmerNode& y) const {
    if (!(path == y.path)) {return false;}
    return true;
}
