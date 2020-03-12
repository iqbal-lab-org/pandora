#include "de_bruijn/node.h"
#include "noise_filtering.h"

using namespace debruijn;

debruijn::Node::Node(
    const uint32_t i, const std::deque<uint_least32_t>& n, const uint32_t r)
    : id(i)
    , hashed_node_ids(n)
    , read_ids({ r })
    , out_nodes({})
    , in_nodes({})
{
}

// Nodes are equal if they correspond to the same sequence of oriented pangraph nodes
// either in the forward or reverse complement direction
bool debruijn::Node::operator==(const Node& y) const
{
    if (y.hashed_node_ids.size() != hashed_node_ids.size()) {
        return false;
    }

    bool match = true;
    for (uint32_t i = 0; i < hashed_node_ids.size(); ++i) {
        match = hashed_node_ids[i] == y.hashed_node_ids[i];
        if (!match)
            break;
    }
    if (match)
        return true;

    auto rc = rc_hashed_node_ids(hashed_node_ids);
    for (uint32_t i = 0; i < hashed_node_ids.size(); ++i) {
        match = rc[i] == y.hashed_node_ids[i];
        if (!match)
            break;
    }
    return match;
}

bool debruijn::Node::operator!=(const Node& y) const { return !(*this == y); }

namespace debruijn {
std::ostream& operator<<(std::ostream& out, const Node& m)
{
    std::string sep = "";
    out << "(";
    for (const auto& n : m.hashed_node_ids) {
        out << sep << n;
        sep = ",";
    }
    out << ")";
    return out;
}
}
