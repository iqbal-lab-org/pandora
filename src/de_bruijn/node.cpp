#include "de_bruijn/node.h"

using namespace debruijn;

Node::Node(const uint16_t i, const deque<uint16_t>& n, const uint32_t r) : id(i), hashed_node_ids(n),
                                                                            read_ids({r}) {}

bool Node::operator == (const Node& y) const
{
    if (y.hashed_node_ids.size() != hashed_node_ids.size())
    {
        return false;
    }
    for (uint i=0; i!=hashed_node_ids.size(); ++i)
    {
        if (hashed_node_ids[i] != y.hashed_node_ids[i]
            and hashed_node_ids[i] != y.hashed_node_ids[hashed_node_ids.size()-i])
        {
            return false;
        }
    }
    return true;
}

bool Node::operator!=(const Node &y) const {
    return !(*this == y);
}
