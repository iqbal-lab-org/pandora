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
    uint fwd_count = 0, rev_count = 0;
    bool fwd_match = false, rev_match = false;
    for (uint i=0; i!=hashed_node_ids.size(); ++i)
    {
        fwd_match = (hashed_node_ids[i] == y.hashed_node_ids[i]);
        rev_match = (hashed_node_ids[i] == y.hashed_node_ids[y.hashed_node_ids.size()-1-i]
                                           + 1*(y.hashed_node_ids[y.hashed_node_ids.size()-1-i]%2 == 0)
                                           - 1*(y.hashed_node_ids[y.hashed_node_ids.size()-1-i]%2 == 1));
        if (fwd_match)
        {
            fwd_count++;
        }
        if (rev_match)
        {
            rev_count++;
        }
        if (!fwd_match and !rev_match)
        {
            return false;
        }
    }
    if (fwd_count != hashed_node_ids.size() and rev_count != hashed_node_ids.size())
    {
        return false;
    }
    return true;
}

bool Node::operator!=(const Node &y) const {
    return !(*this == y);
}
