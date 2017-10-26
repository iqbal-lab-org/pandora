#include <vector>
#include <fstream>
#include <unordered_set>
#include <iostream>
#include <memory>
#include "de_bruijn/graph.h"
#include "de_bruijn/node.h"
#include "utils.h"

using namespace debruijn;

Graph::Graph() : next_id(0) {
    nodes.reserve(20000);
};


Graph::~Graph()
{
    nodes.clear();
}

NodePtr Graph::add_node (const vector<uint16_t>& node_ids, uint32_t read_id)
{
    NodePtr n (make_shared<Node>(next_id, node_ids, read_id));
    for (auto c : nodes)
    {
	    if (*c.second == *n)
	    {
	        c.second->read_ids.insert(read_id);
	        return c.second;
	    }
    }
    cout << "new tuple " << next_id << endl;
    nodes[next_id] = n;
    next_id++;
    return n;
}

void Graph::add_edge (NodePtr from, NodePtr to)
{
    if (find(from->out_nodes.begin(), from->out_nodes.end(), to->id) == from->out_nodes.end())
    {
        from->out_nodes.push_back(to->id);
        to->out_nodes.push_back(from->id);
    }
}

bool Graph::operator == (const Graph& y) const
{
    // want the graphs to have the same nodes, even if
    // the ids given them is different.
    if (nodes.size() != y.nodes.size())
    {
	    return false;
    }
    for (const auto t : nodes) {
        bool found = false;
        for (const auto s : y.nodes) {
            if (*t.second == *s.second) {
                found = true;
                break;
            }
        }
        if (found == false) {
            return false;
        }
    }
    // nodes can't be equal within a graph so don't need to check vice versa
    return true;
}

bool Graph::operator!=(const Graph &y) const {
    return !(*this == y);
}
