#include <vector>
#include <fstream>
#include <unordered_set>
#include <set>
#include <deque>
#include <iostream>
#include <memory>
#include <cassert>
#include "de_bruijn/graph.h"
#include "de_bruijn/node.h"
#include "utils.h"

using namespace debruijn;

Graph::Graph(uint8_t s) : next_id(0), size(s) {
    nodes.reserve(20000);
};


Graph::~Graph()
{
    nodes.clear();
}

// add a node in dbg corresponding to a fixed size vector of pangenome graph
// node/orientation ids and labelled with the read_ids which cover it
NodePtr Graph::add_node (const vector<uint16_t>& node_ids, uint32_t read_id)
{
    assert(node_ids.size() == size);
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

// add an undirected edge
void Graph::add_edge (NodePtr from, NodePtr to)
{
    if (from->out_nodes.find(to->id) == from->out_nodes.end())
    {
        from->out_nodes.insert(to->id);
        to->out_nodes.insert(from->id);
    }
}

// remove de bruijn node with id given
void Graph::remove_node(const uint16_t dbg_node_id)
{
    auto it = nodes.find(dbg_node_id);
    if ( it != nodes.end())
    {
        // remove this node from lists of out nodes from other graph nodes
        for (auto n : it->second->out_nodes)
        {
            nodes[n]->out_nodes.erase(dbg_node_id);
        }

        // and remove from nodes
        nodes.erase(dbg_node_id);
    }
}

// get the dbg node ids corresponding to leaves
unordered_set<uint16_t> Graph::get_leaves()
{
    unordered_set<uint16_t> s;
    for (auto c : nodes)
    {
        if (c.second->out_nodes.size() <= 1)
        {
            s.insert(c.second->id);
        }
    }
    return s;
}

// get deques of dbg node ids corresponding to maximal non-branching paths in dbg
set<deque<uint16_t>> Graph::get_unitigs()
{
    set<deque<uint16_t>> s;
    vector<uint16_t> seen(nodes.size(), 0);
    deque<uint16_t> d;

    for (auto c : nodes)
    {
        if (seen[c.second->id] == 0 and c.second->out_nodes.size() == 2)
        {
            d = {c.second->id};
            extend_unitig(d);
            for (auto i : d)
            {
                seen[i] = 1;
            }
            s.insert(d);
        } else {
            seen[c.second->id] = 1;
        }
    }
    return s;
}

// extend a dbg path on either end to a branch point
void Graph::extend_unitig(deque<uint16_t>& tig)
{
    if (tig.size() == 0)
    {
        return;
    } else if (tig.size() == 1 and nodes[tig.back()]->out_nodes.size() == 2)
    {
        tig.push_front(*nodes[tig.back()]->out_nodes.begin());
        tig.push_back(*nodes[tig.back()]->out_nodes.begin()+1);
    }
    while (nodes[tig.back()]->out_nodes.size() == 2)
    {
        if (*nodes[tig.back()]->out_nodes.begin() == tig[tig.size()-2])
        {
            tig.push_back(*nodes[tig.back()]->out_nodes.begin()+1);
        } else if (*nodes[tig.back()]->out_nodes.begin()+1 == tig[tig.size()-2])
        {
            tig.push_back(*nodes[tig.back()]->out_nodes.begin());
        }
        // else error?
    }
    while (nodes[tig.front()]->out_nodes.size() == 2)
    {
        if (*nodes[tig.front()]->out_nodes.begin() == tig[1])
        {
            tig.push_front(*nodes[tig.front()]->out_nodes.begin()+1);
        } else if (*nodes[tig.front()]->out_nodes.begin()+1 == tig[1])
        {
            tig.push_front(*nodes[tig.front()]->out_nodes.begin());
        }
        // else error?
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
