#include <deque>
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

// add a node in dbg corresponding to a fixed size deque of pangenome graph
// node/orientation ids and labelled with the read_ids which cover it
NodePtr Graph::add_node (const deque<uint16_t>& node_ids, uint32_t read_id)
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
    //cout << "new node " << next_id << endl;
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
        if (seen[c.second->id] == 0 and c.second->out_nodes.size() <= 2)
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
        tig.push_back(*++nodes[tig.back()]->out_nodes.begin());
    } else if (tig.size() == 1 and nodes[tig.back()]->out_nodes.size() == 1)
    {
        tig.push_front(*nodes[tig.back()]->out_nodes.begin());
    }

    while (nodes[tig.back()]->out_nodes.size() == 2)
    {

        if (*nodes[tig.back()]->out_nodes.begin() == tig[tig.size()-2])
        {
            tig.push_back(*++nodes[tig.back()]->out_nodes.begin());
        } else if (*++nodes[tig.back()]->out_nodes.begin() == tig[tig.size()-2])
        {
            tig.push_back(*nodes[tig.back()]->out_nodes.begin());
        } else {
            break;
        }
        // else error?
    }
    cout << "tig front " << tig.front() << endl;
    while (nodes[tig.front()]->out_nodes.size() == 2)
    {
        if (*nodes[tig.front()]->out_nodes.begin() == tig[1])
        {
            tig.push_front(*++nodes[tig.front()]->out_nodes.begin());
        } else if (*++nodes[tig.front()]->out_nodes.begin() == tig[1])
        {
            tig.push_front(*nodes[tig.front()]->out_nodes.begin());
        } else {
            break;
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
        bool out_found, found = false;
        for (const auto s : y.nodes) {
            if (*t.second == *s.second) {
                found = true;

                // also check the outnodes are the same
                if (t.second->out_nodes.size() != s.second->out_nodes.size())
                {
                    return false;
                }
                for (const auto i : t.second->out_nodes) {
                    out_found = false;
                    for (const auto j : s.second->out_nodes) {
                        if (*nodes.at(i) == *y.nodes.at(j))
                        {
                            out_found = true;
                            break;
                        }
                    }
                    if (out_found == false)
                    {
                        return false;
                    }
                }

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

namespace debruijn {
    std::ostream &operator<<(std::ostream &out, const Graph &m) {
        for (const auto &n : m.nodes) {
            out << n.first << endl;
            for (const auto &o : n.second->out_nodes) {
                out << n.first << " -> " << o << endl;
            }
        }
        return out;
    }
}
