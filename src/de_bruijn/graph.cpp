#include <deque>
#include <unordered_set>
#include <set>
#include <iostream>
#include <memory>
#include <algorithm>

#include <boost/log/trivial.hpp>

#include "de_bruijn/graph.h"
#include "noise_filtering.h"

using namespace debruijn;

// Define a debruijn graph with s-mers of genes as nodes
debruijn::Graph::Graph(uint8_t s)
    : next_id(0)
    , size(s)
{
    nodes.reserve(200000);
};

debruijn::Graph::~Graph() { nodes.clear(); }

// Add a node in dbg corresponding to a fixed size deque of pangenome graph
// node/orientation ids and labelled with the read_ids which cover it
OrientedNodePtr debruijn::Graph::add_node(
    const std::deque<uint_least32_t>& node_ids, uint32_t read_id)
{
    const bool correct_number_of_nodes_to_add = node_ids.size() == size;
    if(!correct_number_of_nodes_to_add) {
        fatal_error("Error adding node to de Bruijn Graph: expected node of size ", size,
                    ", received node of size ", node_ids.size());
    }

    if (node_hash.find(node_ids) != node_hash.end()) {
        nodes[node_hash[node_ids]]->read_ids.insert(read_id);
        return make_pair(nodes[node_hash[node_ids]], true);
    } else if (node_hash.find(rc_hashed_node_ids(node_ids)) != node_hash.end()) {
        auto rc = rc_hashed_node_ids(node_ids);
        nodes[node_hash[rc]]->read_ids.insert(read_id);
        return make_pair(nodes[node_hash[rc]], false);
    }

    NodePtr n;
    n = std::make_shared<Node>(next_id, node_ids, read_id);
    nodes[next_id] = n;
    node_hash[node_ids] = next_id;

    if (next_id % 1000 == 0) {
        BOOST_LOG_TRIVIAL(debug) << "added node " << next_id;
    }

    next_id++;
    return make_pair(n, true);
}

// An edge is valid if the kmer of node/orientation ids for from
// overlaps the first k-1 nodes/orientations of to
// Note that forward and reverse complement kmers are treated as
// equal, but an edge is only valid if the orientation they were found
// in in the read allows the overlap
bool edge_is_valid(OrientedNodePtr from, OrientedNodePtr to)
{
    std::deque<uint_least32_t> hashed_node_ids_from = from.first->hashed_node_ids;
    std::deque<uint_least32_t> hashed_node_ids_to = to.first->hashed_node_ids;
    if (!from.second) {
        hashed_node_ids_from = rc_hashed_node_ids(hashed_node_ids_from);
    }

    if (!to.second) {
        hashed_node_ids_to = rc_hashed_node_ids(hashed_node_ids_to);
    }

    return overlap_forwards(hashed_node_ids_from, hashed_node_ids_to);
}

// Add directed edge between from and to
void debruijn::Graph::add_edge(OrientedNodePtr from, OrientedNodePtr to)
{
    bool nodes_are_valid = from.first != nullptr and to.first != nullptr;
    if(!nodes_are_valid) {
        fatal_error("Error adding edge to de Bruijn Graph: from or to node is invalid");
    }

    if (!edge_is_valid(from, to)) {
        fatal_error("Error adding edge to de Bruijn Graph: edge from ", *from.first,
                    " to ", *to.first, " is invalid");
    }

    if (from.second
        and from.first->out_nodes.find(to.first->id) == from.first->out_nodes.end()) {
        from.first->out_nodes.insert(to.first->id);
    } else if (!from.second
        and from.first->in_nodes.find(to.first->id) == from.first->in_nodes.end()) {
        from.first->in_nodes.insert(to.first->id);
    }

    if (to.second
        and to.first->in_nodes.find(from.first->id) == to.first->in_nodes.end()) {
        to.first->in_nodes.insert(from.first->id);
    } else if (!to.second
        and to.first->out_nodes.find(from.first->id) == to.first->out_nodes.end()) {
        to.first->out_nodes.insert(from.first->id);
    }
}

// Remove all mentions of de bruijn node with id given from graph
void debruijn::Graph::remove_node(const uint32_t dbg_node_id)
{
    auto it = nodes.find(dbg_node_id);
    if (it != nodes.end()) {
        // remove this node from lists of out nodes from other graph nodes
        for (const auto& n : it->second->out_nodes) {
            nodes[n]->in_nodes.erase(dbg_node_id);
            nodes[n]->out_nodes.erase(dbg_node_id);
        }
        for (const auto& n : it->second->in_nodes) {
            nodes[n]->out_nodes.erase(dbg_node_id);
            nodes[n]->in_nodes.erase(dbg_node_id);
        }

        // and remove from nodes
        nodes.erase(dbg_node_id);
    }
}

// Remove all copies of read from de bruijn node
void debruijn::Graph::remove_read_from_node(
    const uint32_t read_id, const uint32_t dbg_node_id)
{
    auto it = nodes.find(dbg_node_id);
    bool found_read_intersect;
    if (it != nodes.end()) {
        auto rit = it->second->read_ids.find(read_id);
        if (rit != it->second->read_ids.end()) {
            it->second->read_ids.erase(rit);

            // if there are no more reads covering it, remove the node
            if (it->second->read_ids.empty()) {
                remove_node(dbg_node_id);
            } else {
                // otherwise, remove any outnodes which no longer share a read
                for (std::unordered_set<uint32_t>::iterator nit
                     = it->second->out_nodes.begin();
                     nit != it->second->out_nodes.end();) {
                    found_read_intersect = false;
                    for (const auto& r : it->second->read_ids) {
                        if (nodes[*nit]->read_ids.find(r)
                            != nodes[*nit]->read_ids.end()) {
                            found_read_intersect = true;
                            break;
                        }
                    }
                    if (!found_read_intersect) {
                        nodes[*nit]->in_nodes.erase(dbg_node_id);
                        nit = it->second->out_nodes.erase(nit);
                    } else {
                        nit++;
                    }
                }
                for (std::unordered_set<uint32_t>::iterator nit
                     = it->second->in_nodes.begin();
                     nit != it->second->in_nodes.end();) {
                    found_read_intersect = false;
                    for (const auto& r : it->second->read_ids) {
                        if (nodes[*nit]->read_ids.find(r)
                            != nodes[*nit]->read_ids.end()) {
                            found_read_intersect = true;
                            break;
                        }
                    }
                    if (!found_read_intersect) {
                        nodes[*nit]->out_nodes.erase(dbg_node_id);
                        nit = it->second->in_nodes.erase(nit);
                    } else {
                        nit++;
                    }
                }
            }
        }
    }
}

// Get the dbg node ids corresponding to leaves
std::unordered_set<uint32_t> debruijn::Graph::get_leaves(uint_least32_t covg_thresh)
{
    std::unordered_set<uint32_t> s;
    for (const auto& c : nodes) {
        BOOST_LOG_TRIVIAL(debug)
            << "node " << *c.second << " has " << c.second->out_nodes.size() << " + "
            << c.second->in_nodes.size() << " outnodes";
        if (c.second->read_ids.size() > covg_thresh) {
            continue;
        } else if (c.second->out_nodes.size() + c.second->in_nodes.size() <= 1) {
            s.insert(c.second->id);
        }
    }
    return s;
}

// Get deques of dbg node ids corresponding to maximal non-branching paths in dbg
std::set<std::deque<uint32_t>> debruijn::Graph::get_unitigs()
{
    std::set<std::deque<uint32_t>> all_tigs;
    std::set<uint32_t> seen;

    for (const auto& node_entry : nodes) {
        const auto& id = node_entry.first;
        const auto& node_ptr = node_entry.second;

        bool node_seen = seen.find(id) != seen.end();
        bool at_branch
            = (node_ptr->out_nodes.size() > 1) or (node_ptr->in_nodes.size() > 1);
        if (node_seen or at_branch)
            continue;

        std::deque<uint32_t> tig = { id };
        extend_unitig(tig);
        for (const auto& other_id : tig)
            seen.insert(other_id);
        all_tigs.insert(tig);
    }
    return all_tigs;
}

// Extend a dbg path on either end until reaching a branch point
void debruijn::Graph::extend_unitig(std::deque<uint32_t>& tig)
{
    bool tig_is_empty = (tig.empty());
    bool node_is_isolated = (tig.size() == 1
        and (nodes[tig.back()]->out_nodes.size() + nodes[tig.back()]->in_nodes.size())
            == 0);
    if (tig_is_empty or node_is_isolated) {
        return;
    }

    bool can_extend = nodes[tig.back()]->out_nodes.size() == 1;
    bool use_outnodes = true;
    while (can_extend) {

        if (use_outnodes) {
            tig.push_back(*nodes[tig.back()]->out_nodes.begin());
        } else {
            tig.push_back(*nodes[tig.back()]->in_nodes.begin());
        }

        if (std::find(nodes[tig.back()]->in_nodes.begin(),
                nodes[tig.back()]->in_nodes.end(), *----tig.end())
            != nodes[tig.back()]->in_nodes.end()) {
            can_extend = nodes[tig.back()]->out_nodes.size() == 1
                and nodes[tig.back()]->in_nodes.size() <= 1
                and tig.front() != tig.back();
            use_outnodes = true;
        } else if (std::find(nodes[tig.back()]->out_nodes.begin(),
                       nodes[tig.back()]->out_nodes.end(), *----tig.end())
            != nodes[tig.back()]->out_nodes.end()) {
            can_extend = nodes[tig.back()]->in_nodes.size() == 1
                and nodes[tig.back()]->out_nodes.size() <= 1
                and tig.front() != tig.back();
            use_outnodes = false;
        } else {
            can_extend = false;
        }
    }

    if (tig.size() == 1) {
        can_extend = nodes[tig.front()]->in_nodes.size() == 1
            and nodes[tig.front()]->out_nodes.size() <= 1;
        use_outnodes = false;
    } else {

        if (std::find(nodes[tig.front()]->in_nodes.begin(),
                nodes[tig.front()]->in_nodes.end(), *++tig.begin())
            != nodes[tig.front()]->in_nodes.end()) {
            can_extend = nodes[tig.front()]->out_nodes.size() == 1
                and nodes[tig.front()]->in_nodes.size() <= 1
                and tig.front() != tig.back();
            use_outnodes = true;
        } else if (std::find(nodes[tig.front()]->out_nodes.begin(),
                       nodes[tig.front()]->out_nodes.end(), *++tig.begin())
            != nodes[tig.front()]->out_nodes.end()) {
            can_extend = nodes[tig.front()]->in_nodes.size() == 1
                and nodes[tig.front()]->out_nodes.size() <= 1
                and tig.front() != tig.back();
            use_outnodes = false;
        } else {
            can_extend = false;
        }
    }

    while (can_extend) {
        if (use_outnodes) {
            tig.push_front(*nodes[tig.front()]->out_nodes.begin());
        } else {
            tig.push_front(*nodes[tig.front()]->in_nodes.begin());
        }

        if (std::find(nodes[tig.front()]->in_nodes.begin(),
                nodes[tig.front()]->in_nodes.end(), *++tig.begin())
            != nodes[tig.front()]->in_nodes.end()) {
            can_extend = nodes[tig.front()]->out_nodes.size() == 1
                and nodes[tig.front()]->in_nodes.size() <= 1
                and tig.front() != tig.back();
            use_outnodes = true;
        } else if (std::find(nodes[tig.front()]->out_nodes.begin(),
                       nodes[tig.front()]->out_nodes.end(), *++tig.begin())
            != nodes[tig.front()]->out_nodes.end()) {
            can_extend = nodes[tig.front()]->in_nodes.size() == 1
                and nodes[tig.front()]->out_nodes.size() <= 1
                and tig.front() != tig.back();
            use_outnodes = false;
        } else {
            can_extend = false;
        }
    }

    while (tig.size() > 1 and tig.front() == tig.back()) {
        tig.pop_back();
    }

    std::stringstream tig_ss;
    for (const auto& n : tig)
        tig_ss << n << " ";
    BOOST_LOG_TRIVIAL(debug) << "got tig of length " << tig.size() << ": "
                             << tig_ss.str();
}

// Search the outnodes of node_ptr_to_search for node_ptr_to_find
bool debruijn::Graph::found_in_out_nodes(
    const NodePtr node_ptr_to_search, const NodePtr node_ptr_to_find) const
{
    for (const auto& i : node_ptr_to_search->out_nodes) {
        if (*nodes.at(i) == *node_ptr_to_find) {
            return true;
        }
    }
    return false;
}

// Search the innodes of node_ptr_to_search for node_ptr_to_find
bool debruijn::Graph::found_in_in_nodes(
    const NodePtr node_ptr_to_search, const NodePtr node_ptr_to_find) const
{
    for (const auto& i : node_ptr_to_search->in_nodes) {
        if (*nodes.at(i) == *node_ptr_to_find) {
            return true;
        }
    }
    return false;
}

// Graphs are equal if they have nodes corresponding to the same kmer
// of node/orientations, with outnodes and innodes the same
bool debruijn::Graph::operator==(const Graph& y) const
{
    // want the graphs to have the same nodes, even if
    // the ids given them is different.
    if (nodes.size() != y.nodes.size()) {
        BOOST_LOG_TRIVIAL(debug) << "different num nodes";
        return false;
    }

    for (const auto& t : nodes) {
        bool found = false;
        for (const auto& s : y.nodes) {
            if (*t.second == *s.second) {
                found = true;

                // also check the outnodes are the same
                if (t.second->out_nodes.size() + t.second->in_nodes.size()
                    != s.second->out_nodes.size() + s.second->in_nodes.size()) {
                    BOOST_LOG_TRIVIAL(debug) << "node has different number of outnodes";
                    return false;
                }

                for (const auto& i : t.second->out_nodes) {
                    if (not y.found_in_out_nodes(s.second, nodes.at(i))
                        and not y.found_in_in_nodes(s.second, nodes.at(i))) {
                        BOOST_LOG_TRIVIAL(debug)
                            << "did not find " << i << " in outnode or innodes";
                        return false;
                    }
                }
                for (const auto& i : t.second->in_nodes) {
                    if (not y.found_in_out_nodes(s.second, nodes.at(i))
                        and not y.found_in_in_nodes(s.second, nodes.at(i))) {
                        BOOST_LOG_TRIVIAL(debug)
                            << "did not find " << i << " in outnode or innodes";
                        return false;
                    }
                }
                break;
            }
        }
        if (!found) {
            BOOST_LOG_TRIVIAL(debug)
                << "did not find node " << t.first << " " << *t.second;
            return false;
        }
    }

    // nodes can't be equal within a graph so don't need to check vice versa
    return true;
}

bool debruijn::Graph::operator!=(const Graph& y) const { return !(*this == y); }

namespace debruijn {
std::ostream& operator<<(std::ostream& out, const Graph& m)
{
    for (const auto& n : m.nodes) {
        out << n.first << ": " << *(n.second) << std::endl;
        for (const auto& o : n.second->out_nodes) {
            out << n.first << " -> " << o << std::endl;
        }
    }
    return out;
}
}
