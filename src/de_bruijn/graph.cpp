#include <deque>
#include <limits>
#include <fstream>
#include <unordered_set>
#include <set>
#include <deque>
#include <iostream>
#include <memory>
#include <cassert>
#include "de_bruijn/graph.h"
#include "de_bruijn/node.h"
#include "noise_filtering.h"
#include "utils.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace debruijn;

Graph::Graph(uint8_t s) : next_id(0), size(s) {
    nodes.reserve(200000);
};


Graph::~Graph()
{
    nodes.clear();
}

// add a node in dbg corresponding to a fixed size deque of pangenome graph
// node/orientation ids and labelled with the read_ids which cover it
OrientedNodePtr Graph::add_node (const deque<uint16_t>& node_ids, uint32_t read_id)
{
    assert(node_ids.size() == size);

    if (node_hash.find(node_ids) != node_hash.end()) {
        nodes[node_hash[node_ids]]->read_ids.insert(read_id);
        return make_pair(nodes[node_hash[node_ids]],true);
    } else if (node_hash.find(rc_hashed_node_ids(node_ids)) != node_hash.end()) {
        auto rc = rc_hashed_node_ids(node_ids);
        nodes[node_hash[rc]]->read_ids.insert(read_id);
        return make_pair(nodes[node_hash[rc]],false);
    }

    NodePtr n;
    n = make_shared<Node>(next_id, node_ids, read_id);
    assert(n!=nullptr);
    nodes[next_id] = n;
    node_hash[node_ids] = next_id;

    if (next_id%1000==0)
    {
        cout << "added node " << next_id << endl;
    }

    next_id++;
    assert(next_id < numeric_limits<uint16_t>::max()||assert_msg("WARNING, reached max de bruijn graph node size"));
    return make_pair(n,true);
}

bool edge_is_valid (OrientedNodePtr from, OrientedNodePtr to)
{
    deque<uint16_t> hashed_node_ids_from = from.first->hashed_node_ids;
    deque<uint16_t> hashed_node_ids_to = to.first->hashed_node_ids;
    if (from.second == false) {
        hashed_node_ids_from = rc_hashed_node_ids(hashed_node_ids_from);
	    //cout << "reverse from" << endl;
    }

    if (to.second == false) {
        hashed_node_ids_to = rc_hashed_node_ids(hashed_node_ids_to);
	    //cout << "reverse to" << endl;
    }
    /*cout << from.second << to.second << " compare (";
    for (auto n : hashed_node_ids_from)
    {
	cout << n << " ";
    } 
    cout << ") to (";
    for (auto n : hashed_node_ids_to)
    {
        cout << n << " ";
    }
    cout << ")" << endl;*/

    return overlap_forwards(hashed_node_ids_from, hashed_node_ids_to);
}

// add directed edge
void Graph::add_edge (OrientedNodePtr from, OrientedNodePtr to)
{
    assert(from.first != nullptr and to.first != nullptr);
    assert(edge_is_valid(from,to) or assert_msg("edge from " << *from.first << " to " << *to.first << " is invalid"));

    //uint8_t num_edges_added = 0;
    if (from.second and from.first->out_nodes.find(to.first->id) == from.first->out_nodes.end())
    {
	    from.first->out_nodes.insert(to.first->id);
        cout << "added edge " << from.first->id << " -> " << to.first->id << " so out_nodes.size() == " << from.first->out_nodes.size() << endl;
	    //num_edges_added += 1;
    } else if (!from.second and from.first->in_nodes.find(to.first->id) == from.first->in_nodes.end()){
	    from.first->in_nodes.insert(to.first->id);
        cout << "added edge " << from.first->id << " <- " << to.first->id << " so in_nodes.size() == " << from.first->in_nodes.size() << endl;
	    //num_edges_added += 1;
    }

    if (to.second and to.first->in_nodes.find(from.first->id) == to.first->in_nodes.end())
    {
        to.first->in_nodes.insert(from.first->id);
        cout << "added edge " << to.first->id << " <- " << from.first->id << " so in_nodes.size() == " << to.first->in_nodes.size() << endl;
        //num_edges_added += 1;
    } else if (!to.second and to.first->out_nodes.find(from.first->id) == to.first->out_nodes.end()){
        to.first->out_nodes.insert(from.first->id);
        cout << "added edge " << to.first->id << " -> " << from.first->id << " so out_nodes.size() == " << to.first->out_nodes.size() << endl;
        //num_edges_added += 1;
    }

    //assert(num_edges_added == 2 or assert_msg("did not add edge from " << *from << " to " << *to));
}

// remove de bruijn node with id given
void Graph::remove_node(const uint32_t dbg_node_id)
{
    //cout << "remove node " << dbg_node_id << endl;
    auto it = nodes.find(dbg_node_id);
    if ( it != nodes.end())
    {
        // remove this node from lists of out nodes from other graph nodes
        for (auto n : it->second->out_nodes)
        {
            nodes[n]->in_nodes.erase(dbg_node_id);
	    nodes[n]->out_nodes.erase(dbg_node_id);
        }
        for (auto n : it->second->in_nodes)
        {
            nodes[n]->out_nodes.erase(dbg_node_id);
	    nodes[n]->in_nodes.erase(dbg_node_id);
        }

        // and remove from nodes
        nodes.erase(dbg_node_id);
    }
}

// remove unitig from dbg
/*void Graph::remove_unitig(const deque<uint32_t>& tig)
{
    if (tig.size() < size)
    {
        return;
    }

    // first remove nodes
    unordered_set<uint32_t> read_ids;
    bool orientation = overlap_forwards(nodes[tig[0]]->hashed_node_ids, nodes[tig[1]]->hashed_node_ids);
    for (uint i=1; i<nodes.size()-1; ++i)
    {
        read_ids.insert(nodes[tig[i]]->read_ids.begin(), nodes[tig[i]]->read_ids.end());
        orientation *= overlap_forwards(nodes[tig[i]]->hashed_node_ids, nodes[tig[i+1]]->hashed_node_ids);
        remove_node[tig[i]];
    }

    // then add nodes and edges back to allow continuity
    NodePtr last_node = nodes[tig[0]];
    NodePtr new_node;
    for (uint i=1; i<size; ++i)
    {
        deque<uint16_t> new_tig;
        for (uint j=i; j<size; ++j) {
            new_tig.push_back(nodes[tig.front()]->hashed_node_ids[i]);
        }
        for (uint j=size-i; j<size; ++j) {
            if (orientation)
            {
                new_tig.push_back(nodes[tig.back()]->hashed_node_ids[i]);
            } else {
                new_tig.push_back(rc_hashed_node_ids(nodes[tig.back()]->hashed_node_ids)[i]);
            }
        }
        assert(new_tig.size() == size);
        for(auto r : read_ids)
            new_node = add_node(new_tig, r);
        add_edge(last_node, new_node);
        last_node = new_node;
    }
    add_edge(last_node, nodes[tig.back()]);
}*/

// remove read from de bruijn node
void Graph::remove_read_from_node(const uint32_t read_id, const uint32_t dbg_node_id)
{
    //cout << "remove read " << (uint)read_id << " from node " << (uint)dbg_node_id << endl;
    auto it = nodes.find(dbg_node_id);
    bool found_read_intersect;
    if ( it != nodes.end())
    {
        //cout << "found node " << *(it->second) << endl;
        auto rit = it->second->read_ids.find(read_id);
        if (rit != it->second->read_ids.end())
        {
            //cout << "found read id on node" << endl;
            it->second->read_ids.erase(rit);
            //cout << "removed read id" << endl;

            // if there are no more reads covering it, remove the node

            if (it->second->read_ids.empty())
            {
                //cout << "no more reads, so remove node" << endl;
                remove_node(dbg_node_id);
            } else {
                // otherwise, remove any outnodes which no longer share a read
                //cout << "remove outnodes which no longer share a read" << endl;
                for (unordered_set<uint32_t>::iterator nit = it->second->out_nodes.begin();
                     nit != it->second->out_nodes.end();)
                {
                    //cout << "out node " << *nit << endl;
                    found_read_intersect = false;
                    for (auto r : it->second->read_ids) {
                        if (nodes[*nit]->read_ids.find(r) != nodes[*nit]->read_ids.end()) {
                            found_read_intersect = true;
                            //cout << " shares a read" << endl;
                            break;
                        }
                    }
                    if (found_read_intersect == false)
                    {
                        //cout << " does not share a read" << endl;
                        nodes[*nit]->in_nodes.erase(dbg_node_id);
                        nit = it->second->out_nodes.erase(nit);
                        //cout << "removed" << endl;
                    } else {
                        nit++;
                    }
                }
                for (unordered_set<uint32_t>::iterator nit = it->second->in_nodes.begin();
                     nit != it->second->in_nodes.end();)
                {
                    //cout << "out node " << *nit << endl;
                    found_read_intersect = false;
                    for (auto r : it->second->read_ids) {
                        if (nodes[*nit]->read_ids.find(r) != nodes[*nit]->read_ids.end()) {
                            found_read_intersect = true;
                            //cout << " shares a read" << endl;
                            break;
                        }
                    }
                    if (found_read_intersect == false)
                    {
                        //cout << " does not share a read" << endl;
                        nodes[*nit]->out_nodes.erase(dbg_node_id);
                        nit = it->second->in_nodes.erase(nit);
                        //cout << "removed" << endl;
                    } else {
                        nit++;
                    }
                }
            }
        }
    }
}

// get the dbg node ids corresponding to leaves
unordered_set<uint32_t> Graph::get_leaves(uint16_t covg_thresh)
{
    unordered_set<uint32_t> s;
    for (auto c : nodes)
    {
	cout << "node " << *c.second << " has " << c.second->out_nodes.size() << " + " << c.second->in_nodes.size() << " outnodes" << endl;
        if (c.second->read_ids.size() > covg_thresh) {
            continue;
        } else if (c.second->out_nodes.size() + c.second->in_nodes.size() <= 1) {
            s.insert(c.second->id);
        }
    }
    return s;
}

// get deques of dbg node ids corresponding to maximal non-branching paths in dbg
set<deque<uint32_t>> Graph::get_unitigs() {
    set<deque<uint32_t>> all_tigs;
    set<uint32_t> seen;

    for (auto node_entry : nodes) {
        const auto &id = node_entry.first;
        assert(id <= next_id);

        const auto &node_ptr = node_entry.second;
        assert(node_ptr != nullptr);

        bool node_seen = seen.find(id) != seen.end();
        bool at_branch = (node_ptr->out_nodes.size() > 1) or (node_ptr->in_nodes.size() > 1);
        if (node_seen or at_branch)
            continue;

        deque<uint32_t> tig = {id};
        extend_unitig(tig);
        for (auto other_id: tig)
            seen.insert(other_id);
        all_tigs.insert(tig);
    }
    return all_tigs;
}

// extend a dbg path on either end to a branch point
void Graph::extend_unitig(deque<uint32_t>& tig) {
    bool tig_is_empty = (tig.size() == 0);
    bool node_is_isolated = (tig.size() == 1
                             and nodes[tig.back()]->out_nodes.size() + nodes[tig.back()]->in_nodes.size() == 0);
    if (tig_is_empty or node_is_isolated) {
        //cout << "node is isolated or tig empty" << endl;
        return;
    }

    bool can_extend = nodes[tig.back()]->out_nodes.size() == 1;// and nodes[tig.back()]->in_nodes.size() <= 1;
    bool use_outnodes = true;
    while (can_extend) {

        /*cout << "tig in progress before b: ";
        for (auto n : tig) {
            cout << n << " ";
        }
        cout << endl;*/

        if (use_outnodes)// and find(tig.begin(), tig.end(), *nodes[tig.back()]->out_nodes.begin())!=tig.end())
            tig.push_back(*nodes[tig.back()]->out_nodes.begin());
        else //if (find(tig.begin(), tig.end(), *nodes[tig.back()]->in_nodes.begin())!=tig.end())
            tig.push_back(*nodes[tig.back()]->in_nodes.begin());

        if (find(nodes[tig.back()]->in_nodes.begin(), nodes[tig.back()]->in_nodes.end(),
                 *----tig.end()) != nodes[tig.back()]->in_nodes.end()) {
            can_extend = nodes[tig.back()]->out_nodes.size() == 1
                         and nodes[tig.back()]->in_nodes.size() <= 1
                         and tig.front() != tig.back();
            use_outnodes = true;
            //cout << "A";
        } else if (find(nodes[tig.back()]->out_nodes.begin(), nodes[tig.back()]->out_nodes.end(),
                        *----tig.end()) != nodes[tig.back()]->out_nodes.end()) {
            can_extend = nodes[tig.back()]->in_nodes.size() == 1
                         and nodes[tig.back()]->out_nodes.size() <= 1
                         and tig.front() != tig.back();
            use_outnodes = false;
            //cout << "B";
        } else {
            can_extend = false;
            //cout << "C";
        }

        /*cout << "tig in progress b: ";
        for (auto n : tig) {
            cout << n << " ";
        }
        cout << endl;*/
    }

    if (tig.size() == 1) {
        can_extend = nodes[tig.front()]->in_nodes.size() == 1 and nodes[tig.front()]->out_nodes.size() <= 1;
        use_outnodes = false;
        //cout << "D";
    } else {

        if (find(nodes[tig.front()]->in_nodes.begin(), nodes[tig.front()]->in_nodes.end(),
                 *++tig.begin()) != nodes[tig.front()]->in_nodes.end()) {
            can_extend = nodes[tig.front()]->out_nodes.size() == 1
                         and nodes[tig.front()]->in_nodes.size() <= 1
                         and tig.front() != tig.back();
            use_outnodes = true;
            //cout << "E";
        } else if (find(nodes[tig.front()]->out_nodes.begin(), nodes[tig.front()]->out_nodes.end(),
                        *++tig.begin()) != nodes[tig.front()]->out_nodes.end()) {
            can_extend = nodes[tig.front()]->in_nodes.size() == 1
                         and nodes[tig.front()]->out_nodes.size() <= 1
                         and tig.front() != tig.back();
            use_outnodes = false;
            //cout << "F";
        } else {
            can_extend = false;
        }
    }

    while (can_extend)
    {
        /*cout << "tig in progress before f: ";
        for (auto n : tig) {
            cout << n << " ";
        }
        cout << endl;*/

        if (use_outnodes)// and find(tig.begin(), tig.end(), *nodes[tig.front()]->out_nodes.begin())!=tig.end())
            tig.push_front(*nodes[tig.front()]->out_nodes.begin());
        else //if (find(tig.begin(), tig.end(), *nodes[tig.front()]->in_nodes.begin())!=tig.end())
            tig.push_front(*nodes[tig.front()]->in_nodes.begin());

        if (find(nodes[tig.front()]->in_nodes.begin(), nodes[tig.front()]->in_nodes.end(),
                 *++tig.begin()) != nodes[tig.front()]->in_nodes.end()) {
            can_extend = nodes[tig.front()]->out_nodes.size() == 1
                         and nodes[tig.front()]->in_nodes.size() <= 1
                         and tig.front() != tig.back();
            use_outnodes = true;
            //cout << "G";
        } else if (find(nodes[tig.front()]->out_nodes.begin(), nodes[tig.front()]->out_nodes.end(),
                        *++tig.begin()) != nodes[tig.front()]->out_nodes.end()) {
            can_extend = nodes[tig.front()]->in_nodes.size() == 1
                         and nodes[tig.front()]->out_nodes.size() <= 1
                         and tig.front() != tig.back();
            use_outnodes = false;
            //cout << "H";
        } else {
            can_extend = false;
            //cout << "I";
        }
        
        /*cout << "tig in progress f: ";
        for (auto n : tig)
        {   
            cout << n << " ";
        }
        cout << endl;*/
    }

    while (tig.size() > 1 and tig.front() == tig.back())//find(tig.begin(), --tig.end(), tig.back()) == tig.end())
        tig.pop_back();

    cout << "got tig of length " << tig.size() << ": ";
    for (auto n : tig)
    {
        cout << n << " ";
    }
    cout << endl;
}

bool Graph::found_in_out_nodes(const NodePtr node_ptr_to_search, const NodePtr node_ptr_to_find) const
{
    for (const auto i : node_ptr_to_search->out_nodes) {
        if (*nodes.at(i) == *node_ptr_to_find)
        {
            return true;
        } else {
            cout << *nodes.at(i) << " != " << *node_ptr_to_find << endl;
        }
    }
    return false;
}

bool Graph::found_in_in_nodes(const NodePtr node_ptr_to_search, const NodePtr node_ptr_to_find) const
{
    for (const auto i : node_ptr_to_search->in_nodes) {
        if (*nodes.at(i) == *node_ptr_to_find)
        {
            return true;
        } else {
            cout << *nodes.at(i) << " != " << *node_ptr_to_find << endl;
        }
    }
    return false;
}

bool Graph::operator == (const Graph& y) const
{
    // want the graphs to have the same nodes, even if
    // the ids given them is different.
    if (nodes.size() != y.nodes.size())
    {
        cout << "different num nodes" << endl;
	    return false;
    }

    for (const auto t : nodes) {
        bool found = false;
        for (const auto s : y.nodes) {
            if (*t.second == *s.second) {
                //cout << t.first << " " << *t.second << " == " << *s.second << " " << s.first << endl;

                found = true;

                // also check the outnodes are the same
                if (t.second->out_nodes.size() + t.second->in_nodes.size() != s.second->out_nodes.size()+s.second->in_nodes.size())
                {
                    cout << "node has different number of outnodes" << endl;
                    return false;
                }

                for (const auto i : t.second->out_nodes) {
                    if (not y.found_in_out_nodes(s.second, nodes.at(i))
                        and not y.found_in_in_nodes(s.second, nodes.at(i))){
                        cout << "did not find " << i << " in outnode or innodes " << endl;
                        return false;
                    }
                }
                for (const auto i : t.second->in_nodes) {
                    if (not y.found_in_out_nodes(s.second, nodes.at(i))
                        and not y.found_in_in_nodes(s.second, nodes.at(i))){
                        cout << "did not find " << i << " in outnode or innodes " << endl;
                        return false;
                    }
                }
                break;
            }
        }
        if (found == false) {
            cout << "did not find node " << t.first << " " << *t.second << endl;
            return false;
        }
    }

    // nodes can't be equal within a graph so don't need to check vice versa
    //cout << "found all nodes" << endl;
    return true;
}

bool Graph::operator!=(const Graph &y) const {
    return !(*this == y);
}

namespace debruijn {
    std::ostream &operator<<(std::ostream &out, const Graph &m) {
        for (const auto &n : m.nodes) {
            out << n.first << ": " << *(n.second) << endl;
            for (const auto &o : n.second->out_nodes) {
                out << n.first << " -> " << o << endl;
            }
        }
        return out;
    }
}
