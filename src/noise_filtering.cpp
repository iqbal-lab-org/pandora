#include <iostream>
#include <map>
#include <unordered_set>
#include <vector>
#include <fstream>
#include "utils.h"
#include "pangenome/pangraph.h"
#include "pangenome/panread.h"
#include "de_bruijn/graph.h"
#include <cassert>
#include "minihit.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

uint16_t node_plus_orientation_to_num (const uint16_t node_id, const bool orientation)
{
    assert(node_id < 32768);
    uint16_t r = 2*node_id;
    if (orientation)
    {
        r+=1;
    }
    return r;
}

void num_to_node_plus_orientation (uint16_t& node_id, bool& orientation, const uint16_t num)
{
    if (num % 2 == 1)
    {
        orientation = true;
        node_id = (num - 1)/2;
    } else {
        orientation = false;
        node_id = num/2;
    }
}

void hashed_node_ids_to_ids_and_orientations(const deque<uint16_t>& hashed_node_ids,
                                             vector<uint16_t>& node_ids,
                                             vector<bool>& node_orients)
{
    uint16_t node_id;
    bool orientation;
    for (auto i : hashed_node_ids)
    {
        num_to_node_plus_orientation(node_id, orientation, i);
        node_ids.push_back(node_id);
        node_orients.push_back(orientation);
    }
}

debruijn::Graph construct_debruijn_graph_from_pangraph(uint8_t size, const pangenome::Graph & pg)
{
    cout << now() << "Construct de Bruijn Graph from PanGraph with size " << size << endl;
    debruijn::Graph dbg(size);

    debruijn::NodePtr prev, current;
    deque<uint16_t> hashed_ids;

    for (auto r : pg.reads)
    {
        if (r.second->nodes.size() < size)
        {
            // can't add anything for this read
            continue;
        }

        prev = nullptr;
        hashed_ids.clear();

        for (uint i=0; i<r.second->nodes.size(); ++i)
        {
            hashed_ids.push_back(node_plus_orientation_to_num(r.second->nodes[i]->node_id,
                                                              r.second->node_orientations[i]));

            if (hashed_ids.size() == size)
            {
                current = dbg.add_node(hashed_ids, r.first);
                if (prev != nullptr)
                {
                    dbg.add_edge(prev, current);
                }
                prev = current;
                hashed_ids.pop_front();
            }
        }
    }
    return dbg;
}

void remove_leaves(pangenome::Graph & pg, debruijn::Graph & dbg)
{
    bool leaves_exist = true;
    unordered_set<uint16_t> leaves;
    vector<uint16_t> node_ids;
    vector<bool> node_orients;
    uint pos;
    pangenome::NodePtr node;

    while (leaves_exist) {

        leaves = dbg.get_leaves();
        cout << "there are " << leaves.size() << " leaves" << endl;

        if (leaves.empty())
        {
            leaves_exist = false;
        }
        cout << "leaves exist is " << leaves_exist << endl;

        for (auto i : leaves) {
            cout << endl << "looking at leaf " << i << ": ";
            for (const auto j : dbg.nodes[i]->hashed_node_ids)
            {
                cout << j << " ";
            }
            cout << endl;

            // look up the node ids and orientations associated with this node
            hashed_node_ids_to_ids_and_orientations(dbg.nodes[i]->hashed_node_ids, node_ids, node_orients);

            // remove the last node from corresponding reads
            assert(dbg.nodes[i]->read_ids.size()>0);
            for (auto r : dbg.nodes[i]->read_ids)
            {
                if (pg.reads[r]->nodes.size() == dbg.size)
                {
                    pg.remove_read(r);
                } else {
                    pos = pg.reads[r]->find_position(node_ids, node_orients);
                    assert(pos == 0 or pos + node_ids.size() == pg.reads[r]->nodes.size());
                    if (pos == 0) {
                        node = pg.reads[r]->nodes[0];
                        pg.reads[r]->remove_node(pg.reads[r]->nodes.begin());
                        node->remove_read(pg.reads[r]);
                    } else if (pos + node_ids.size() == pg.reads[r]->nodes.size()) {
                        node = pg.reads[r]->nodes.back();
                        pg.reads[r]->remove_node(--pg.reads[r]->nodes.end());
                        node->remove_read(pg.reads[r]);
                    }
                }
            }
            if (node->covg == 0)
            {
                pg.remove_node(node);
            }

            //clean node ids/orients
            node_ids.clear();
            node_orients.clear();

            // remove dbg node
            dbg.remove_node(i);

            cout << "pg is now: " << endl << pg << endl;
        }
    }
}

void filter_unitigs(pangenome::Graph & pg, debruijn::Graph & dbg, const uint16_t& threshold)
{
    set<deque<uint16_t>> unitigs = dbg.get_unitigs();
    for (auto d : unitigs)
    {
        // if a tig has low coverage?
        // for middle nodes
        // reverse look up ids
        // remove from reads
    }
}

void detangle_pangraph_with_debruijn_graph(pangenome::Graph & pg, debruijn::Graph & dbg)
{
    vector<uint16_t> node_ids;
    vector<bool> node_orients;
    bool all_reads_tig;
    unordered_set<pangenome::ReadPtr> reads_along_tig;

    set<deque<uint16_t>> unitigs = dbg.get_unitigs();
    for (auto d : unitigs)
    {
        // look up the node ids and orientations associated with this node
        hashed_node_ids_to_ids_and_orientations(d, node_ids, node_orients);

        // for each node on tig, for each read covering that node, if read doesn't lie along whole tig, create a new node
        for (auto n : node_ids)
        {
            all_reads_tig = true;
            for (auto r : pg.nodes[n]->reads)
            {
                if (reads_along_tig.find(r) == reads_along_tig.end()) {
                    if (r->find_position(node_ids, node_orients) == std::numeric_limits<uint>::max()) {
                        all_reads_tig = false;
                    } else {
                        reads_along_tig.insert(r);
                    }
                }
            }
            if (!all_reads_tig)
            {
                pg.split_node_by_reads(reads_along_tig, node_ids, node_orients, n);
            }
        }
    }
}

void clean_pangraph_with_debruijn_graph(pangenome::Graph & pg, debruijn::Graph & dbg, const uint16_t threshold)
{
    remove_leaves(pg, dbg);
    filter_unitigs(pg, dbg, threshold);
    detangle_pangraph_with_debruijn_graph(pg,dbg);
}