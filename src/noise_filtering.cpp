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

uint16_t rc_num(const uint16_t& num)
{
    return num + 1*(num%2 == 0) - 1*(num%2 == 1);
}

void hashed_node_ids_to_ids_and_orientations(const deque<uint16_t>& hashed_node_ids,
                                             vector<uint16_t>& node_ids,
                                             vector<bool>& node_orients)
{
    node_ids.clear();
    node_orients.clear();

    uint16_t node_id;
    bool orientation;
    for (auto i : hashed_node_ids)
    {
        num_to_node_plus_orientation(node_id, orientation, i);
        node_ids.push_back(node_id);
        node_orients.push_back(orientation);
    }
}

bool overlap_forwards(const deque<uint16_t>& node1, const deque<uint16_t>& node2)
{
    // second deque should extend first by 1
    assert(node1.size()>=node2.size());
    uint i=node1.size()-node2.size()+1;
    uint j=0;
    while(i<node1.size() and j<node2.size())
    {
        if (node1[i] != node2[j])
        {
            return false;
        }
        i++;
        j++;
    }
    return true;
}

bool overlap_backwards(const deque<uint16_t>& node1, const deque<uint16_t>& node2)
{
    for (uint i=1; i<min(node1.size()+1, node2.size()); ++i)
    {
        if (node2[i] != node1[i-1])
        {
            return false;
        }
    }
    return true;
}

deque<uint16_t> reverse_hashed_node(const deque<uint16_t>& node)
{
    deque<uint16_t> d;
    for (auto i : node)
    {
        d.push_front(rc_num(i));
    }
    return d;
}


void dbg_node_ids_to_ids_and_orientations(const debruijn::Graph & dbg,
                                          const deque<uint32_t>& dbg_node_ids,
                                          vector<uint16_t>& node_ids,
                                          vector<bool>& node_orients)
{
    cout << "convert dbg node ids to read ids" << endl;
    node_ids.clear();
    node_orients.clear();

    if (dbg_node_ids.empty())
    {
        return;
    }

    deque<uint16_t> hashed_pg_node_ids = dbg.nodes.at(dbg_node_ids.at(0))->hashed_node_ids;
    deque<uint16_t> rev_node;
    //cout << "hashed pg length " << hashed_pg_node_ids.size() << endl;
    for (uint i=1; i<dbg_node_ids.size(); ++i)
    {
        /*cout << "hashed dbg is currently ";
        for (auto j : hashed_pg_node_ids)
        {
            cout << j << " ";
        }
        cout << endl;
        cout << "want to add ";
        for (auto j : dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids)
        {
            cout << j << " ";
        }
        cout << endl;*/
        rev_node = reverse_hashed_node(dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids);
        if (overlap_backwards(hashed_pg_node_ids, dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids))
        {
            hashed_pg_node_ids.push_front(dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids[0]);
        } else if (overlap_backwards(hashed_pg_node_ids, rev_node))
        {
            hashed_pg_node_ids.push_front(rc_num(dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids.back()));
        } else if (overlap_forwards(hashed_pg_node_ids, dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids))
        {
            hashed_pg_node_ids.push_back(dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids.back());
        } else if (overlap_forwards(hashed_pg_node_ids, rev_node ))
        {
            hashed_pg_node_ids.push_back(rc_num(dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids[0]));
        } else {
            cout << "ERROR" << endl;
        }

        //cout << "length is now " << hashed_pg_node_ids.size() << endl;
    }
    hashed_node_ids_to_ids_and_orientations(hashed_pg_node_ids, node_ids, node_orients);
}

debruijn::Graph construct_debruijn_graph_from_pangraph(uint8_t size, const pangenome::Graph* pg)
{
    cout << now() << "Construct de Bruijn Graph from PanGraph with size " << (uint)size << endl;
    debruijn::Graph dbg(size);

    debruijn::NodePtr prev, current;
    deque<uint16_t> hashed_ids;

    for (auto r : pg->reads)
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

void remove_leaves(pangenome::Graph* pg, debruijn::Graph & dbg, uint16_t covg_thresh)
{
    cout << now() << "Remove leaves of debruijn graph from pangraph" << endl;
    cout << "Start with " << pg->nodes.size() << " pg.nodes, " << pg->reads.size() << " pg.reads, and " << dbg.nodes.size() << " dbg.nodes" << endl;
    bool leaves_exist = true;
    unordered_set<uint32_t> leaves;
    vector<uint16_t> node_ids;
    vector<bool> node_orients;
    uint pos;
    pangenome::NodePtr node;

    while (leaves_exist) {

        leaves = dbg.get_leaves(covg_thresh);
        cout << "there are " << leaves.size() << " leaves" << endl;

        if (leaves.empty())
        {
            leaves_exist = false;
        }
        cout << "leaves exist is " << leaves_exist << endl;

        for (auto i : leaves) {
            cout << endl << "looking at leaf " << i << ": ";
            for (auto j : dbg.nodes[i]->hashed_node_ids)
            {
                cout << j << " ";
            }
            cout << endl;

            // look up the node ids and orientations associated with this node
            hashed_node_ids_to_ids_and_orientations(dbg.nodes[i]->hashed_node_ids, node_ids, node_orients);
	    cout << "looked up node ids" << endl;

            // remove the last node from corresponding reads
            assert(dbg.nodes[i]->read_ids.size()>0);
            for (auto r : dbg.nodes[i]->read_ids)
            {
		cout << "remove from read " << r << ": ";
		for (auto n : pg->reads[r]->nodes)
		{
		    cout << n->node_id << " ";
		}
		cout << endl;
                if (pg->reads[r]->nodes.size() == dbg.size)
                {
		    cout << "remove read";
                    pg->remove_read(r);
		    cout << " done" << endl;
                } else {
		    cout << "remove from read ";
                    pos = pg->reads[r]->find_position(node_ids, node_orients);
		    cout << "pos " << pos;
                    assert(pos == 0 or pos + node_ids.size() == pg->reads[r]->nodes.size());
                    if (pos == 0) {
                        node = pg->reads[r]->nodes[0];
                        pg->reads[r]->remove_node(pg->reads[r]->nodes.begin());
                        node->remove_read(pg->reads[r]);
                    } else if (pos + node_ids.size() == pg->reads[r]->nodes.size()) {
                        node = pg->reads[r]->nodes.back();
                        pg->reads[r]->remove_node(--pg->reads[r]->nodes.end());
                        node->remove_read(pg->reads[r]);
                    }
		    cout << "done" << endl;
                }
            }
            if (node and node->covg == 0)
            {
                pg->remove_node(node);
            }

            // remove dbg node
            dbg.remove_node(i);

            //cout << "pg is now: " << endl << pg << endl;
        }
    }
    cout << "There are now " << pg->nodes.size() << " pg.nodes, " << pg->reads.size() << " pg.reads, and " << dbg.nodes.size() << " dbg.nodes" << endl;
}

void find_reads_along_tig(const debruijn::Graph & dbg,
                     deque<uint32_t>& dbg_node_ids,
                     const pangenome::Graph* pg,
                     vector<uint16_t>& pg_node_ids,
                     vector<bool>& pg_node_orients,
                     unordered_set<pangenome::ReadPtr>& reads_along_tig,
                     bool& all_reads_along_tig)
{
    // collect the reads covering that tig
    for (auto n : dbg_node_ids)
    {
        for (auto r : dbg.nodes.at(n)->read_ids)
        {
            reads_along_tig.insert(pg->reads.at(r));
        }
    }
    /*cout << "candidate reads ";
    for (auto r : reads_along_tig)
    {
        cout << r->id << " ";
    }
    cout << endl;*/

    // filter out some which don't really overlap the unitig, keeping those
    // which overlap at least consecutive 2 dbg nodes or only one node
    all_reads_along_tig = true;
    //cout << "kept reads along tig: ";
    for (unordered_set<pangenome::ReadPtr>::iterator r=reads_along_tig.begin(); r!=reads_along_tig.end();)
    {
        if ((*r)->nodes.size() > dbg.size and (*r)->find_position(pg_node_ids, pg_node_orients, dbg.size+1) == std::numeric_limits<uint>::max())
        {
            r = reads_along_tig.erase(r);
            all_reads_along_tig = false;
        } else {
            //cout << (*r)->id << " ";
            ++r;
        }
    }
    //cout << endl;
}

// Remove the internal nodes of low coverage unitigs e.g.
// suppose when dbg kmer size is 3, we have a low covg tig
// 012 -> 126 -> 263 -> 634 -> 345
// then we would remove the 3 internal kmers from the dbg
// and node 6 from the pg->
// If the tig is smaller than k+2 long, currently does nothing
void filter_unitigs(pangenome::Graph* pg, debruijn::Graph & dbg, const uint16_t& threshold)
{
    cout << now() << "Filter unitigs using threshold " << threshold << endl;
    vector<uint16_t> node_ids;
    vector<bool> node_orients;
    unordered_set<pangenome::ReadPtr> reads_along_tig;
    uint pos;
    bool all_reads_tig;

    set<deque<uint32_t>> unitigs = dbg.get_unitigs();
    cout << "have " << unitigs.size() << " tigs" << endl;
    for (auto d : unitigs)
    {
        // look up the node ids and orientations associated with this node
        dbg_node_ids_to_ids_and_orientations(dbg, d, node_ids, node_orients);
        cout << "tig: ";
        for (auto n : node_ids)
        {
            cout << n << " ";
        }
        cout << endl;

        // collect the reads covering that tig
        find_reads_along_tig(dbg, d, pg, node_ids, node_orients, reads_along_tig, all_reads_tig);

        // now if the number of reads covering tig falls below threshold, remove the
        // middle nodes of this tig from the reads
        if (reads_along_tig.size() <= threshold)
        {
            cout << "not enough reads, so remove the tig from the reads" << endl;
            for (auto r : reads_along_tig)
            {
                cout << "read " << r->id << " was ";
                for (auto n : r->nodes)
                {
                    cout << n->node_id << " ";
                }
                cout << endl;
                pos = r->find_position(node_ids, node_orients);
                //cout << "remove from " << pos+dbg.size << " to max of this and " << (int)(pos+d.size()-dbg.size+2) << " which is " << pos+max((int)dbg.size,(int)d.size()-dbg.size+2) << "." << endl;
                for (auto it = r->nodes.begin()+(pos+dbg.size); it != r->nodes.begin()+pos+max((int)dbg.size,(int)d.size()-dbg.size+2); ++it)
                {
                    if (it == r->nodes.end())
                    {
                        break;
                    }
                    if ((*it)->reads.size()==1 and (*it)->reads.find(r)!=(*it)->reads.end())
                    {
                        pg->remove_node(*it);
                    } else {
                        r->remove_node(it);
                    }
                }
                cout << "read " << r->id << " is now ";
                for (auto n : r->nodes)
                {
                    cout << n->node_id << " ";
                }
                cout << endl;
            }
            // also remove read_ids from each of the corresponding nodes of dbg
            //cout << "now remove read from dbg nodes" << endl;
            for (uint i=1; i<d.size()-1; ++i)
            {
                for (auto r :reads_along_tig)
                {
                    dbg.remove_read_from_node(r->id, d[i]);
                }
            }
        //} else {
        //    cout << "tig had enough reads" << endl;
        }
        reads_along_tig.clear();
    }
}

void detangle_pangraph_with_debruijn_graph(pangenome::Graph* pg, debruijn::Graph & dbg)
{
    cout << now() << "Detangle pangraph with debruijn graph" << endl;
    vector<uint16_t> node_ids;
    vector<bool> node_orients;
    bool all_reads_tig;
    unordered_set<pangenome::ReadPtr> reads_along_tig;

    set<deque<uint32_t>> unitigs = dbg.get_unitigs();
    for (auto d : unitigs)
    {
        // look up the node ids and orientations associated with this node
        dbg_node_ids_to_ids_and_orientations(dbg, d, node_ids, node_orients);
        /*cout << "tig: ";
        for (auto n : node_ids)
        {
            cout << n << " ";
        }
        cout << endl;*/

        // collect the reads covering that tig
        find_reads_along_tig(dbg, d, pg, node_ids, node_orients, reads_along_tig, all_reads_tig);

        // for each node on tig, for each read covering that node,
        // if we find a read which doesn't lie along whole tig,
        // split that node by reads and create a new node on the tig
        if (!all_reads_tig and !reads_along_tig.empty())
        {
            //cout << "not all reads contain tig" << endl;
            for (uint i=0; i<node_ids.size(); ++i)
            {
                //cout << "for node " << pg->nodes[node_ids[i]]->node_id << endl;
                for (auto r : pg->nodes[node_ids[i]]->reads)
                {
                    //cout << "for read " << r->id << endl;
                    if (reads_along_tig.find(r) == reads_along_tig.end())
                    {
                        //cout << "split node" << endl;
                        pg->split_node_by_reads(reads_along_tig, node_ids, node_orients, node_ids[i]);
                        //cout << "done" << endl;
                        break;
                    }
                }
                //cout << pg << endl;
            }
        }
    }
}

void clean_pangraph_with_debruijn_graph(pangenome::Graph* pg, const uint16_t size, const uint16_t threshold)
{
    debruijn::Graph dbg = construct_debruijn_graph_from_pangraph(size,pg);
    remove_leaves(pg, dbg, threshold);
    filter_unitigs(pg, dbg, threshold);

    // update dbg now that have removed leaves and some inner nodes
    dbg = construct_debruijn_graph_from_pangraph(size,pg);

    detangle_pangraph_with_debruijn_graph(pg,dbg);
}
