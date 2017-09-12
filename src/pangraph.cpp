#include <iostream>
#include <cstring>
#include <map>
#include <vector>
#include <fstream>
#include "utils.h"
#include "pannode.h"
#include "pangraph.h"
#include "panread.h"
#include "panedge.h"
#include <cassert>
#include "minihit.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

PanGraph::PanGraph() : next_id(0) {}

PanGraph::~PanGraph()
{
    for (auto c: reads)
    { 
        delete c.second;
    }
    for (auto c: edges)
    {   
        delete c;
    }
    for (auto c: nodes)
    {
        delete c.second;
    }
}

void PanGraph::add_node (const uint32_t prg_id, const string prg_name, const uint32_t read_id, const set<MinimizerHit*, pComp>& cluster)
{
    // check sensible things in new cluster
    for (set<MinimizerHit*, pComp>::iterator it = cluster.begin(); it != cluster.end(); ++it)
    {
            assert(read_id == (*it)->read_id); // the hits should correspond to the read we are saying...
	    assert(prg_id == (*it)->prg_id);
    }

    // add new node if it doesn't exist
    PanNode *n;
    map<uint32_t, PanNode*>::iterator it=nodes.find(prg_id);
    if(it==nodes.end())
    {
        n = new PanNode(prg_id, prg_id, prg_name);
	//cout << "add node " << *n << endl;
        nodes[prg_id] = n;
    } else {
	n = it->second;
	it->second->covg += 1;
	//cout << "node " << *n << " already existed " << endl;
    }

    // add a new read if it doesn't exist
    map<uint32_t, PanRead*>::iterator rit=reads.find(read_id);
    if (rit==reads.end())
    {
	//cout << "new read " << read_id << endl;
	PanRead* r;
        r = new PanRead(read_id);
	r->add_hits(prg_id, cluster);
        reads[read_id] = r;
	n->reads.insert(r);	
    } else {
	//cout << "read " << read_id  << " already existed " << endl;
	rit->second->add_hits(prg_id, cluster);
	pointer_values_equal<PanRead> eq = { rit->second };
        if (find_if(n->reads.begin(), n->reads.end(), eq) == n->reads.end())
	{
	    n->reads.insert(rit->second);
	}
    }

    return;
}

PanEdge* PanGraph::add_edge (const uint32_t& from, const uint32_t& to, const uint& orientation)
{
    // NB this adds an edge from node_id from to node_id to
    //
    // checks
    map<uint32_t, PanNode*>::iterator from_it=nodes.find(from);
    map<uint32_t, PanNode*>::iterator to_it=nodes.find(to);
    assert((from_it!=nodes.end()) and (to_it!=nodes.end()));
    assert((orientation < (uint)4) || assert_msg("tried to add an edge with a rubbish orientation " << orientation << " which should be < 4"));
    PanNode *f = (nodes.find(from)->second);
    PanNode *t = (nodes.find(to)->second);

    // update edges with new edge or increase edge coverage, and also add edge to read
    PanEdge *e;
    e = new PanEdge(f, t, orientation);
    pointer_values_equal<PanEdge> eq = { e };
    auto it = find_if(edges.begin(), edges.end(), eq);
    if (it == edges.end())
    {
	//cout << "add edge " << *e << endl;
	edges.push_back(e);
	f->edges.push_back(e);
        t->edges.push_back(e);
	return e;
    } else {
	(*it)->covg += 1;
	//cout << "edge " << **it << " already in graph" << endl;
    	delete e;
	return *it;
    }
}

void PanGraph::add_edge (const uint32_t& from, const uint32_t& to, const uint& orientation, const uint& read_id)
{
    PanEdge* e;
    e = add_edge(from, to, orientation);
 
    PanRead *r;
    map<uint32_t, PanRead*>::iterator read_it=reads.find(read_id);
    if(read_it==reads.end())
    {
        r = new PanRead(read_id);
        reads[read_id] = r;
    } else {
        r = (read_it->second);
    }   

    r->edges.push_back(e);
    e->reads.insert(r);
    //cout << *e << " has read size " << e->reads.size() << endl;
    //cout << "added edge " << *e << endl;

}

vector<PanEdge*>::iterator PanGraph::remove_edge (PanEdge* e)
{
    cout << "Remove graph edge " << *e << endl;
    // remove mentions of edge from nodes
    e->from->edges.erase(std::remove(e->from->edges.begin(), e->from->edges.end(), e), e->from->edges.end());
    e->to->edges.erase(std::remove(e->to->edges.begin(), e->to->edges.end(), e), e->to->edges.end());
    // remove mentions from reads. Note that if the edge is present in any reads, the vector of
    // edges along that read will become inconsistent
    for (unordered_set<PanRead*>::iterator r = e->reads.begin(); r!=e->reads.end();)
    {
	cout << "read was " << **r << endl;
        r = (*r)->remove_edge(e, r);
	//cout << "read is now " << **r << endl;
    }

    vector<PanEdge*>::iterator it = find(edges.begin(), edges.end(),e);
    delete e;
    it = edges.erase(it);
    //return iterator to next edge in pangraph
    return it;
}

map<uint32_t, PanNode*>::iterator PanGraph::remove_node(PanNode* n)
{
    cout << "Remove graph node " << *n << endl;
    // removes node n and for consistency all edges including n
    // this will remove those edges from reads which is likely to lead to
    // inconsistent vectors of edges in reads unless handled.
    for(uint i=n->edges.size(); i>0; --i)
    {
	remove_edge(n->edges[i-1]);
    }
    map<uint32_t, PanNode*>::iterator it = nodes.find(n->node_id);
    delete n;
    if (it!=nodes.end())
    {
        it = nodes.erase(it);
    }
    return it;
}

uint combine_orientations(uint f, uint t)
{
    uint nice = f + t -3*(f>1);
    uint fix = (f%2==1) + 2*(t>1);
    if (nice != fix)
    {
	cout << "Warning, the edges were incompatible, but we handled it" << endl;
    }
    return fix;
}
PanEdge* PanGraph::add_shortcut_edge(const vector<PanEdge*>::iterator prev, const vector<PanEdge*>::iterator current, PanRead* r)
{
    PanEdge* e = nullptr;
    PanNode* n = nullptr;

    //have edges A->B and B->C, create edge A->C
    if ((*prev)->to->node_id == (*current)->from->node_id and (*prev)->from->node_id != (*current)->to->node_id)
    {
        e = add_edge((*prev)->from->node_id, (*current)->to->node_id, combine_orientations((*prev)->orientation, (*current)->orientation));
	e->covg -= 1;
        n = (*prev)->to;
    } else if ((*prev)->to->node_id == (*current)->to->node_id and (*prev)->from->node_id != (*current)->from->node_id)
    {
        e = add_edge((*prev)->from->node_id, (*current)->from->node_id, combine_orientations((*prev)->orientation, rev_orient((*current)->orientation)));
        e->covg -= 1;
	n = (*prev)->to;
    } else if ((*prev)->from->node_id == (*current)->to->node_id and (*prev)->to->node_id != (*current)->from->node_id)
    {
        e = add_edge((*prev)->to->node_id, (*current)->from->node_id, combine_orientations(rev_orient((*prev)->orientation), rev_orient((*current)->orientation)));
        e->covg -= 1;
	n = (*prev)->from;
    } else if ((*prev)->from->node_id == (*current)->from->node_id and (*prev)->to->node_id != (*current)->to->node_id)
    {
        e = add_edge((*prev)->to->node_id, (*current)->to->node_id, combine_orientations(rev_orient((*prev)->orientation), (*current)->orientation));
        e->covg -= 1;
	n = (*prev)->from;
    } else {
        vector<PanEdge*>::iterator it = r->get_previous_edge(*prev);
	if (it==r->edges.end())
	{
	    it = r->get_next_edge(*current);
        }
	if (it!=r->edges.end() and ((*it)->from == (*prev)->from or (*it)->to == (*prev)->from))
	{
	    // then to is the node to remove
	    n = (*prev)->to;
	} else if (it!=r->edges.end() and ((*it)->from == (*prev)->to or (*it)->to == (*prev)->to))
	{
	    n = (*prev)->from;
	}
    }

    // if n is not null, remove node from read r
    if (n != nullptr)
    {
	r->remove_node(n, r);
    }

    return e;
}

vector<PanEdge*>::iterator PanGraph::split_node_by_edges(PanNode* n_original, PanEdge* e_original1, PanEdge* e_original2)
{
    // results in two copies of n_original, one keeping the reads covered by e_original1
    // and the other without them.
    cout << "Split node " << *n_original << " by edge " << *e_original1 << " (and " << *e_original2 << ")" << endl;
    assert(n_original == e_original1->from or n_original == e_original1->to);

    // give a unique id
    while( nodes.find(next_id) != nodes.end() )
    {
        next_id++;
    }

    // create new node
    PanNode *n;
    n = new PanNode(n_original->prg_id, next_id, n_original->name);
    n->covg = 0;
    nodes[next_id] = n;

    // create new edges
    PanEdge *e1, *e2;
    //cout << "add edge equivalent to " << *e_original1 << endl;
    if (n_original == e_original1->from)
    {
	//cout << "try to add edge " << n->node_id << " -> " << e_original1->to->node_id << " " << e_original1->orientation << endl;
   	e1 = add_edge(n->node_id, e_original1->to->node_id, e_original1->orientation);
    } else if (n_original == e_original1->to)
    {
	//cout << "try to add edge " << e_original1->from->node_id << " -> " << n->node_id << " " << e_original1->orientation << endl;
	e1 = add_edge(e_original1->from->node_id, n->node_id, e_original1->orientation);
    }
    //cout << "add edge equivalent to " << *e_original2 << endl;
    if (n_original == e_original2->from)
    {
	//cout << "try to add edge " << n->node_id << " -> " << e_original2->to->node_id << " " << e_original2->orientation << endl;   
        e2 = add_edge(n->node_id, e_original2->to->node_id, e_original2->orientation);
    } else if (n_original == e_original2->to)
    {
	//cout << "try to add edge " << e_original2->from->node_id << " -> " << n->node_id << " " << e_original2->orientation << endl;
        e2 = add_edge(e_original2->from->node_id, n->node_id, e_original2->orientation);
    }
    // we are going to add covg for each read, so subtract the covg added by creating the edge
    e1->covg -= 1;
    e2->covg -= 1;
    
    // amend reads
    //cout << "now fix reads" << endl;
    vector<PanEdge*>::iterator it, it2;
    unordered_set<PanRead*>::iterator r = e_original1->reads.begin();
    while (r!=e_original1->reads.end())
    {
	// add read to new node n and remove from n_original
	(*r)->replace_node(n_original, n, *r);

        // replace e_original2 in read
        (*r)->replace_edge(e_original2, e2, r);
	
	// replace e_original1 in read
	r = (*r)->replace_edge(e_original1, e1, r);

    }

    //cout << "delete edges if they have 0 covg" << endl;
    // if e_original1 or e_original2 now have 0 covg, delete both from nodes and graph
    if (e_original2->covg == 0)
    {
        assert(e_original2->reads.size() == 0);
	remove_edge(e_original2);
    }

    //cout << "same for original edge" << endl;
    assert(e_original1->covg == 0);
    assert(e_original1->reads.size() == 0);
    if (e_original1->to == n_original)
    {
	it = find(n_original->edges.begin(), n_original->edges.end(), e_original1);
	it = n_original->edges.erase(it);
    } else {
	it = find(n_original->edges.begin(), n_original->edges.end(), e_original1);
        it = n_original->edges.erase(it);
    }
    remove_edge(e_original1);
    
    return it;	    
}

void PanGraph::split_nodes_by_reads(const uint& thresh)
{
    cout << "Start breaking up pannodes by reads" << endl;
    for (auto n : nodes)
    {
        if (n.second->edges.size() > 2 and n.second->covg > thresh)
        {
	    vector<PanEdge*>::iterator it=n.second->edges.begin();
            while (it!= n.second->edges.end() and n.second->edges.size() > 2)
            {
                // see if all reads along this edge enter/leave the node along the same edge
                unordered_set<PanEdge*> s;
                for (auto r : (*it)->reads)
                {
		    //cout << "look at node " << n.first << " edge " << **it << " and read " << *r << endl;
                    vector<PanEdge*>::iterator found = r->get_other_edge(*it, n.second);
                    if ( found != r->edges.end())
                    {
			//cout << "found edge " << **found << endl;
                        s.insert(*found);
                    }
                    if (s.size() > 1)
                    {
                        break;
                    /*} else {
			cout << "set contains: ";
			for (auto i : s)
			{
			    cout << *i << endl;
			}*/
		    }
                }

                // if they do, clone the node
                if (s.size() == 1)
                {
                    it = split_node_by_edges(n.second, *it, *s.begin());
                } else {
		    ++it;
		}
            }
        }
    }
    cout << "Finished splitting nodes by reads. PanGraph now has  " << nodes.size() << " nodes and " << edges.size() << " edges" << endl;
    return;
}

void PanGraph::read_clean(const uint& thresh)
{
    cout << now() << "Start read cleaning with threshold " << thresh << " and " << edges.size() << " edges" << endl;
    PanEdge* e;

    for(map<uint32_t, PanRead*>::iterator read=reads.begin(); read!=reads.end(); ++read)
    {
        if (read->second->edges.size() < 2)
	{
	    continue;
	}
	vector<PanEdge*>::iterator prev = read->second->edges.begin();
	for (vector<PanEdge*>::iterator current=++read->second->edges.begin(); current!=read->second->edges.end();)
	{
	    //cout << "consider edges" << **prev << " and " << **current << endl;
	    if ((*prev)->covg <= thresh and (*current)->covg <= thresh)
	    {
		cout << "read " << read->first << " edges " << **prev << " and " << **current << " have low covg " << endl;
		e = add_shortcut_edge(prev, current, read->second);

		if (e != nullptr)
		{
		    current = read->second->remove_edge(*current, read->second);
		    prev = read->second->replace_edge(*prev, e, read->second);
		    assert((*prev == e));
		} else {
		    // this can only happen if we had something circular like A->B and B->A, in which case we want to delete both
		    if (current + 1 != read->second->edges.end() and current + 2 != read->second->edges.end())
		    {
			current = read->second->remove_edge(*current, read->second);
			prev = read->second->remove_edge(*prev, read->second);
			prev = current;
                        ++current;
		    } else {
		   	read->second->remove_edge(*current, read->second);
                        read->second->remove_edge(*prev, read->second);
			break;
		    }
		}
		cout << "read is now: " << *(read->second) << endl;
	    } else {
		prev = current;
		++current;
	    }
	}
    }
    cout << now() << "Finished read cleaning. PanGraph now has " << edges.size() << " edges" << endl;
}

void PanGraph::remove_low_covg_nodes(const uint& thresh)
{
    cout << now() << "remove nodes with covg <= " << thresh << endl;
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end();)
    {
	//cout << "look at node " << *(it->second) << endl;
	if (it->second->covg <= thresh or it->second->edges.size() == 0)
        {
            //cout << "delete node " << it->second->name;
            it = remove_node(it->second);
            cout << " so pangraph now has " << nodes.size() << " nodes" << endl;
        } else {
            ++it;
        }
    }
}

void PanGraph::remove_low_covg_edges(const uint& thresh)
{
    cout << now() << "remove edges with covg <= " << thresh << endl;
    for(vector<PanEdge*>::iterator it=edges.begin(); it!=edges.end();)
    {
	//cout << "look at edge " << **it << endl;
        if ((*it)->covg <= thresh)
        {
            //cout << "delete edge " << **it << endl;
            it = remove_edge(*it);
        } else {
            ++it;
        }
    }
}

/*reroute_edge(PanEdge* e, PanEdge* d1, PanEdge* d2)
{
    if (combine_orientations((*out)->orientation, (*to)->orientation) == it->orientation)
    {
    }

void PanGraph::remove_triangles(const uint& thresh)
{
    // removes triangles where A-lo->C and A-hi->B-hi->C
    cout << now() << "remove triangles using threshhold " << thresh << endl;
    for(vector<PanEdge*>::iterator it=edges.begin(); it!=edges.end();)
    {
        //cout << "look at edge " << **it << endl;
        if ((*it)->covg <= thresh)
        {
            //cout << "delete edge " << **it << endl;
	    for(vector<PanEdge*>::iterator out=from->edges.begin(); out!=from->edges.end();)
	    {
		for(vector<PanEdge*>::iterator in=to->edges.begin(); in!=to->edges.end();)
            	{
		    if ((*out)->covg > 2*thresh and (*in)->covg > 2*thresh)
		    {
			if ((*out)->from == (*in)->from)
			{
			    for 
			} else if ((*out)->from == (*in)->to)
                        {
                        } else if ((*out)->to == (*in)->from)
			    reroute_reads(*it, *out, *in)(combine_orientations((*out)->orientation, (*to)->orientation) == it->orientation or )
                        {
			    for (auto r: reads)
			    {
				
			    }
                        } else if ((*out)->to == (*in)->to)
                        {
                        }
		    }
            (*it)->from->edges.erase(std::remove((*it)->from->edges.begin(), (*it)->from->edges.end(), (*it)), (*it)->from->edges.end());
            (*it)->to->edges.erase(std::remove((*it)->to->edges.begin(), (*it)->to->edges.end(), (*it)), (*it)->to->edges.end());
            delete *it;
            it = edges.erase(it);
        } else {
            ++it;
        }
    }
}*/

void PanGraph::clean(const uint32_t& coverage)
{
    // get total edge covg / total node covg (will be < 1 is short reads and closer to 1 if long reads)
    uint edge_covg=0, node_covg=0;
    for (uint i=0; i!=edges.size(); ++i)
    {
	edge_covg += edges[i]->covg;
    }
    cout << edge_covg << endl;
    for (auto n : nodes)
    {
	node_covg += n.second->covg;
    }
    cout << node_covg << endl;
    float edges_per_node = coverage * edge_covg / node_covg;
    cout << "factor " << edges_per_node << endl;
    read_clean(0.025*edges_per_node);
    read_clean(0.05*edges_per_node);
    read_clean(0.1*edges_per_node);
    read_clean(0.2*edges_per_node);
    remove_low_covg_edges(0);
    remove_low_covg_nodes(0.05*coverage);

    split_nodes_by_reads(edges_per_node);
    read_clean(0.2*edges_per_node);

    remove_low_covg_edges(0.05*edges_per_node);
    remove_low_covg_nodes(0.05*coverage);
}

void PanGraph::add_hits_to_kmergraphs(const vector<LocalPRG*>& prgs)
{
    uint num_hits[2];
    for (auto pnode : nodes)
    {
        // copy kmergraph
	pnode.second->kmer_prg = prgs[pnode.second->prg_id]->kmer_prg;
	num_hits[0] = 0;
	num_hits[1] = 0;
    
        // add hits
        for (auto read : pnode.second->reads)
        {
            for (set<MinimizerHit*, pComp_path>::iterator mh = read->hits[pnode.second->node_id].begin(); mh != read->hits[pnode.second->node_id].end(); ++mh)
            {
		//bool added = false;
                // update the covg in the kmer_prg
		pnode.second->kmer_prg.nodes[(*mh)->knode_id]->covg[(*mh)->strand] += 1;
		num_hits[(*mh)->strand] += 1;
                /*for (uint i=0; i!=pnode.second->kmer_prg.nodes.size(); ++i)
                {
                    if (pnode.second->kmer_prg.nodes[i]->path == (*mh)->prg_path)
                    {
                        pnode.second->kmer_prg.nodes[i]->covg[(*mh)->strand] += 1;
			num_hits[(*mh)->strand] += 1;
                        added = true;
                        break;
                    }
                }
                assert(added == true || assert_msg("could not find kmernode corresponding to " << **mh));*/
            }
        }
        cout << now() << "Added " << num_hits[1] << " hits in the forward direction and " << num_hits[0] << " hits in the reverse" << endl;
    }
}

bool PanGraph::operator == (const PanGraph& y) const {
    // false if have different numbers of nodes
    if (y.nodes.size() != nodes.size()) {
	cout << "different num nodes " << nodes.size() << "!=" << y.nodes.size() << endl;
        return false;
    }

    // false if have different numbers of edges
    if (y.edges.size() != edges.size()) {
        cout << "different num edges " << edges.size() << "!=" << y.edges.size() << endl;
        return false;
    }

    // false if have different nodes
    for ( const auto c: nodes)
    {
        // if node id doesn't exist 
        map<uint32_t, PanNode*>::const_iterator it=y.nodes.find(c.first);
        if(it==y.nodes.end()) {
	    cout << "can't find node " << c.first << endl;
            return false;
	}
    }
    // or different edges
    for (const auto e: edges)
    {
	pointer_values_equal<PanEdge> eq = { e };
    	auto it = find_if(y.edges.begin(), y.edges.end(), eq);

	if (it == y.edges.end()) {
	    cout << "couldn't find " << *e << endl;
            return false;
	}

    }
    // otherwise is true
    return true;
}

bool PanGraph::operator != (const PanGraph& y) const {
    return !(*this == y);
}

void PanGraph::write_gfa (const string& filepath)
{
    ofstream handle;
    handle.open (filepath);
    handle << "H\tVN:Z:1.0" << endl;
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
	
        handle << "S\t" << it->second->name << "." << it->first << "\t*\tFC:i:" << it->second->covg << endl;
    }

    for (uint i=0; i!=edges.size(); ++i)
    {
	if (edges[i]->orientation == 0)
	{
	    handle << "L\t" << edges[i]->from->name << "." << edges[i]->from->node_id << "\t-\t" << edges[i]->to->name << "." << edges[i]->to->node_id << "\t-\t0M\tRC:i:" << edges[i]->covg << endl;
	} else if (edges[i]->orientation == 1)
	{
	    handle << "L\t" << edges[i]->from->name << "." << edges[i]->from->node_id << "\t+\t" << edges[i]->to->name << "." << edges[i]->to->node_id << "\t-\t0M\tRC:i:" << edges[i]->covg << endl;
        } else if (edges[i]->orientation == 2)
        {
            handle << "L\t" << edges[i]->from->name << "." << edges[i]->from->node_id << "\t-\t" << edges[i]->to->name << "." << edges[i]->to->node_id << "\t+\t0M\tRC:i:" << edges[i]->covg << endl;
        } else if (edges[i]->orientation == 3)
        {
            handle << "L\t" << edges[i]->from->name << "." << edges[i]->from->node_id << "\t+\t" << edges[i]->to->name << "." << edges[i]->to->node_id << "\t+\t0M\tRC:i:" << edges[i]->covg << endl;
        }
    }
    handle.close();
}

std::ostream& operator<< (std::ostream & out, PanGraph const& m) {
    //cout << "printing pangraph" << endl;
    for (auto n : m.nodes)
    {
	cout << n.second->prg_id << endl;
    }
    for (auto e : m.edges)
    {
        cout << *e << endl;
    }
    return out ;
}
