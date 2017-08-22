#include <iostream>
#include <cstring>
#include <map>
#include <vector>
#include <fstream>
#include "utils.h"
#include "pannode.h"
#include "pangraph.h"
#include <cassert>
#include "minihit.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

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
        n = new PanNode(prg_id, prg_name);
        nodes[prg_id] = n;
    } else {
	n = *it;
	it->covg += 1;
    }

    // add a new read if it doesn't exist
    map<uint32_t, PanRead*>::iterator rit=reads.find(read_id);
    if(rit==reads.end())
	PanRead *r;
        r = new PanRead(read_id);
	r->add_hits(prg_id, cluster);
        reads[read_id] = r;
	n->reads.insert(r);	
    } else {
	rit->add_hits(prg_id, cluster);
	n->reads.insert(*rit);
    }

    return;
}

PanEdge* PanGraph::add_edge (const uint32_t& from, const uint32_t& to, const uint& orientation)
{
    // checks
    map<uint32_t, PanNode*>::iterator from_it=nodes.find(from);
    map<uint32_t, PanNode*>::iterator to_it=nodes.find(to);
    assert((from_it!=nodes.end()) and (to_it!=nodes.end()));
    assert((orientation < (uint)4) || assert_msg("tried to add an edge with a rubbish orientation " << orientation << " which should be < 4"));
    PanNode *f = (nodes.find(from)->second);
    PanNode *t = (nodes.find(to)->second);

    // update edges with new edge or increase edge coverage, and also add edge to read
    PanEdge *e, *rev_e;
    e = new PanEdge(f, t, orientation);
    rev_e = new PanEdge(t, f, rev_orient(orientation));
    pointer_values_equal<PanEdge> eq = { e };
    auto it = find_if(edges.begin(), edges.end(), eq);
    pointer_values_equal<PanEdge> r_eq = { rev_e };
    auto rit = find_if(edges.begin(), edges.end(), r_eq);
    if (it == edges.end() and rit == edges.end())
    {
	edges.push_back(e);
	delete rev_e;
	return e;
    } else if (rit == edges.end()) {
	it->covg += 1;
    	delete e;
        delete rev_e;
	return *it;
    } else {
	rit->covg += 1;
        delete e;
        delete rev_e;
	return *rit;
    }
}

void PanGraph::add_edge (const uint32_t& from, const uint32_t& to, const uint& orientation, const uint& read_id)
{
    PanEdge *e = add_edge(const uint32_t& from, const uint32_t& to, const uint& orientation);
 
    PanRead *r;
    map<uint32_t, PanRead*>::iterator read_it=reads.find(read_id);
    if(read_it==reads.end())
    {
        r = new PanRead(read_id);
        reads[read_id] = r;
    } else {
        r = (read_it-second);
    }   

    r->edges.push_back(e);

    /*PanEdge *e;
    e = new PanEdge(f, t, orientation);
    pointer_values_equal<PanEdge> eq = { e };
    auto it = find_if(edges.begin(), edges.end(), eq);
    if (it != edges.end()) {
        r->edges.push_back(*it);
        delete e;
        delete rev_e;
    } else {
	PanEdge rev_e;
	rev_e = new PanEdge(t, f, rev_orient(orientation));
        pointer_values_equal<PanEdge> r_eq = { rev_e };
	it = find_if(edges.begin(), edges.end(), r_eq);
	assert(it != edges.end());
        r->edges.push_back(*it);
        delete rev_e;
    }
    delete e;*/

}

uint PanGraph::rev_orient(const uint& orientation)
{
    // 3 A  -> B  = B- -> A- 0
    // 2 A- -> B  = B- -> A  2
    // 0 A- -> B- = B  -> A  3
    // 1 A  -> B- = B  -> A- 1

    uint r_orientation = orientation;
    if (orientation == 0)
    {
        r_orientation = 3;
    } else if (orientation == 3)
    {
        r_orientation = 0;
    }
    return r_orientation;
}

void PanGraph::read_clean(const uint& thresh)
{
    PanEdge* e;

    for(map<uint32_t, PanRead*>::iterator read=reads.begin(); read!=reads.end(); ++read)
    {
	vector<PanEdge*>::iterator prev = (*read)->edges.begin();
	for (vector<PanEdge*>::iterator current=++(*read)->edges.begin(); current!=(*read)->edges.end();)
	{
	    if ((*prev)->covg < thresh and (*current)->covg < thresh)
	    {
		//replace these two edge with an new edge
		if ((*prev)->to->id == (*current)->from->id)
		{
		    e = add_edge((*prev)->from->id, (*current)->to->id, (*prev)->orientation + (*current)->orientation -3*((*prev)->orientation > 1));
		    (*prev)->to->covg -= 1;
		} else if ((*prev)->to->id == (*current)->to->id)
		{
		    e = add_edge((*prev)->from->id, (*current)->from->id, (*prev)->orientation + rev_orient((*current)->orientation) -3*((*prev)->orientation > 1));
		    (*prev)->to->covg -= 1;
		} else if ((*prev)->from->id == (*current)->to->id)
                {
                    e = add_edge((*prev)->to->id, (*current)->from->id, rev_orient((*prev)->orientation) + rev_orient((*current)->orientation) -3*(rev_orient((*prev)->orientation) > 1));
		    (*prev)->from->covg -= 1;
                } else if ((*prev)->from->id == (*current)->from->id)
                {
		    e = add_edge((*prev)->to->id, (*current)->to->id, rev_orient((*prev)->orientation) + (*current)->orientation -3*(rev_orient((*prev)->orientation) > 1));
		    (*prev)->from->covg -= 1;
                }
		(*prev)->covg -= 1;
		prev = (*read)->edges.insert(prev, e);
		(*read)->edges.erase(prev+1);
		(*current)->covg -= 1;
		current = (*read)->edges.erase(current);
	    } else {
		prev = current;
		++current;
	    }
	}
    }
}

void PanGraph::remove_low_covg_nodes(const uint& thresh)
{
    cout << now() << "remove nodes with covg <= " << thresh << endl;
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end();)
    {
	if (it->second->covg <= thresh)
        {
            cout << "delete node " << it->second->name;
	    for (uint i=0; i!=it->second->edges.size(); ++i)
            {
		edges.erase(find(edges.begin(), edges.end(), it->second->edges[i]));
	    }
            delete it->second;
            nodes.erase(it++);
            cout << " so pangraph now has " << nodes.size() << " nodes" << endl;
        } else {
            ++it;
        }
    }
}

void PanGraph::remove_low_covg_edges(const uint& thresh)
{
    //WARNING, leaves dead pointers in reads, so can't do read cleaning after
    cout << now() << "remove edges with covg <= " << thresh << endl;
    for(vector<PanEdge*>::iterator it=edges.begin(); it!=edges.end();)
    {
        if ((*it)->covg <= thresh)
        {
            cout << "delete edge " << **it;
            edges.erase(it++);
        } else {
            ++it;
        }
    }
}

void PanGraph::clean(const uint32_t& covg)
{
    uint thresh = 0.025*covg;
    read_clean(thresh);   
    remove_low_covg_edges(thresh);
    remove_low_covg_nodes(0);
}

bool PanGraph::operator == (const PanGraph& y) const {
    // false if have different numbers of nodes
    if (y.nodes.size() != nodes.size()) {
        return false;
    }

    // false if have different nodes
    for ( const auto c: nodes)
    {
        // if node id doesn't exist 
        map<uint32_t, PanNode*>::const_iterator it=y.nodes.find(c.first);
        if(it==y.nodes.end()) {
            return false;
	}
    }
    // or different edges
    for (const auto e: edges)
    {
	pointer_values_equal<PanEdge> eq = { e };
    	auto it = find_if(y.edges.begin(), y.edges.end(), eq);

	PanEdge *re;
	re = new PanEdge(e->to, e->from, rev_orient(e->orientation));
        pointer_values_equal<PanEdge> eq2 = { re };
        auto it2 = find_if(y.edges.begin(), y.edges.end(), eq2);

	if (it == y.edges.end() and it2 == y.edges.end()) {
            return false;
	}

    }
    // otherwise is true
    return true;
}

void PanGraph::write_gfa (const string& filepath)
{
    ofstream handle;
    handle.open (filepath);
    handle << "H\tVN:Z:1.0" << endl;
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
	
        handle << "S\t" << it->second->name << "\t*" << endl; //\tRC:i:" << it->second->foundHits.size() * k << endl;
    }

    for (uint i=0; i!=edges.size(); ++i)
    {
	if (edges[i]->orientation == 0)
	{
	    handle << "L\t" << edges[i]->from->name << "\t-\t" << edges[i]->to->name << "\t-\t0M\tRC:i:" << edges[i]->covg << endl;
	} else if (edges[i]->orientation == 1)
	{
	    handle << "L\t" << edges[i]->from->name << "\t+\t" << edges[i]->to->name << "\t-\t0M\tRC:i:" << edges[i]->covg << endl;
        } else if (edges[i]->orientation == 2)
        {
            handle << "L\t" << edges[i]->from->name << "\t-\t" << edges[i]->to->name << "\t+\t0M\tRC:i:" << edges[i]->covg << endl;
        } else if (edges[i]->orientation == 3)
        {
            handle << "L\t" << edges[i]->from->name << "\t+\t" << edges[i]->to->name << "\t+\t0M\tRC:i:" << edges[i]->covg << endl;
        }
    }
    handle.close();
}

/*std::ostream& operator<< (std::ostream & out, PanGraph const& m) {
    for(map<uint32_t, PanNode*>::const_iterator it=m.nodes.begin(); it!=m.nodes.end(); ++it)
    {
        for (uint32_t j=0; j<it->second->outNodes[3].size(); ++j)
        {
            out << "+" << it->second->id << " -> +" << it->second->outNodes[3][j]->id << endl;
        }
	for (uint32_t j=0; j<it->second->outNodes[1].size(); ++j)
        {
            out << "+" << it->second->id << " -> -" << it->second->outNodes[1][j]->id << endl;
        }
	
    }
    out << endl;
    return out ;
}*/
