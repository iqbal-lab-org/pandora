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
    //cout << "added edge " << *e << endl;

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
PanEdge* PanGraph::add_shortcut_edge(const vector<PanEdge*>::iterator prev, const vector<PanEdge*>::iterator current)
{
    PanEdge* e;

    //have edges A->B and B->C, create edge A->C
    if ((*prev)->to->node_id == (*current)->from->node_id and (*prev)->from->node_id != (*current)->to->node_id)
    {
        e = add_edge((*prev)->from->node_id, (*current)->to->node_id, combine_orientations((*prev)->orientation, (*current)->orientation));
        //cout << "decrease covg on node " << (*prev)->to->node_id << endl;
        (*prev)->to->covg -= 1;
    } else if ((*prev)->to->node_id == (*current)->to->node_id and (*prev)->from->node_id != (*current)->from->node_id)
    {
        e = add_edge((*prev)->from->node_id, (*current)->from->node_id, combine_orientations((*prev)->orientation, rev_orient((*current)->orientation)));
        //cout << "decrease covg on node " << (*prev)->to->node_id << endl;
        (*prev)->to->covg -= 1;
    } else if ((*prev)->from->node_id == (*current)->to->node_id and (*prev)->to->node_id != (*current)->from->node_id)
    {
        e = add_edge((*prev)->to->node_id, (*current)->from->node_id, combine_orientations(rev_orient((*prev)->orientation), rev_orient((*current)->orientation)));
        //cout << "decrease covg on node " << (*prev)->from->node_id << endl;
        (*prev)->from->covg -= 1;
    } else if ((*prev)->from->node_id == (*current)->from->node_id and (*prev)->to->node_id != (*current)->to->node_id)
    {
        e = add_edge((*prev)->to->node_id, (*current)->to->node_id, combine_orientations(rev_orient((*prev)->orientation), (*current)->orientation));
        //cout << "decrease covg on node " << (*prev)->from->node_id << endl;
        (*prev)->from->covg -= 1;
    } else {
	e = nullptr;
    }

    return e;
}

void PanGraph::read_clean(const uint& thresh)
{
    cout << now() << "Start read cleaning with threshold " << thresh << " and " << edges.size() << " edges" << endl;
    PanEdge* e;

    for(map<uint32_t, PanRead*>::iterator read=reads.begin(); read!=reads.end(); ++read)
    {
	cout << "read " << read->first << endl;
        if (read->second->edges.size() == 0)
	{
	    continue;
	} else if (read->second->edges.size() == 1)
	{
            // read has too few nodes, but in the special case where has 1 edge and 2 nodes can 
	    // remove the edge if it has too low coverage
	    if (read->second->edges[0]->covg <= thresh)
	    {
		read->second->edges.pop_back();
	    }	
	    continue;
	}
	vector<PanEdge*>::iterator prev = read->second->edges.begin();
	for (vector<PanEdge*>::iterator current=++read->second->edges.begin(); current!=read->second->edges.end();)
	{
	    //cout << "consider edges" << **prev << " and " << **current << endl;
	    if ((*prev)->covg <= thresh and (*current)->covg <= thresh)
	    {
		cout << "edges " << **prev << " and " << **current << " have low covg";
		e = add_shortcut_edge(prev, current);

		if (e != nullptr)
		{
		    (*current)->covg -= 1;
                    current = read->second->edges.erase(current);

		    (*prev)->covg -= 1;
		    prev = read->second->edges.erase(prev);
		    prev = read->second->edges.insert(prev, e);
		    assert((*prev == e));

		    cout << " so replace with " << *e << endl;
		} else {
		    // this can only happen if we had something circular like A->B and B->A, in which case we want to delete both
		    cout << " but want to avoid self edges, so delete both" << endl;
		    if (current + 1 != read->second->edges.end() and current + 2 != read->second->edges.end())
                    {
			(*current)->covg -= 1;
		        current = read->second->edges.erase(current);
			(*prev)->covg -= 1;
		        prev = read->second->edges.erase(prev);
		        prev = current;
			++current;
		    } else {
			//cout << "current was end so delete and break" << endl;
			(*current)->covg -= 1;
			read->second->edges.erase(current);
			(*prev)->covg -= 1;
			read->second->edges.erase(prev);
			break;
		    }
		}
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
            cout << "delete node " << it->second->name;
	    for(uint i=0; i < it->second->edges.size(); i++)
 	    {
                auto iter = std::find(edges.begin(),edges.end(),it->second->edges[i]);
  		if(iter != edges.end())
  		{
     		    edges.erase(iter);
  		}
 	    }
            delete it->second;
            it = nodes.erase(it);
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
	//cout << "look at edge " << **it << endl;
        if ((*it)->covg <= thresh)
        {
            //cout << "delete edge " << **it << endl;
	    (*it)->from->edges.erase(std::remove((*it)->from->edges.begin(), (*it)->from->edges.end(), (*it)), (*it)->from->edges.end());
	    (*it)->to->edges.erase(std::remove((*it)->to->edges.begin(), (*it)->to->edges.end(), (*it)), (*it)->to->edges.end());
	    delete *it;
            it = edges.erase(it);
        } else {
            ++it;
        }
    }
}

void PanGraph::clean(const uint32_t& coverage)
{
    read_clean(0.025*coverage);
    read_clean(0.05*coverage);
    read_clean(0.1*coverage);
    read_clean(0.2*coverage);

    remove_low_covg_edges(0.025*coverage);
    remove_low_covg_nodes(0);
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
