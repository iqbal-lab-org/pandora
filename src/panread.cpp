#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <unordered_set>
#include "panread.h"
#include "panedge.h"
#include "pannode.h"
#include "minihits.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

PanRead::PanRead (const uint32_t i): id(i) {}

void PanRead::add_hits(const uint32_t node_id, const set<MinimizerHit*, pComp>& c)
{
    hits[node_id].insert(c.begin(), c.end());
}

vector<PanEdge*>::iterator PanRead::get_edge(const PanEdge* e)
{
    return find(edges.begin(), edges.end(), e);
}

vector<PanEdge*>::iterator PanRead::get_next_edge(const PanEdge* e)
{
    vector<PanEdge*>::iterator found = find(edges.begin(), edges.end(), e);
    assert(found != edges.end() || assert_msg("couldn't find edge " << *e << " in read"));
    return ++found;   
}

vector<PanEdge*>::iterator PanRead::get_previous_edge(const PanEdge* e)
{
    vector<PanEdge*>::iterator found = find(edges.begin(), edges.end(), e);
    assert(found != edges.end() || assert_msg("couldn't find edge " << *e << " in read"));
    if (found == edges.begin())
    {
	return edges.end(); // handle the case where e is first edge, and there is nothing before
    }
    return --found;
}

vector<PanEdge*>::iterator PanRead::get_other_edge(const PanEdge* e, const PanNode* n)
{
    //cout << "get other edge containing node " << *n << " which isn't " << *e << endl;
    // get the other edge in read containing node n
    vector<PanEdge*>::iterator found = get_next_edge(e);
    if (found != edges.end() and ((*found)->from == n or (*found)->to == n))
    {
	return found;
    }
    found = get_previous_edge(e);
    if (found != edges.end() and ((*found)->from == n or (*found)->to == n))
    {
        return found;
    }
    return edges.end();
}

vector<PanEdge*>::iterator PanRead::replace_edge(PanEdge* e_original, PanEdge* e, PanRead* me)
{
    vector<PanEdge*>::iterator it = get_edge(e_original);
    if (it != edges.end())
    {   
        it = edges.erase(it);
        it = edges.insert(it, e);
        e->reads.insert(me);
        e->covg += 1;
        e_original->reads.erase(me);
        e_original->covg -= 1;
    }
    return it;
}

unordered_set<PanRead*>::iterator PanRead::replace_edge(PanEdge* e_original, PanEdge* e, unordered_set<PanRead*>::iterator me)
{
    unordered_set<PanRead*>::iterator rt = e_original->reads.find(*me);
    vector<PanEdge*>::iterator it = get_edge(e_original);
    if (it != edges.end())
    {
        it = edges.erase(it);
        it = edges.insert(it, e);
        e->reads.insert(*me);
        e->covg += 1;
	assert(rt!=e_original->reads.end());
        rt = e_original->reads.erase(rt);
        e_original->covg -= 1;
    }
    return rt;
}

vector<PanEdge*>::iterator PanRead::remove_edge(PanEdge* e_original, PanRead* me)
{
    vector<PanEdge*>::iterator it = get_edge(e_original);
    if (it != edges.end())
    {
        it = edges.erase(it);
        e_original->reads.erase(me);
        e_original->covg -= 1;
    }
    return it;
}

unordered_set<PanRead*>::iterator PanRead::remove_edge(PanEdge* e_original, unordered_set<PanRead*>::iterator me)
{
    unordered_set<PanRead*>::iterator rt = e_original->reads.find(*me);
    vector<PanEdge*>::iterator it = get_edge(e_original);
    if (it != edges.end())
    {
        it = edges.erase(it);
        rt = e_original->reads.erase(rt);
        e_original->covg -= 1;
    }
    return rt;
}

void PanRead::replace_node(PanNode* n_original, PanNode* n, PanRead* me)
{
    // does not change the edges including the node in read
    n_original->reads.erase(me);
    n_original->covg -= 1;
    n->reads.insert(me);
    n->covg += 1;
    hits[n->node_id] = hits[n_original->node_id]; //NB we don't remove the n_original hits in case of reads long enough to pass through node twice
} 

void PanRead::remove_node(PanNode* n_original, PanRead* me)
{
    // does not change the edges including the node in read
    n_original->reads.erase(me);
    n_original->covg -= 1;
}

bool PanRead::operator == (const PanRead& y) const {
    if (id != y.id) {return false;}
    /*if (edges.size() != y.edges.size()) {return false;}
    for (uint i=0; i!=edges.size(); ++i)
    {
	if (edges[i] != y.edges[i]) {return false;}
    }*/
	
    return true;
}

bool PanRead::operator != (const PanRead& y) const {
    return !(*this == y);
}

bool PanRead::operator < (const PanRead& y) const {
    return (id < y.id);
}

std::ostream& operator<< (std::ostream & out, PanRead const& r) {
    out << r.id << "\t";
    for (uint i=0; i!=r.edges.size(); ++i)
    {
	out << *r.edges[i] << ", ";
    }
    return out ;
}
