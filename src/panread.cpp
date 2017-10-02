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
    assert(e->from == n or e->to == n || assert_msg("trying to find other edge which contains " << *n << " which isn't " << *e << " but e doesn't belong to n"));
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

vector<PanEdge*>::iterator PanRead::replace_edge(PanEdge* e_original, PanEdge* e)
{
    uint size = edges.size();
    vector<PanEdge*>::iterator it = get_edge(e_original);
    if (it != edges.end())
    {   
        it = edges.erase(it);
        assert(edges.size() == size - 1);
        it = edges.insert(it, e);
        assert(edges.size() == size);
        e->reads.insert(this);
        e->covg += 1;
        e_original->reads.erase(e_original->reads.find(this));
        e_original->covg -= 1;
    }
    return it;
}

unordered_multiset<PanRead*>::iterator PanRead::replace_edge(PanEdge* e_original, PanEdge* e, unordered_multiset<PanRead*>::iterator me)
{
    // Note that the iterator must be over e_original reads
    uint size = e_original->reads.size();
    vector<PanEdge*>::iterator it = get_edge(e_original);
    if (it != edges.end())
    {
        it = edges.erase(it);
        it = edges.insert(it, e);
        e->reads.insert(*me);
        e->covg += 1;
	assert(me!=e_original->reads.end());
        me = e_original->reads.erase(me);
	assert(e_original->reads.size() == size - 1);
        e_original->covg -= 1;
    }
    return me;
}

vector<PanEdge*>::iterator PanRead::remove_edge(PanEdge* e_original)
{
    vector<PanEdge*>::iterator it = get_edge(e_original);
    if (it != edges.end())
    {
        it = edges.erase(it);
        e_original->reads.erase(e_original->reads.find(this));
        e_original->covg -= 1;
    }
    return it;
}

unordered_multiset<PanRead*>::iterator PanRead::remove_edge(PanEdge* e_original, unordered_multiset<PanRead*>::iterator me)
{
    vector<PanEdge*>::iterator it = get_edge(e_original);
    if (it != edges.end())
    {
        it = edges.erase(it);
        me = e_original->reads.erase(me);
        e_original->covg -= 1;
    }
    return me;
    //return rt;
}

void PanRead::replace_node(PanNode* n_original, PanNode* n)
{
    // does not change the edges including the node in read
    if (n_original->reads.find(this)!=n_original->reads.end())
    {
        n_original->reads.erase(n_original->reads.find(this));
        n_original->covg -= 1;
        n->reads.insert(this);
        n->covg += 1;
        hits[n->node_id] = hits[n_original->node_id]; //NB we don't remove the n_original hits in case of reads long enough to pass through node twice
    }
} 

void PanRead::remove_node(PanNode* n_original)
{
    // does not change the edges including the node in read
    if (n_original->reads.find(this)!=n_original->reads.end())
    {
        n_original->reads.erase(n_original->reads.find(this));
        n_original->covg -= 1;
    }
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
