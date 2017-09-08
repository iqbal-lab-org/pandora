#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include "panread.h"
#include "panedge.h"
#include "pannode.h"
#include "minihits.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

PanRead::PanRead (const uint32_t i): id(i) {}

void PanRead::add_hits(const uint32_t node_id, const set<MinimizerHit*, pComp>& c)
{
    /*map<uint32_t, std::set<MinimizerHit*, pComp>>::iterator it=hits.find(node_id);
    if(it==hits.end())
    {
	hits[node_id] = {};
    }*/
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
