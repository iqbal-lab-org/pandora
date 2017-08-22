#include <iostream>
#include <string>
#include <fstream>
#include "panread.h"
#include "minihits.h"

using namespace std;

PanRead::PanRead (const uint32_t i): id(i) {}

void PanRead::add_hits(const uint32_t prg_id, const set<MinimizerHit*, pComp>& c)
{
    /*map<uint32_t, std::set<MinimizerHit*, pComp>>::iterator it=hits.find(prg_id);
    if(it==hits.end())
    {
	hits[prg_id] = {};
    }*/
    hits[prg_id].insert(c.begin(), c.end());
}

bool PanRead::operator == (const PanRead& y) const {
    if (id!= y.id) {return false;}
    if (edges.size() != y.edges.size()) {return false;}
    for (uint i=0; i!=edges.size(); ++i)
    {
	if !(edges[i] == y.edges[i]) {return false;}
    }
	
    return true;
}

std::ostream& operator<< (std::ostream & out, PanRead const& r) {
    out << n.id << "\t";
    for (uint i=0; i!=edges.size(); ++i)
    {
	out << edges[i]->from->name << "->"
    }
    return out ;
}
