#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <unordered_set>
#include "pansample.h"
#include "panedge.h"
#include "pannode.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

PanSample::PanSample (const string& s): name(s) {}

void PanSample::add_path(const uint32_t node_id, const vector<KmerNodePtr>& c)
{
    if (paths.find(node_id) == paths.end())
    {
	paths[node_id] = {c};
    } else {
	paths[node_id].push_back(c);
    }
}

bool PanSample::operator == (const PanSample& y) const {
    if (name != y.name) {return false;}
    return true;
}

bool PanSample::operator != (const PanSample& y) const {
    return !(*this == y);
}

bool PanSample::operator < (const PanSample& y) const {
    return (name < y.name);
}

std::ostream& operator<< (std::ostream & out, PanSample const& s) {
    out << s.name << ":\t";
    for(auto p : s.paths)
    {
	for (uint i=0; i!=p.second.size(); ++i)
	{
	    out << p.first << "\t";
	}
    }
    return out ;
}
