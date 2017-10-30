#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <unordered_set>
#include "pangenome/pansample.h"
#include "pangenome/pannode.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace pangenome;

Sample::Sample (const string& s): name(s) {}

void Sample::add_path(const uint32_t node_id, const vector<KmerNodePtr>& c)
{
    if (paths.find(node_id) == paths.end())
    {
	    paths[node_id] = {c};
    } else {
	    paths[node_id].push_back(c);
    }
}

bool Sample::operator == (const Sample& y) const {
    if (name != y.name) {return false;}
    return true;
}

bool Sample::operator != (const Sample& y) const {
    return !(*this == y);
}

bool Sample::operator < (const Sample& y) const {
    return (name < y.name);
}

std::ostream& pangenome::operator<< (std::ostream & out, pangenome::Sample const& s) {
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
