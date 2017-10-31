#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <unordered_set>
#include "pangenome/panread.h"
#include "pangenome/pannode.h"
#include "minihits.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace pangenome;

Read::Read (const uint32_t i): id(i) {}

void Read::add_hits(const uint32_t node_id, const set<MinimizerHitPtr, pComp>& c)
{
    hits[node_id].insert(c.begin(), c.end());
}

void Read::remove_node(NodePtr n_original)
{
    //removes all copies of node
    auto it = find(nodes.begin(), nodes.end(), n_original);
    while (it != nodes.end())
    {
        uint d = distance(nodes.begin(), it);
        nodes.erase(it);
        node_orientations.erase(node_orientations.begin() + d);
        it = find(nodes.begin(), nodes.end(), n_original);
    }
}

/*void Read::replace_node(NodePtr n_original, NodePtr n) //?orientation
{
    if (n_original->reads.find(this)!=n_original->reads.end())
    {
        n_original->reads.erase(n_original->reads.find(this));
        n_original->covg -= 1;
        n->reads.insert(this);
        n->covg += 1;

        hits[n->node_id] = hits[n_original->node_id]; //NB we don't remove the n_original hits in case of reads long enough to pass through node twice
    }
    assert(n_original->covg == n_original->reads.size());
    assert(n->covg == n->reads.size());
}*/

bool Read::operator == (const Read& y) const {
    if (id != y.id) {return false;}

    return true;
}

bool Read::operator != (const Read& y) const {
    return !(*this == y);
}

bool Read::operator < (const Read& y) const {
    return (id < y.id);
}


std::ostream & pangenome::operator<<(std::ostream &out, const pangenome::Read &r) {
    out << r.id << "\t";

    return out;
}

