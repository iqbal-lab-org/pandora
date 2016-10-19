#include <iostream>
#include <string>
#include <fstream>
#include<vector>
#include "pannode.h"
#include "seq.h"
//#include "minimizerhit.h"
//#include "utils.h"

using namespace std;

PanNode::PanNode (uint32_t i): id(i) {}

void PanNode::add_read(uint32_t j)
{
    foundReads.push_back(j);
}

/*void Node::add_hits(set<MinimizerHit*> c)
{
    foundHits.insert(c.begin(), c.end());
}*/

bool PanNode::operator == (const PanNode& y) const {
    return (id==y.id);
}

std::ostream& operator<< (std::ostream & out, PanNode const& m) {
    out << m.id;
    return out ;
}
