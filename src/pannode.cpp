#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "pannode.h"
//#include "seq.h"
#include "utils.h"

using namespace std;

PanNode::PanNode (const uint32_t i): id(i) {}

void PanNode::add_read(const uint32_t j)
{
    foundReads.push_back(j);
}

void PanNode::add_hits(const set<MinimizerHit*, pComp>& c)
{
    foundHits.insert(c.begin(), c.end());
}

bool PanNode::operator == (const PanNode& y) const {
    if (id!= y.id) {return false;}
    for (uint32_t i=0; i!=outNodes.size(); ++i)
    {
        pointer_values_equal<PanNode> eq = { outNodes[i] };
        if ( find_if(y.outNodes.begin(), y.outNodes.end(), eq) == y.outNodes.end() )
        {//cout << "the out edge points to a different node" << endl; 
            return false;}
    }
    for (uint32_t i=0; i!=y.outNodes.size(); ++i)
    {
        pointer_values_equal<PanNode> eq = { y.outNodes[i] };
        if ( find_if(outNodes.begin(), outNodes.end(), eq) == outNodes.end() )
        {//cout << "the out edge points to a different node" << endl; 
            return false;}
    }
    return true;
}

std::ostream& operator<< (std::ostream & out, PanNode const& m) {
    out << m.id << " -> ";
    for (uint32_t i=0; i!=m.outNodes.size(); ++i)
    { out << m.outNodes[i]->id << " ";}
    return out ;
}
