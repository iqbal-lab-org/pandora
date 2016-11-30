#include <iostream>
#include <string>
#include <algorithm>
#include "localnode.h"
#include "interval.h"
#include "utils.h"

using namespace std;

LocalNode::LocalNode (string s, Interval p, uint32_t i, uint32_t n): seq(s), pos(p), id(i), nested_level(n), covg(0), num_minis(0) {}

std::ostream& operator<< (std::ostream & out, LocalNode const& n) {
    out << "(" << n.id << " " << n.pos << " " << n.seq << ")";
    return out ;
}

/*bool comparator::operator()(LocalNode const* ln)
{
    return *toFine == *ln;
} LocalNode* toFind;*/

bool LocalNode::operator == (const LocalNode& y) const {
    if (seq != y.seq) {//cout << "different seq" << endl; 
	return false;}
    if (!(pos == y.pos)) {//cout << "different interval" << endl; 
	return false;}
    if (id != y.id) {//cout << "different id" << endl; 
	return false;}
    if (outNodes.size() != y.outNodes.size()) {//cout << "differnet numbers of out edges" << endl; 
	return false;}
    for (uint32_t i=0; i!=outNodes.size(); ++i)
    {
        pointer_values_equal<LocalNode> eq = { outNodes[i] };
        if ( find_if(y.outNodes.begin(), y.outNodes.end(), eq) == y.outNodes.end() )
	{//cout << "the out edge points to a different node" << endl; 
	    return false;}
    }
    return true;
}
