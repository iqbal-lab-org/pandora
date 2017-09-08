#include <iostream>
#include <fstream>
#include <cassert>
#include "panedge.h"
#include "pannode.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

PanEdge::PanEdge (PanNode* f, PanNode* t, uint i): from(f), to(t), orientation(i), covg(1) {
    assert(i<4);
}

/*bool is_out_edge(const PanNode* n)
{
    // is e and out edge of node n?
    assert(from == n or to == n || assert_msg("expected from or to to be " << *n << " but they were " << *from << " and " << *to));
    if (from == n and orientation % 2 == 1)
    {
	return true;
    } else if (to == n and orientation < 2)
    {
	return true;
    }
    return false;
}*/

bool PanEdge::operator == (const PanEdge& y) const {
    if (from->node_id==y.from->node_id and to->node_id==y.to->node_id and orientation==y.orientation) { return true;}
    if (from->node_id==y.to->node_id and to->node_id==y.from->node_id and orientation==rev_orient(y.orientation)) { return true;}
    return false;
}

bool PanEdge::operator != (const PanEdge& y) const {
    return !(*this == y);
}

bool PanEdge::operator < (const PanEdge& y) const {
    if (from->node_id < y.from->node_id) {return true;}
    if (from->node_id > y.from->node_id) {return false;}
    if (to->node_id < y.to->node_id) {return true;}
    if (to->node_id > y.to->node_id) {return false;}
    if (orientation < y.orientation) {return true;}
    return false;
}
std::ostream& operator<< (std::ostream & out, PanEdge const& e) {
    //out << e.from->prg_id << "->" << e.to->prg_id << " " << e.orientation << " covg: " << e.covg;
    out << e.from->node_id << "->" << e.to->node_id << " " << e.orientation << " covg: " << e.covg;
    return out ;
}

uint rev_orient(const uint& orientation)
{
    // 3 A  -> B  = B- -> A- 0
    // 2 A- -> B  = B- -> A  2
    // 0 A- -> B- = B  -> A  3
    // 1 A  -> B- = B  -> A- 1

    uint r_orientation = orientation;
    if (orientation == 0)
    {
        r_orientation = 3;
    } else if (orientation == 3)
    {
        r_orientation = 0;
    }
    return r_orientation;
}

