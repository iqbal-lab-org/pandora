#include <iostream>
#include <fstream>
#include "panedge.h"

using namespace std;

PanEdge::PanEdge (const PanNode* f, const PanNode* t, uint i): from(f), to(t), orientation(i), covg(1) {}

bool PanEdge::operator == (const PanEdge& y) const {
    if (from->id!= y.from->id) {return false;}
    if (to->id!= y.to->id) {return false;}
    if (orientation!= y.orientation) {return false;}
    return true;
}

std::ostream& operator<< (std::ostream & out, PanEdge const& e) {
    out << e.from->id << "->" << e.to->id << " " << e.orientation;
    return out ;
}
