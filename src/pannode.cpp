#include <iostream>
#include <string>
#include <fstream>
#include "pannode.h"
#include "utils.h"

using namespace std;

PanNode::PanNode (const uint32_t i, const string n): id(i), name(n), covg(1) {}

bool PanNode::operator == (const PanNode& y) const {
    if (id!= y.id) {return false;}
    return true;
}

std::ostream& operator<< (std::ostream & out, PanNode const& n) {
    out << n.id;
    return out ;
}
