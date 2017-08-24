#include <iostream>
#include <string>
#include <fstream>
#include "pannode.h"
#include "utils.h"

using namespace std;

PanNode::PanNode (const uint32_t i, const string n): id(i), name(n), covg(1) {}

bool PanNode::operator == (const PanNode& y) const {
    return (id == y.id);
}

bool PanNode::operator != (const PanNode& y) const {
    return (id != y.id);
}

bool PanNode::operator < (const PanNode& y) const {
    return (id < y.id);
}

std::ostream& operator<< (std::ostream & out, PanNode const& n) {
    out << n.id;
    return out ;
}
