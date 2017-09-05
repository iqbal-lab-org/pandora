#include <iostream>
#include <string>
#include <fstream>
#include "pannode.h"
#include "utils.h"

using namespace std;

PanNode::PanNode (const uint32_t i, const uint32_t j, const string n): prg_id(i), node_id(j), name(n), covg(1) {}

bool PanNode::operator == (const PanNode& y) const {
    return (node_id == y.node_id);
}

bool PanNode::operator != (const PanNode& y) const {
    return (node_id != y.node_id);
}

bool PanNode::operator < (const PanNode& y) const {
    return (node_id < y.node_id);
}

std::ostream& operator<< (std::ostream & out, PanNode const& n) {
    out << n.prg_id << " covg: " << n.covg;
    return out ;
}
