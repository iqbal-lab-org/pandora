#include <iostream>
#include <string>
#include "localnode.h"
#include "interval.h"

using namespace std;

LocalNode::LocalNode (string s, Interval p, uint32_t i): seq(s), pos(p), id(i) {}

std::ostream& operator<< (std::ostream & out, LocalNode const& n) {
    out << "(" << n.id << " " << n.pos << " " << n.seq << ")";
    return out ;
}

