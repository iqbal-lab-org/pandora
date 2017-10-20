#include "tuple.h"
#include "panedge.h"

using namespace std;

Tuple::Tuple(const vector<PanEdge*>& n, PanRead* r) : edges(n), reads({r}) {}

bool Tuple::operator==(const Tuple &y) const 
{
    if (edges.size() != y.edges.size()) {
        return false;
    }
    for (uint i=0; i!=edges.size(); ++i)
    {
	if (edges[i] != y.edges[i] and edges[i] != y.edges[edges.size()-i])
	{
	    return false;
	}
    }
    return true;
}

bool Tuple::operator!=(const Tuple &y) const {
    return !(*this == y);
}

std::ostream &operator<<(std::ostream &out, Tuple const &t) {
    out << "(";
    for (const auto n : t.edges)
    { 
	out << *n << ", ";
    }
    out << ")";
    return out;
}

