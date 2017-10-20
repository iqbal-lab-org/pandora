#include <vector>
#include <unordered_set>
#include <iostream>
#include "tuplegraph.h"
#include "tuple.h"
#include "utils.h"

using namespace std;

TupleGraph::TupleGraph(){};

TupleGraph::~TupleGraph()
{
    for (auto t : tuples)
    {
	delete t;
    }
    tuples.clear();
}

void TupleGraph::add_tuple (const vector<PanEdge*>& nodes, PanRead* read)
{
    Tuple* t;
    t = new Tuple(nodes, read);
    for (auto c : tuples)
    {
	if (*c == *t)
	{
	    c->reads.insert(read);
	    delete t;
	    return;
	}
    }
    tuples.insert(t);
}

bool TupleGraph::operator == (const TupleGraph& y) const
{
    if (tuples.size() != y.tuples.size())
    {
	return false;
    }
    for (const auto t : tuples)
    {
	if (y.tuples.find(t) == y.tuples.end())
	{
	    return false;
	}
    }
    for (const auto t : y.tuples)
    {
        if (tuples.find(t) == tuples.end())
        {
            return false;
        }
    }
    return true;
}

bool TupleGraph::operator!=(const TupleGraph &y) const {
    return !(*this == y);
}

