#include <vector>
#include <fstream>
#include <unordered_set>
#include <iostream>
#include "tuplegraph.h"
#include "tuple.h"
#include "utils.h"

using namespace std;

TupleGraph::TupleGraph() : next_id(0) {
    tuples.reserve(100000);
};


TupleGraph::~TupleGraph()
{
    for (auto t : tuples)
    {
	delete t;
    }
    tuples.clear();
}

Tuple* TupleGraph::add_tuple (const vector<PanEdge*>& nodes, PanRead* read)
{
    Tuple* t;
    t = new Tuple(next_id, nodes, read);
    for (auto c : tuples)
    {
	if (*c == *t)
	{
	    c->reads.insert(read);
	    delete t;
	    return c;
	}
    }
    cout << "new tuple " << next_id << endl;
    tuples.insert(t);
    next_id++;
    return t;
}

void TupleGraph::add_edge (Tuple* from, Tuple* to)
{
    if (find(from->outTuples.begin(), from->outTuples.end(), to) == from->outTuples.end())
    {
        from->outTuples.push_back(to);
    }
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

void TupleGraph::save(const string& filepath)
{
    cout << now() << "Save TupleGraph to " << filepath << endl;
    ofstream handle;
    handle.open(filepath);
    handle << *this;
    handle.close();
}

std::ostream &operator<<(std::ostream &out, const TupleGraph& tg) {
    out << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << endl;
    for (const auto t : tg.tuples)
    {
        out << "S\t" << t->id << "\t" << *t << "\n";
	for (const auto s : t->outTuples)
	{
	    out << "L\t" << t->id << "\t+\t" << s->id << "\t+\t0M\n";
	}
    }
    return out;
}
