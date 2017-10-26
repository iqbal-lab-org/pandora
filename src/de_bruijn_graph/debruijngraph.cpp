#include <vector>
#include <fstream>
#include <unordered_set>
#include <iostream>
#include "tuplegraph.h"
#include "tuple.h"
#include "utils.h"

using namespace debruijn;

DeBruijnGraph::DeBruijnGraph() : next_id(0) {
    nodes.reserve(20000);
};


DeBruijnGraph::~DeBruijnGraph()
{
    nodes.clear();
}

Node* DeBruijnGraph::add_tuple (const vector<PanEdge*>& nodes, PanRead* read)
{
    Node* t;
    t = new Node(next_id, nodes, read);
    for (auto c : nodes)
    {
	if (*c == *t)
	{
	    c->reads.insert(read);
	    delete t;
	    return c;
	}
    }
    cout << "new tuple " << next_id << endl;
    nodes.insert(t);
    next_id++;
    return t;
}

void DeBruijnGraph::add_edge (Node* from, Node* to)
{
    if (find(from->outNodes.begin(), from->outNodes.end(), to) == from->outNodes.end())
    {
        from->outNodes.push_back(to);
    }
}

bool DeBruijnGraph::operator == (const DeBruijnGraph& y) const
{
    if (nodes.size() != y.nodes.size())
    {
	return false;
    }
    for (const auto t : nodes)
    {
	if (y.nodes.find(t) == y.nodes.end())
	{
	    return false;
	}
    }
    for (const auto t : y.nodes)
    {
        if (nodes.find(t) == nodes.end())
        {
            return false;
        }
    }
    return true;
}

bool DeBruijnGraph::operator!=(const DeBruijnGraph &y) const {
    return !(*this == y);
}

/*void DeBruijnGraph::save(const string& filepath)
{
    cout << now() << "Save DeBruijnGraph to " << filepath << endl;
    ofstream handle;
    handle.open(filepath);
    handle << *this;
    handle.close();
}

std::ostream &operator<<(std::ostream &out, const DeBruijnGraph& tg) {
    out << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << endl;
    for (const auto t : tg.nodes)
    {
        out << "S\t" << t->id << "\t" << *t << "\n";
	for (const auto s : t->outNodes)
	{
	    out << "L\t" << t->id << "\t+\t" << s->id << "\t+\t0M\n";
	}
    }
    return out;
}*/
