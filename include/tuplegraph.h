#ifndef __TUPLEGRAPH_H_INCLUDED__   // if tuplegraph.h hasn't been included yet...
#define __TUPLEGRAPH_H_INCLUDED__

class Tuple;
class PanRead;
class PanEdge;

#include <vector>
#include <unordered_set>
#include <iostream>

class TupleGraph {
    std::unordered_set<Tuple*> tuples; // representing nodes in graph
  public:
    TupleGraph();
    ~TupleGraph();
    void add_tuple (const std::vector<PanEdge*>&, PanRead*);

    bool operator == (const TupleGraph& y) const;
    bool operator != (const TupleGraph& y) const;
    friend std::ostream& operator<< (std::ostream& out, const TupleGraph& tg);
};

#endif
