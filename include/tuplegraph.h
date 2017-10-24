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
    uint32_t next_id;
  public:
    TupleGraph();
    ~TupleGraph();
    Tuple* add_tuple (const std::vector<PanEdge*>&, PanRead*);
    void add_edge (Tuple*, Tuple*);

    bool operator == (const TupleGraph& y) const;
    bool operator != (const TupleGraph& y) const;
    void save(const std::string&);
    friend std::ostream& operator<< (std::ostream& out, const TupleGraph& tg);
};

#endif
