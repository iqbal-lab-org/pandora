#ifndef __TUPLE_H_INCLUDED__   // if tuple.h hasn't been included yet...
#define __TUPLE_H_INCLUDED__

#include <vector>
#include <unordered_set>
#include <iostream>

class PanRead;
class PanEdge;
class TupleGraph;

class Tuple {
    const uint32_t id;
    std::vector<PanEdge*> edges;
    std::vector<Tuple*> outTuples;
  public:
    std::unordered_multiset<PanRead*> reads;

    Tuple(const uint32_t&, const std::vector<PanEdge*>&, PanRead*);

    bool operator == (const Tuple& y) const;
    bool operator != (const Tuple& y) const;
    friend std::ostream& operator<< (std::ostream& out, const Tuple& t);
    friend std::ostream& operator<< (std::ostream& out, const TupleGraph& tg);
    friend class TupleGraph;
};

#endif
