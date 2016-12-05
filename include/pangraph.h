#ifndef __PANGRAPH_H_INCLUDED__   // if pangraph.h hasn't been included yet...
#define __PANGRAPH_H_INCLUDED__

class PanNode;
//class MinimizerHit;

#include <cstring>
#include <map>
#include <ostream>
#include <functional>

using namespace std;

class PanGraph {
  public:
    map<uint32_t, PanNode*> nodes; // representing nodes in graph
    PanGraph() {};
    ~PanGraph();
    void add_node (uint32_t, uint32_t, const set<MinimizerHit*, pComp>&);
    void add_edge (const uint32_t&, const uint32_t&);
    void write_gfa (const string&);
    bool operator == (const PanGraph& y) const;
    friend ostream& operator<< (ostream& out, const PanGraph& m);
};

#endif
