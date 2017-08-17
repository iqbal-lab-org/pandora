#ifndef __PANGRAPH_H_INCLUDED__   // if pangraph.h hasn't been included yet...
#define __PANGRAPH_H_INCLUDED__

class PanNode;
struct MinimizerHit;

#include <cstring>
#include <map>
#include <ostream>
#include <functional>
#include <minihits.h>

class PanGraph {
  public:
    std::map<uint32_t, PanNode*> nodes; // representing nodes in graph
    PanGraph() {};
    ~PanGraph();
    void add_node (const uint32_t, const std::string, uint32_t, const std::set<MinimizerHit*, pComp>&);
    void add_edge (const uint32_t&, const uint32_t&, const uint&, const uint&);
    uint rev_orient(const uint&);
    void check_graph_symmetry();
    void remove_low_covg_nodes(const uint&);
    void remove_low_covg_edges(const uint&);
    void remove_isolated_nodes();
    void delete_edge(PanNode*, PanNode*, const uint&);
    void clean(const uint32_t&);
    bool operator == (const PanGraph& y) const;
    void write_gfa (const std::string&);
    friend std::ostream& operator<< (std::ostream& out, const PanGraph& m);
};

#endif
