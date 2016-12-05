#ifndef __LOCALGRAPH_H_INCLUDED__   // if localgraph.h hasn't been included yet...
#define __LOCALGRAPH_H_INCLUDED__

class LocalNode;

#include <cstring>
#include <map>
#include <vector>
#include "interval.h"
#include "path.h"
#include "localnode.h"

using namespace std;

class LocalGraph {
  public:
    map<uint32_t, LocalNode*> nodes; // representing nodes in graph
    LocalGraph() {}
    ~LocalGraph();
    void add_node (const uint32_t& id, const string& seq, const Interval& pos, uint32_t nested_level=0);
    void add_edge (const uint32_t&, const uint32_t&);
    void write_gfa (const string&);
    vector<Path> walk(const uint32_t&, const uint32_t&, const uint32_t&);
    bool operator == (const LocalGraph& y) const;
    void add_read_support_node (LocalNode*);
    void add_read_support_edge (LocalNode*, LocalNode*);
    vector<deque<LocalNode*>> node_step_forwards(vector<deque<LocalNode*>>&);
    vector<deque<LocalNode*>> node_step_back(vector<deque<LocalNode*>>&);
    void infer_read_supported_graph();
    vector<deque<LocalNode*>> get_read_supported_graph_paths(deque<LocalNode*>&);
    void write_read_supported_graph_paths(const string&);
};

#endif
