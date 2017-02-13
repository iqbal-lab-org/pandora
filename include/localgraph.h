#ifndef __LOCALGRAPH_H_INCLUDED__   // if localgraph.h hasn't been included yet...
#define __LOCALGRAPH_H_INCLUDED__

class LocalNode;

#include <cstring>
#include <map>
#include <vector>
#include <iostream>
#include "interval.h"
#include "path.h"
#include "localnode.h"

class LocalGraph {
  public:
    std::map<uint32_t, LocalNode*> nodes; // representing nodes in graph
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> index; // varsite index
	// For each nesting level, has a vector of node pairs, each corresponding to a pre-varsite and post-varsite node id. 
	// Every varsite is represented somewhere in this index.
    LocalGraph();
    ~LocalGraph();
    void add_node (const uint32_t& id, const std::string& seq, const Interval& pos);
    void add_edge (const uint32_t&, const uint32_t&);
    void add_varsite (const uint8_t, const uint32_t, const uint32_t);
    void write_gfa (const std::string&);
    void read_gfa (const std::string&);
    std::vector<Path> walk(const uint32_t&, const uint32_t&, const uint32_t&);
    std::vector<LocalNode*> nodes_along_string(const std::string&);
    bool operator == (const LocalGraph& y) const;
    friend std::ostream& operator<< (std::ostream & out, LocalGraph const& data);
};

#endif
