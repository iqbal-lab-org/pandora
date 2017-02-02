#ifndef __LOCALNODE_H_INCLUDED__   // if localnode.h hasn't been included yet...
#define __LOCALNODE_H_INCLUDED__

#include <string>
#include <vector>
#include <ostream>
#include <algorithm>
#include "interval.h"

class LocalNode {
  public:
    std::string seq;
    Interval pos;
    uint32_t id;
    //uint32_t nested_level;
    uint32_t covg; // covg by hits
    uint32_t sketch_next; // used by minimizer_sketch function in localPRG.cpp
    //uint32_t num_minis; // number of minimizers covering node in sketch of PRG
    //bool supported_by_reads;
    //uint32_t supported_in_deg;
    //uint32_t supported_out_deg;
    std::vector<LocalNode*> outNodes; // representing edges from this node to the nodes in the vector
    //vector<LocalNode*> supported_outNodes; // representing edges from this node to another supportd nodes in the vector
    LocalNode(std::string, Interval, uint32_t);
    bool operator == (const LocalNode& y) const;
  friend std::ostream& operator<< (std::ostream& out, const LocalNode& n);  
  friend class LocalGraph;
};
#endif
