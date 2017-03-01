#ifndef __LOCALNODE_H_INCLUDED__   // if localnode.h hasn't been included yet...
#define __LOCALNODE_H_INCLUDED__

#include <string>
#include <vector>
#include <set>
#include <ostream>
#include <algorithm>
#include "interval.h"
#include "path.h"

class LocalNode {
    std::set<Path> prev_kmer_paths;
  public:
    std::string seq;
    Interval pos;
    uint32_t id;
    uint32_t covg; // covg by hits
    uint32_t sketch_next; // used by minimizer_sketch function in localPRG.cpp

    std::vector<LocalNode*> outNodes; // representing edges from this node to the nodes in the vector

    LocalNode(std::string, Interval, uint32_t);
    bool operator == (const LocalNode& y) const;
  friend std::ostream& operator<< (std::ostream& out, const LocalNode& n);  
  friend class LocalGraph;
  friend class LocalPRG;
};
#endif
