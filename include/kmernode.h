#ifndef __KMERNODE_H_INCLUDED__   // if kmernode.h hasn't been included yet...
#define __KMERNODE_H_INCLUDED__

#include <vector>
#include <ostream>
#include "path.h"

class KmerNode {
    uint32_t id;
    Path path;
    std::vector<KmerNode*> outNodes; // representing edges from this node to the nodes in the vector
    std::vector<KmerNode*> inNodes; // representing edges from other nodes to this node
    uint32_t covg; // covg by hits

  public:
    KmerNode(uint32_t, const Path&);
    bool operator == (const KmerNode& y) const;

  friend std::ostream& operator<< (std::ostream& out, const KmerNode& n);  
  friend class KmerGraph;
  friend struct condition;
};
#endif
