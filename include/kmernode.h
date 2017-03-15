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
    std::vector<uint32_t> covg; // covg by hits in fwd, rev dir

  public:
    KmerNode(uint32_t, const Path&);
    bool operator == (const KmerNode& y) const;

  friend std::ostream& operator<< (std::ostream& out, const KmerNode& n);  
  friend class KmerGraph;
  friend struct condition;
  friend class LocalPRG;
  friend int pandora_check_kmergraph(int argc, char *argv[]);
  friend class KmerGraphTest_addNode_Test;
  friend class KmerGraphTest_addEdge_Test;
  friend class KmerGraphTest_save_Test;
  friend class KmerGraphTest_load_Test;
  friend class KmerNodeTest_create_Test;
  friend class KmerNodeTest_equals_Test;
  friend class KmerGraphTest_findMaxPathSimple_Test;
  friend class KmerGraphTest_findMaxPath2Level_Test;
};
#endif
