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
    uint64_t khash; //the kmer hash value
    uint8_t num_AT; // the number of As and Ts in this kmer

  public:
    KmerNode(uint32_t, const Path&);
    bool operator == (const KmerNode& y) const;

  friend bool equal_except_null_nodes (const KmerNode& x, const KmerNode& y);
  friend std::ostream& operator<< (std::ostream& out, const KmerNode& n);  
  friend class KmerGraph;
  friend struct condition;
  friend struct pCompKmerNode;
  friend class LocalPRG;
  friend int pandora_check_kmergraph(int argc, char *argv[]);
  friend class KmerGraphTest_addNode_Test;
  friend class KmerGraphTest_addEdge_Test;
  friend class KmerGraphTest_save_Test;
  friend class KmerGraphTest_load_Test;
  friend class KmerNodeTest_create_Test;
  friend class KmerNodeTest_equals_Test;
  friend class KmerGraphTest_sortTopologically_Test;
  friend class KmerGraphTest_findMaxPathSimple_Test;
  friend class KmerGraphTest_findMaxPath2Level_Test;
  friend class LocalPRGTest_minimizerSketch_Test;
  friend class LocalPRGTest_minimizerSketchSameAsSeqw1_Test;
  friend class LocalPRGTest_minimizerSketchSameAsSeqw5_Test;
  friend class LocalPRGTest_minimizerSketchSameAsSeqw10_Test;
  friend class LocalPRGTest_minimizerSketchSameAsSeqw15_Test;
  friend class LocalPRGTest_updateCovgWithHit_Test;
  friend class LocalPRGTest_find_path_and_variants_Test;
};

#endif
