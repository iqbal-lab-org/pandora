#ifndef __PANNODE_H_INCLUDED__   // if pannode.h hasn't been included yet...
#define __PANNODE_H_INCLUDED__

#include <string>
#include <vector>
#include "minihits.h"
#include "seq.h"

class PanNode {
  public:
    uint32_t id; // corresponding the the LocalPRG id
    std::vector<PanNode*> outNodes; // representing edges from this node to the nodes in the vector
    std::vector<uint32_t> foundReads; // representing read ids for those reads intersecting this node
    std::set<MinimizerHit*, pComp_path> foundHits;
    PanNode(const uint32_t);

    void add_read(const uint32_t);
    void add_hits(const std::set<MinimizerHit*, pComp>&);
    bool operator == (const PanNode& y) const;
    friend std::ostream& operator<< (std::ostream& out, const PanNode& m);
};

#endif
