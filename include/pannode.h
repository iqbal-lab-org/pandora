#ifndef __PANNODE_H_INCLUDED__   // if pannode.h hasn't been included yet...
#define __PANNODE_H_INCLUDED__

#include <string>
#include <vector>
#include "minihits.h"
#include "seq.h"

using namespace std;
using std::vector;

class PanNode {
  public:
    uint32_t id; // corresponding the the LocalPRG id
    vector<PanNode*> outNodes; // representing edges from this node to the nodes in the vector
    vector<uint32_t> foundReads; // representing read ids for those reads intersecting this node
    set<MinimizerHit*, pComp_path> foundHits;
    PanNode(const uint32_t);

    void add_read(const uint32_t);
    void add_hits(const set<MinimizerHit*, pComp>&);
    bool operator == (const PanNode& y) const;
    friend ostream& operator<< (ostream& out, const PanNode& m);
};

#endif
