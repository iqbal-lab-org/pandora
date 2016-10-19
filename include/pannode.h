#ifndef __PANNODE_H_INCLUDED__   // if pannode.h hasn't been included yet...
#define __PANNODE_H_INCLUDED__

#include <string>
#include <vector>
//#include "minimizerhit.h"
#include "seq.h"

using namespace std;
using std::vector;


class PanNode {
  public:
    uint32_t id; // corresponding the the LocalPRG id
    vector<PanNode*> outNodes; // representing edges from this node to the nodes in the vector
    vector<uint32_t> foundReads; // representing read ids for those reads intersecting this node
    //set<MinimizerHit*> foundHits;
    PanNode(uint32_t);

    void add_read(uint32_t);
    //void add_hits(set<MinimizerHit*>);
    bool operator == (const PanNode& y) const;
    friend ostream& operator<< (ostream& out, const PanNode& m);
};

#endif
