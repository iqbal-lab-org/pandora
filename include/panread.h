#ifndef __PANREAD_H_INCLUDED__   // if panread.h hasn't been included yet...
#define __PANREAD_H_INCLUDED__

#include <vector>
#include <unordered_map>
#include "minihits.h"

class PanEdge;

class PanRead {
  public:
    const uint32_t id; // corresponding the the read id
    std::vector<PanEdge*> edges;
    std::unordered_map<uint32_t,std::set<MinimizerHit*, pComp_path>> hits; // from node id to cluster of hits against that node in this read

    PanRead(const uint32_t);
    void add_hits(const uint32_t, const std::set<MinimizerHit*, pComp>&);

    bool operator == (const PanRead& y) const;
    bool operator != (const PanRead& y) const;
    bool operator < (const PanRead& y) const;
    friend std::ostream& operator<< (std::ostream& out, const PanRead& r);
};

#endif
