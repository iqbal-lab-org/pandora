#ifndef __PANNODE_H_INCLUDED__   // if pannode.h hasn't been included yet...
#define __PANNODE_H_INCLUDED__

#include <string>

class PanEdge;
class PanRead;

class PanNode {
  public:
    const uint32_t id; // corresponding the the LocalPRG id
    const std::string name;
    uint32_t covg;

    std::vector<PanEdge*> edges;
    std::unordered_set<PanRead*> reads;

    PanNode(const uint32_t, const std::string);

    bool operator == (const PanNode& y) const;
    friend std::ostream& operator<< (std::ostream& out, const PanNode& n);
};

#endif
