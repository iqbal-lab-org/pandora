#ifndef __PANSAMPLE_H_INCLUDED__   // if pansample.h hasn't been included yet...
#define __PANSAMPLE_H_INCLUDED__

#include <vector>
#include <unordered_map>

class PanEdge;
class KmerNode;

class PanSample {
  public:
    const std::string name; // corresponding the the read id
    std::vector<PanEdge*> edges;
    std::unordered_map<uint32_t, std::vector<std::vector<KmerNode*>>> paths; // from prg id (or unique id) to kmernnode path(s) through each node

    PanSample(const std::string&);
    void add_path(const uint32_t, const std::vector<KmerNode*>&);

    bool operator == (const PanSample& y) const;
    bool operator != (const PanSample& y) const;
    bool operator < (const PanSample& y) const;
    friend std::ostream& operator<< (std::ostream& out, const PanSample& r);
};

#endif
