#ifndef __PANNODE_H_INCLUDED__   // if pannode.h hasn't been included yet...
#define __PANNODE_H_INCLUDED__

#include <string>
#include <unordered_set>
#include <vector>
#include "kmergraph.h"

class PanEdge;
class PanRead;
class PanSample;
class KmerNode;
class LocalPRG;

class PanNode {
  public:
    const uint32_t prg_id; // corresponding the the LocalPRG id
    const uint32_t node_id; // unique node id, so can have multiple copies of a localPRG in graph
    const std::string name;
    mutable uint32_t covg;
    KmerGraph kmer_prg;

    std::vector<PanEdge*> edges;
    std::unordered_multiset<PanRead*> reads;
    std::unordered_set<PanSample*> samples;

    PanNode(const uint32_t, const uint32_t, const std::string);
    //PanNode(const PanNode&);
    //PanNode& operator=(const PanNode&);
    
    std::string get_name();
    void add_path(const std::vector<KmerNode*>&);
    void output_samples(const LocalPRG*, const std::string&, const uint);

    bool operator == (const PanNode& y) const;
    bool operator != (const PanNode& y) const;
    bool operator < (const PanNode& y) const;
    friend std::ostream& operator<< (std::ostream& out, const PanNode& n);
};

#endif
