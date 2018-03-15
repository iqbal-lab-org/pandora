#ifndef __PANNODE_H_INCLUDED__   // if pannode.h hasn't been included yet...
#define __PANNODE_H_INCLUDED__

#include <string>
#include <unordered_set>
#include <vector>
#include "kmergraph.h"
#include "pangenome/ns.cpp"


class KmerNode;

class LocalPRG;

class pangenome::Node {
public:
    std::unordered_multiset<ReadPtr> reads;
    std::unordered_set<SamplePtr> samples;
    const uint32_t prg_id; // corresponding the the LocalPRG id
    const uint32_t node_id; // unique node id, so can have multiple copies of a localPRG in graph
    const std::string name;
    mutable uint32_t covg;
    KmerGraph kmer_prg;

    Node(const uint32_t, const uint32_t, const std::string);
    //Node(const Node&);
    //Node& operator=(const Node&);

    void remove_read(ReadPtr);

    std::string get_name() const;

    void add_path(const std::vector<KmerNodePtr> &);

    void output_samples(const LocalPRG *, const std::string &, const uint, const std::string &);

    bool operator==(const Node &y) const;

    bool operator!=(const Node &y) const;

    bool operator<(const Node &y) const;

    friend std::ostream &operator<<(std::ostream &out, const Node &n);

    friend class pangenome::Graph;

    friend class pangenome::Read;
};

#endif
