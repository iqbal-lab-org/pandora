#ifndef __PANNODE_H_INCLUDED__   // if pannode.h hasn't been included yet...
#define __PANNODE_H_INCLUDED__

#include <string>
#include <cstdint>
#include <unordered_set>
#include <vector>
#include "kmergraph.h"
#include "localPRG.h"
#include "pangenome/ns.cpp"
#include "vcf.h"


class LocalPRG;

class pangenome::Node {
public:
    std::unordered_multiset<ReadPtr> reads;
    std::vector<SamplePtr> samples;
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

    void add_path(const std::vector<KmerNodePtr> &, const uint32_t &sample_id);

    void get_read_overlap_coordinates(std::vector<std::vector<uint32_t>> &);

    void construct_multisample_vcf(VCF &master_vcf, const std::vector<LocalNodePtr> &,
                                   const std::shared_ptr<LocalPRG> &,
                                   const uint32_t);

    bool operator==(const Node &y) const;

    bool operator!=(const Node &y) const;

    bool operator<(const Node &y) const;

    friend std::ostream &operator<<(std::ostream &out, const Node &n);

    friend class pangenome::Graph;

    friend class pangenome::Read;
};

#endif
