#ifndef __KMERNODE_H_INCLUDED__   // if kmernode.h hasn't been included yet...
#define __KMERNODE_H_INCLUDED__

class KmerNode;

#include <cstring>
#include <cstdint>
#include <vector>
#include <ostream>
#include <memory>

#include "prg/path.h"
#include "pangenome/ns.cpp"


typedef std::shared_ptr<KmerNode> KmerNodePtr;
typedef prg::Path Path;

class KmerNode {

public:
    uint32_t id;
    Path path;
    std::vector<KmerNodePtr> outNodes; // representing edges from this node to the nodes in the vector
    std::vector<KmerNodePtr> inNodes; // representing edges from other nodes to this node
    std::vector<std::pair<uint32_t, uint32_t>> covg_new; // sample covg by hits in fwd, rev dir
    uint64_t khash; //the kmer hash value
    uint8_t num_AT; // the number of As and Ts in this kmer

    KmerNode(uint32_t, const Path &);

    KmerNode(const KmerNode &);

    KmerNode &operator=(const KmerNode &);

    bool operator==(const KmerNode &y) const;

    void increment_covg(const bool &, const uint32_t &sample_id = 0);

    uint32_t get_covg(const bool &strand, const uint32_t &sample_id);

    void set_covg(const uint32_t &, const bool &, const uint32_t &);

    friend std::ostream &operator<<(std::ostream &out, const KmerNode &n);

    friend class KmerGraph;

    friend struct condition;
    friend struct pCompKmerNode;

    friend class LocalPRG;

    friend class pangenome::Graph;

    friend class pangenome::Node;

    friend int pandora_check_kmergraph(int argc, char *argv[]);

    /*
    friend void estimate_parameters(pangenome::Graph *, const std::string &, const uint32_t, float &, const uint32_t);
     */

};

#endif
