#ifndef __LOCALNODE_H_INCLUDED__ // if localnode.h hasn't been included yet...
#define __LOCALNODE_H_INCLUDED__

#include "interval.h"
#include "kmernode.h"
#include "prg/path.h"
#include <cstdint>
#include <cstring>
#include <memory>
#include <ostream>
#include <unordered_set>
#include <vector>

class LocalNode;

typedef std::shared_ptr<LocalNode> LocalNodePtr;

class LocalNode {
    std::unordered_set<KmerNodePtr> prev_kmer_paths;

public:
    std::string seq;
    Interval pos; // pos in the prg
    uint32_t id;
    uint32_t covg; // covg by hits - initially has the size of the interval
    uint32_t sketch_next; // used by minimizer_sketch function in localPRG.cpp - the
                          // next position to sketch?
    bool skip; // used by minimizer_sketch function in localPRG.cpp

    std::vector<LocalNodePtr>
        outNodes; // representing edges from this node to the nodes in the vector

    LocalNode(std::string, Interval, uint32_t);

    bool operator==(const LocalNode& y) const;

    friend std::ostream& operator<<(std::ostream& out, const LocalNode& n);

    friend class LocalGraph;

    friend class LocalPRG;
};

#endif
