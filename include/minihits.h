#ifndef __MINIHITS_H_INCLUDED__   // if minihits.h hasn't been included yet...
#define __MINIHITS_H_INCLUDED__

struct MinimizerHit;

#include <set>
#include <unordered_set>
#include <cstdint>
#include <memory>
#include "minimizer.h"
#include "minirecord.h"

typedef std::shared_ptr<MinimizerHit> MinimizerHitPtr;

struct pComp {
    bool operator()(const MinimizerHitPtr &lhs, const MinimizerHitPtr &rhs);
};

struct pEq {
    bool operator()(const MinimizerHitPtr &lhs, const MinimizerHitPtr &rhs) const;
};

struct Hash {
    size_t operator()(const MinimizerHit *mh) const;
};

struct pComp_path {
    bool operator()(const MinimizerHitPtr &lhs, const MinimizerHitPtr &rhs);
};

struct clusterComp {
    bool operator()(std::set<MinimizerHitPtr, pComp> lhs, std::set<MinimizerHitPtr, pComp> rhs);
};

struct clusterComp_size {
    bool operator()(std::set<MinimizerHitPtr, pComp> lhs, std::set<MinimizerHitPtr, pComp> rhs);
};

class MinimizerHits {
public:
    MinimizerHits(const uint &num_hits = 30000);

    ~MinimizerHits();

    void clear();

    //std::unordered_set<MinimizerHit*, Hash, pEq> uhits;
    std::unordered_set<MinimizerHitPtr> uhits;
    std::set<MinimizerHitPtr, pComp> hits;

    void add_hit(const uint32_t i, const Minimizer& m, const MiniRecord *r);

    void sort();

    friend std::ostream &operator<<(std::ostream &out, const MinimizerHits &m);
};

#endif
