#ifndef __MINIHITS_H_INCLUDED__ // if minihits.h hasn't been included yet...
#define __MINIHITS_H_INCLUDED__

#include <set>
#include <unordered_set>
#include <memory>
#include "minimizer.h"
#include "minirecord.h"

struct MinimizerHit;
struct pComp;
typedef std::shared_ptr<MinimizerHit> MinimizerHitPtr;
typedef std::set<MinimizerHitPtr, pComp> MinimizerHitCluster;

struct pComp {
    bool operator()(const MinimizerHitPtr& lhs, const MinimizerHitPtr& rhs);
};

struct pEq {
    bool operator()(const MinimizerHitPtr& lhs, const MinimizerHitPtr& rhs) const;
};

struct Hash {
    size_t operator()(const MinimizerHit* mh) const;
};

struct pComp_path {
    bool operator()(const MinimizerHitPtr& lhs, const MinimizerHitPtr& rhs);
};

struct clusterComp {
    bool operator()(const MinimizerHitCluster lhs, const MinimizerHitCluster rhs);
};

struct clusterComp_size {
    bool operator()(const MinimizerHitCluster lhs, const MinimizerHitCluster rhs);
};

class MinimizerHits {
public:
    MinimizerHits() = default;
    ~MinimizerHits() = default;

    std::set<MinimizerHitPtr, pComp> hits;

    void add_hit(const uint32_t i, const Minimizer& minimizer_from_read,
        const MiniRecord& minimizer_from_PRG);

    void clear() { hits.clear(); }

    // friend std::ostream &operator<<(std::ostream &out, const MinimizerHits &m);
};

#endif
