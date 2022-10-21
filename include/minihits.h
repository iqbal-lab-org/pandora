#ifndef __MINIHITS_H_INCLUDED__ // if minihits.h hasn't been included yet...
#define __MINIHITS_H_INCLUDED__

#include "forward_declarations.h"
#include <set>
#include <unordered_set>
#include <memory>
#include "minimizer.h"
#include "minirecord.h"
#include <memory>

struct pComp {
    bool operator()(const MinimizerHitPtr& lhs, const MinimizerHitPtr& rhs) const;
};

struct pCompReadPositionFirst {
    bool operator()(const MinimizerHitPtr& lhs, const MinimizerHitPtr& rhs) const;
};

struct pEq {
    bool operator()(const MinimizerHitPtr& lhs, const MinimizerHitPtr& rhs) const;
};


class MinimizerHits {
private:
    std::set<MinimizerHitPtr, pComp> hits;

public:
    MinimizerHits() = default;
    ~MinimizerHits() = default;

    inline void insert(const MinimizerHitPtr minimizer_hit) {
        hits.insert(minimizer_hit);
    }

    inline void insert(decltype(hits.begin()) begin, decltype(hits.end()) end) {
        while (begin != end) {
            this->insert(*begin);
            begin++;
        }
    }

    void insert(const uint32_t i, const Minimizer& minimizer_from_read,
        const MiniRecord& minimizer_from_PRG);

    inline auto size() const {
        return hits.size();
    }

    inline auto empty() const {
        return hits.empty();
    }

    inline auto begin () const {
        return hits.begin();
    }

    inline auto end () const {
        return hits.end();
    }

    inline void clear() { hits.clear(); }

    bool operator<(const MinimizerHits &rhs) const;

    std::pair<uint32_t, uint32_t> get_strand_counts() const;
};
#endif
