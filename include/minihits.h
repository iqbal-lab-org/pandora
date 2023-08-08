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
    uint32_t number_of_equal_read_minimizers;
    std::vector<uint32_t> *prg_max_path_lengths;  // required to calculate target coverage
public:
    MinimizerHits(std::vector<uint32_t> *prg_max_path_lengths) :
        hits(), number_of_equal_read_minimizers(0), prg_max_path_lengths(prg_max_path_lengths) {};
    ~MinimizerHits() = default;

    // copy constructors
    MinimizerHits(const MinimizerHits& other) = default;
    MinimizerHits& operator=(const MinimizerHits& other) = default;

    // move constructors
    MinimizerHits(MinimizerHits&& other) = default;
    MinimizerHits& operator=(MinimizerHits&& other) = default;

    inline void insert(const MinimizerHitPtr &minimizer_hit) {
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

    inline auto rbegin () const {
        return hits.rbegin();
    }

    inline const MinimizerHitPtr& front() const {
        return *hits.begin();
    }

    inline auto end () const {
        return hits.end();
    }

    inline auto rend () const {
        return hits.rend();
    }

    inline const MinimizerHitPtr& back() const {
        return *hits.rbegin();
    }

    uint32_t read_span_size() const;

    inline void clear() {
        hits.clear();
        number_of_equal_read_minimizers=0;
    }

    inline void increment_number_of_equal_read_minimizers() {
        number_of_equal_read_minimizers++;
    }

    inline uint32_t get_number_of_equal_read_minimizers() const {
        return number_of_equal_read_minimizers;
    }

    inline size_t get_number_of_unique_mini_in_cluster() const {
        return size() - number_of_equal_read_minimizers;
    }

    bool operator<(const MinimizerHits &rhs) const;

    std::pair<uint32_t, uint32_t> get_strand_counts() const;

    /**
     * Returns a proportion denoting the amount of overlap between this cluster and the
     * cluster passed as an argument. The proportion is over the smallest cluster.
     */
    double overlap_amount(const MinimizerHits& cluster) const;

    double target_coverage() const;

    bool is_preferred_to(const MinimizerHits& cluster,
        double minimisers_tolerance) const;

};
#endif
