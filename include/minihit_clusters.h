#ifndef PANDORA_MINIHIT_CLUSTERS_H
#define PANDORA_MINIHIT_CLUSTERS_H

#include "minihits.h"
#include <algorithm>
#include <random>


/**
 * Class that represents a set of minimizer hit clusters to be further filtered out.
 * Firstly, clusters have to be inserted into this set, and then insertions must be
 * explicitly finalised to allow for queries.
 * When insertions are finalised, we put the clusters inserted in random order and
 * stably sort them, which will make the the clusters sorted, but in places where we
 * have the same mappings, the order will be random. The cluster filtering algorithm
 * favours the first cluster, this will cause multimapping clusters to be randomly selected,
 * avoiding deterministic multimapping.
 * Querying can only happen after insertions are finalised.
 */
class MinimizerHitClusters {
public:
    explicit MinimizerHitClusters(const uint32_t rng_seed) : insertion_phase(true), clusters(){
        if (bool deterministic_run = rng_seed > 0; deterministic_run) {
            rng = std::mt19937(rng_seed);
        } else{
            rng = std::mt19937(std::random_device()());
        }
    }
    ~MinimizerHitClusters() = default;  // Note: this is not virtual as there is no need to

    // enable copy/move
    MinimizerHitClusters(const MinimizerHitClusters &) = default;  // copy constructor
    MinimizerHitClusters& operator=(const MinimizerHitClusters &) = default; // copy assignment operator
    MinimizerHitClusters(MinimizerHitClusters &&) = default; // move constructor
    MinimizerHitClusters& operator=(MinimizerHitClusters &&) = default; // move assignment operator

    inline void insert(const MinimizerHits& cluster) {
        check_if_can_insert();
        clusters.emplace_back(cluster);
    }

    inline void finalise_insertions() {
        check_if_can_insert();
        insertion_phase=false;

        std::shuffle(clusters.begin(), clusters.end(), rng);
        std::stable_sort(clusters.begin(), clusters.end());

        check_if_can_query();
    }

    inline auto size() const {
        check_if_can_query();
        return clusters.size();
    }

    inline auto empty() const {
        check_if_can_query();
        return clusters.empty();
    }

    inline auto begin () const {
        check_if_can_query();
        return clusters.begin();
    }
    inline auto begin () {
        check_if_can_query();
        return clusters.begin();
    }

    inline auto end () const {
        check_if_can_query();
        return clusters.end();
    }
    inline auto end () {
        check_if_can_query();
        return clusters.end();
    }

private:
    bool insertion_phase;
    std::vector<MinimizerHits> clusters;
    std::mt19937 rng;

    inline void check_if_can_insert() const {
        if (!insertion_phase) {
            fatal_error("Tried to insert a cluster in a MinimizerHitClusters when insertion phase was finished");
        }
    }

    inline void check_if_can_query() const {
        if (insertion_phase) {
            fatal_error("Tried to query a MinimizerHitClusters when insertion phase was not finished");
        }
    }
};

#endif // PANDORA_MINIHIT_CLUSTERS_H
