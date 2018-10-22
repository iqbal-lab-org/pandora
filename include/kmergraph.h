#ifndef __KMERGRAPH_H_INCLUDED__   // if kmergraph.h hasn't been included yet...
#define __KMERGRAPH_H_INCLUDED__

class LocalPRG;

#include <cstdint>
#include <vector>
#include <iostream>
#include "prg/path.h"
#include "kmernode.h"
#include "pangenome/ns.cpp"


class KmerGraph {
    uint32_t reserved_size;
    uint32_t k;
    float p;
    float nb_p;
    float nb_r;
    int thresh;
public:
    uint32_t exp_depth_covg;
    uint32_t num_reads;
    uint32_t shortest_path_length;
    std::vector<KmerNodePtr> nodes;
    std::vector<KmerNodePtr> sorted_nodes; // representing ordering of the nodes compatible with dp

    KmerGraph();

    KmerGraph(const KmerGraph &);

    KmerGraph &operator=(const KmerGraph &);

    ~KmerGraph();

    void clear();

    KmerNodePtr add_node(const prg::Path &);

    KmerNodePtr add_node_with_kh(const prg::Path &, const uint64_t &, const uint8_t &num = 0);

    void add_edge(KmerNodePtr, KmerNodePtr);

    void remove_shortcut_edges();

    void sort_topologically();

    void check();

    void set_exp_depth_covg(const uint32_t);

    void set_p(const float);

    void set_nb(const float &, const float &);

    float nb_prob(uint32_t);

    float prob(uint32_t);

    float prob(uint32_t, uint32_t);

    float find_max_path(std::vector<KmerNodePtr> &);

    float find_nb_max_path(std::vector<KmerNodePtr> &);

    std::vector<std::vector<KmerNodePtr>> find_max_paths(uint32_t);

    void save_covg_dist(const std::string &);

    uint32_t min_path_length();

    std::vector<std::vector<KmerNodePtr>> get_random_paths(uint32_t);

    float prob_path(const std::vector<KmerNodePtr> &);

    float prob_paths(const std::vector<std::vector<KmerNodePtr>> &);

    void save(const std::string &, const shared_ptr<LocalPRG> = nullptr);

    void load(const std::string &);

    void append_coverages_to_file(const std::string &, const std::string &);

    bool operator==(const KmerGraph &y) const;

    friend std::ostream &operator<<(std::ostream &out, KmerGraph const &data);

    friend void
    estimate_parameters(pangenome::Graph *, const std::string &, const uint32_t, float &, const uint32_t, const bool);

    friend struct condition;

    friend class KmerGraphTest_set_p_Test;

    friend class KmerGraphTest_prob_Test;

    friend class KmerGraphTest_findMaxPathSimple_Test;

    friend class KmerGraphTest_findMaxPath2Level_Test;

    friend class KmerGraphTest_find_max_paths_2Level_Test;

    friend class KmerGraphTest_path_prob_Test;

    friend class KmerGraphTest_path_probs_Test;
};

struct condition {
    prg::Path q;

    condition(const prg::Path &);

    bool operator()(const KmerNodePtr) const;
};

struct pCompKmerNode {
    bool operator()(KmerNodePtr, KmerNodePtr);
};

#endif
