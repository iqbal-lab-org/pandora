#ifndef __KMERGRAPH_H_INCLUDED__   // if kmergraph.h hasn't been included yet...
#define __KMERGRAPH_H_INCLUDED__

class KmerNode;
class PanGraph;
class LocalPRG;

#include <cstring>
#include <map>
#include <vector>
#include <deque>
#include <unordered_map>
#include <iostream>
#include "path.h"
#include "kmernode.h"

class KmerGraph {
    uint reserved_size;
    uint32_t next_id;
    uint32_t k;
    float p;
    int thresh;
  public:
    uint32_t num_reads;
    uint32_t shortest_path_length;
    std::unordered_map<uint32_t, KmerNodePtr> nodes;
    std::vector<KmerNodePtr> sorted_nodes; // representing ordering of the nodes compatible with dp
    std::vector<std::vector<std::vector<bool>>> covgs;


    KmerGraph();
    KmerGraph(const KmerGraph&);
    KmerGraph& operator=(const KmerGraph&);
    ~KmerGraph();
    void clear();

    KmerNodePtr add_node (const Path&);
    KmerNodePtr add_node_with_kh (const Path&, const uint64_t&, const uint8_t& num=0);
    void add_edge (const Path&, const Path&);
    void add_edge (KmerNodePtr, KmerNodePtr);

    void sort_topologically();
    void check();

    void get_prev(const uint16_t, const uint8_t, const uint16_t, uint16_t&, std::vector<std::deque<KmerNodePtr>>&);
    void get_prev(const uint16_t, const uint8_t, uint16_t&, std::vector<std::deque<KmerNodePtr>>&);
    void get_next(const uint16_t, const uint8_t, const uint16_t, uint16_t&, std::vector<std::deque<KmerNodePtr>>&);
    void get_next(const uint16_t, const uint8_t, uint16_t&, std::vector<std::deque<KmerNodePtr>>&);
    void extend_paths_back(std::vector<std::deque<KmerNodePtr>>&, const std::vector<std::deque<KmerNodePtr>>&);
    void extend_paths_forward(std::vector<std::deque<KmerNodePtr>>&, const std::vector<std::deque<KmerNodePtr>>&);
    //void find_compatible_paths(const uint16_t, std::vector<std::deque<KmerNodePtr>>&);
    void find_compatible_paths(const uint8_t, std::vector<std::deque<KmerNodePtr>>&);
    void find_all_compatible_paths(std::vector<std::deque<KmerNodePtr>>&, std::vector<std::vector<std::pair<uint16_t, uint16_t>>>&, const uint8_t thresh=16);

    void set_p(const float);
    float prob(uint);
    float prob(uint, uint);

    float find_max_path(std::vector<KmerNodePtr>&);
    float find_min_path(std::vector<KmerNodePtr>&);
    std::vector<std::vector<KmerNodePtr>> find_max_paths(uint);
    void save_covg_dist(const std::string&);
    uint min_path_length();
    std::vector<std::vector<KmerNodePtr>> get_random_paths(uint);
    float prob_path(const std::vector<KmerNodePtr>&);
    float prob_paths(const std::vector<std::vector<KmerNodePtr>>&);
    void save (const std::string&);
    void load (const std::string&);
    bool operator == (const KmerGraph& y) const;
    friend std::ostream& operator<< (std::ostream & out, KmerGraph const& data);
    friend void estimate_parameters(PanGraph*, const std::string&, const uint32_t, float&, const uint);
    friend struct condition;
    friend class KmerGraphTest_set_p_Test;
    friend class KmerGraphTest_prob_Test;
    friend class KmerGraphTest_findMaxPathSimple_Test;
    friend class KmerGraphTest_findMaxPath2Level_Test;
    friend class KmerGraphTest_find_max_paths_2Level_Test;
    friend class KmerGraphTest_path_prob_Test;
    friend class KmerGraphTest_path_probs_Test;
};

struct condition
{
    Path q;
    condition(const Path&);
    bool operator()(const std::pair<uint32_t,KmerNodePtr>&) const;
};

struct pCompKmerNode
{
  bool operator()(KmerNodePtr, KmerNodePtr);
};

#endif
