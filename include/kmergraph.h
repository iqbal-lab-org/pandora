#ifndef __KMERGRAPH_H_INCLUDED__   // if kmergraph.h hasn't been included yet...
#define __KMERGRAPH_H_INCLUDED__

class KmerNode;
class PanGraph;
class LocalPRG;

#include <cstring>
#include <map>
#include <vector>
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
    std::unordered_map<uint32_t, KmerNode*> nodes;
    std::vector<KmerNode*> sorted_nodes; // representing ordering of the nodes compatible with dp

    KmerGraph();
    KmerGraph(const KmerGraph&);
    KmerGraph& operator=(const KmerGraph&);
    ~KmerGraph();
    void clear();

    KmerNode* add_node (const Path&);
    KmerNode* add_node_with_kh (const Path&, const uint64_t&, const uint8_t& num=0);
    void add_edge (const Path&, const Path&);
    void add_edge (KmerNode*, KmerNode*);

    void sort_topologically();
    void check();

    void set_p(const float);
    float prob(uint);
    float prob(uint, uint);
    float find_max_path(std::vector<KmerNode*>&);
    float find_min_path(std::vector<KmerNode*>&);
    std::vector<std::vector<KmerNode*>> find_max_paths(uint);
    void save_covg_dist(const std::string&);
    uint min_path_length();
    std::vector<std::vector<KmerNode*>> get_random_paths(uint);
    float prob_path(const std::vector<KmerNode*>&);
    float prob_paths(const std::vector<std::vector<KmerNode*>>&);
    void save (const std::string&);
    void load (const std::string&);
    bool operator == (const KmerGraph& y) const;
    friend std::ostream& operator<< (std::ostream & out, KmerGraph const& data);
    friend void estimate_parameters(PanGraph*, std::string&, uint32_t, float&);
    friend struct condition;
    friend class KmerGraphTest_findMaxPathSimple_Test;
    friend class KmerGraphTest_findMaxPath2Level_Test;
    friend class KmerGraphTest_find_max_paths_2Level_Test;
    friend class KmerGraphTest_path_prob_Test;
};

struct condition
{
    Path q;
    condition(const Path&);
    bool operator()(const std::pair<uint32_t,KmerNode*>&) const;
};

struct pCompKmerNode
{
  bool operator()(KmerNode*, KmerNode*);
};

#endif
