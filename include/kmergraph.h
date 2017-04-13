#ifndef __KMERGRAPH_H_INCLUDED__   // if kmergraph.h hasn't been included yet...
#define __KMERGRAPH_H_INCLUDED__

class KmerNode;

#include <cstring>
#include <map>
#include <vector>
#include <iostream>
#include "path.h"
#include "kmernode.h"

class KmerGraph {
    uint32_t next_id;
    uint32_t k;
    float p;
  public:
    uint32_t num_reads;
    uint32_t shortest_path_length;
    std::vector<KmerNode*> nodes; // representing nodes in graph
    KmerGraph();
    ~KmerGraph();
    void clear();
    KmerNode* add_node (const Path&);
    KmerNode* add_node_with_kh (const Path&, const uint64_t&);
    //void add_edge (const uint32_t&, const uint32_t&);
    void add_edge (const Path&, const Path&);
    void add_edge (KmerNode*, KmerNode*);
    //void copy_innodes (const Path&, const Path&);
    std::set<Path> get_innodes (const Path&);
    void check (uint);
    void sort_topologically();
    //std::vector<KmerNode*> get_node_order();
    float prob(uint);
    float find_max_path(float, std::vector<KmerNode*>&);
    //float find_max_path_backward(int, float, std::vector<KmerNode*>&);
    //float find_max_path_forward(int, float, std::vector<KmerNode*>&);
    //float find_max_path_coverage(int, float, std::vector<KmerNode*>&);
    uint min_path_length();
    void save (const std::string&);
    void load (const std::string&);
    bool operator == (const KmerGraph& y) const;
    friend std::ostream& operator<< (std::ostream & out, KmerGraph const& data);
    friend struct condition;
    friend class KmerGraphTest_findMaxPathSimple_Test;
    friend class KmerGraphTest_findMaxPath2Level_Test;
};

struct condition
{
    Path q;
    condition(const Path&);
    bool operator()(const KmerNode*) const;
};

struct pCompKmerNode
{
  bool operator()(KmerNode*, KmerNode*);
};

#endif
