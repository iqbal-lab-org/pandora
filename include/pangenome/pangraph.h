#ifndef __PANGRAPH_H_INCLUDED__   // if pangraph.h hasn't been included yet...
#define __PANGRAPH_H_INCLUDED__

struct MinimizerHit;
class LocalPRG;
class KmerNode;

#include <cstring>
#include <unordered_map>
#include <ostream>
#include <functional>
#include <memory>
#include "minihits.h"
#include "pangenome/ns.cpp"


typedef std::shared_ptr<KmerNode> KmerNodePtr;

class pangenome::Graph {
protected:
    std::unordered_map<std::string, SamplePtr> samples;
    uint32_t next_id;
public:
    std::unordered_map<uint32_t, ReadPtr> reads;
    std::unordered_map<uint32_t, NodePtr> nodes;

    Graph();
    ~Graph();
    void clear();

    // graph additions/removals
    void add_node(const uint32_t, const std::string, uint32_t,
                  const std::set<MinimizerHitPtr, pComp> &); // used by pandora map
    void add_node(const uint32_t, const std::string &, const std::string &, const std::vector<KmerNodePtr> &,
                  const LocalPRG *); // used by pandora compare

    std::unordered_map<uint32_t, NodePtr>::iterator remove_node(NodePtr);
    void remove_read(const uint32_t);

    void remove_low_covg_nodes(const uint &);

    void split_node_by_reads(const unordered_set<ReadPtr>&, vector<uint16_t>&, const vector<bool>&, const uint16_t);
    //unordered_set<ReadPtr> find_reads_on_node_path(const std::vector<uint16_t>, const std::vector<bool> );
    void add_hits_to_kmergraphs(const std::vector<LocalPRG *> &);

    // graph comparison
    bool operator==(const Graph &y) const;
    bool operator!=(const Graph &y) const;

    // graph read/write
    void save_matrix(const std::string &);

    friend std::ostream &operator<<(std::ostream &out, const Graph &m);

};

struct same_prg_id
{
    uint32_t q;
    same_prg_id(const pangenome::NodePtr&);
    bool operator()(const std::pair<uint32_t,pangenome::NodePtr>&) const;
};
#endif

