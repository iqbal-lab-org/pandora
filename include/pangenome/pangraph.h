#ifndef __PANGRAPH_H_INCLUDED__   // if pangraph.h hasn't been included yet...
#define __PANGRAPH_H_INCLUDED__

struct MinimizerHit;

class LocalPRG;

class KmerNode;

#include <string>
#include <cstdint>
#include <unordered_map>
#include <ostream>
#include <vector>
#include "minihits.h"
#include "pangenome/ns.cpp"


using KmerNodePtr = std::shared_ptr<KmerNode>;
using ReadId = uint32_t;
using NodeId = uint32_t;


class pangenome::Graph {
protected:
    std::unordered_map<std::string, SamplePtr> samples;
    uint32_t next_id;
public:
    std::unordered_map<ReadId, ReadPtr> reads;
    std::unordered_map<NodeId, NodePtr> nodes;

    Graph();

    ~Graph();

    void clear();

    // graph additions/removals
    void reserve_num_reads(uint32_t &);

    ReadPtr get_read(const uint32_t &read_id);

    NodePtr add_coverage(ReadPtr &read_ptr,
                         const NodeId &node_id,
                         const uint32_t &prg_id,
                         const std::string &prg_name);

    void add_node(const uint32_t, const std::string &, uint32_t,
                  std::set<MinimizerHitPtr, pComp> &); // used by pandora map
    void add_node(const uint32_t, const std::string &, const std::string &, const std::vector<KmerNodePtr> &,
                  const std::shared_ptr<LocalPRG> &); // used by pandora compare

    std::unordered_map<uint32_t, NodePtr>::iterator remove_node(NodePtr);

    void remove_read(const uint32_t);

    std::vector<NodePtr>::iterator remove_node_from_read(std::vector<NodePtr>::iterator, ReadPtr);

    void remove_low_covg_nodes(const uint32_t &);

    void split_node_by_reads(std::unordered_set<ReadPtr> &, std::vector<uint_least32_t> &, const std::vector<bool> &,
                             const uint_least32_t);

    //unordered_set<ReadPtr> find_reads_on_node_path(const std::vector<uint16_t>, const std::vector<bool> );
    void add_hits_to_kmergraphs(const std::vector<std::shared_ptr<LocalPRG>> &);

    // graph comparison
    bool operator==(const Graph &y) const;

    bool operator!=(const Graph &y) const;

    // graph read/write
    void save_matrix(const std::string &);

    void save_mapped_read_strings(const std::string &read_filepath, const std::string &outprefix, const int buff = 0);

    void save_kmergraph_coverages(const std::string &, const std::string &);

    friend std::ostream &operator<<(std::ostream &out, const Graph &m);

};

struct same_prg_id {
    uint32_t q;

    same_prg_id(const pangenome::NodePtr &);

    bool operator()(const std::pair<uint32_t, pangenome::NodePtr> &) const;
};

#endif

