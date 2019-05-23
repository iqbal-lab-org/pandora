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
#include "localPRG.h"
#include "pangenome/ns.cpp"


using KmerNodePtr = std::shared_ptr<KmerNode>;
using ReadId = uint32_t;
using NodeId = uint32_t;


class pangenome::Graph {
protected:
    std::unordered_map<std::string, SamplePtr> samples;
    uint32_t next_id;
public:
    std::map<ReadId, ReadPtr> reads;
    std::unordered_map<NodeId, NodePtr> nodes;

    Graph();

    ~Graph();

    void clear();

    // graph additions/removals
    void reserve_num_reads(uint32_t &);

    ReadPtr get_read(const uint32_t &);

    NodePtr get_node(const NodeId &,
                     const uint32_t &,
                     const std::string &);

    SamplePtr get_sample(const std::string &, const uint32_t &);

    NodePtr add_coverage(ReadPtr &read_ptr,
                         const NodeId &node_id,
                         const uint32_t &prg_id,
                         const std::string &prg_name);

    void add_node(const uint32_t, const std::string &, uint32_t,
                  std::set<MinimizerHitPtr, pComp> &); // used by pandora map
    void add_node(const uint32_t, const std::string &, const std::string &, const uint32_t &sample_id,
                      const std::shared_ptr<LocalPRG> &, const std::vector<KmerNodePtr> &); // used by pandora compare

    std::unordered_map<uint32_t, NodePtr>::iterator remove_node(NodePtr);

    void remove_read(const uint32_t);

    std::vector<NodePtr>::iterator remove_node_from_read(std::vector<NodePtr>::iterator, ReadPtr);

    void remove_low_covg_nodes(const uint32_t &);

    void split_node_by_reads(std::unordered_set<ReadPtr> &, std::vector<uint_least32_t> &, const std::vector<bool> &,
                             const uint_least32_t);

    void setup_kmergraphs(const std::vector<std::shared_ptr<LocalPRG>> &prgs,
                          const uint64_t &total_number_samples = 1);

    //unordered_set<ReadPtr> find_reads_on_node_path(const std::vector<uint16_t>, const std::vector<bool> );
    void add_hits_to_kmergraphs(const std::vector<std::shared_ptr<LocalPRG>> &, const uint32_t &sample_id = 0);

    void copy_coverages_to_kmergraphs(const Graph &, const uint32_t &);

    std::vector<LocalNodePtr>
    infer_node_vcf_reference_path(const Node &, const std::shared_ptr<LocalPRG> &, const uint32_t &,
                                  const std::unordered_map<std::string, std::string> &);

    std::vector<LocalNodePtr>
    get_node_closest_vcf_reference(const Node &, const uint32_t &, const LocalPRG &);

    // graph comparison
    bool operator==(const Graph &y) const;

    bool operator!=(const Graph &y) const;

    // graph read/write
    void save_matrix(const std::string &, const std::vector<std::string> &);

    void save_mapped_read_strings(const std::string &read_filepath, const std::string &outprefix, const int buff = 0);

    friend std::ostream &operator<<(std::ostream &out, const Graph &m);

};

struct same_prg_id {
    uint32_t q;

    same_prg_id(const pangenome::NodePtr &);

    bool operator()(const std::pair<uint32_t, pangenome::NodePtr> &) const;
};

#endif

