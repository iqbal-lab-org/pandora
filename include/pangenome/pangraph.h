#ifndef __PANGRAPH_H_INCLUDED__ // if pangraph.h hasn't been included yet...
#define __PANGRAPH_H_INCLUDED__

#include "forward_declarations.h"
#include <string>
#include <cstdint>
#include <unordered_map>
#include <ostream>
#include <vector>
#include <boost/filesystem.hpp>

#include "minihits.h"
#include "localPRG.h"

namespace fs = boost::filesystem;

using KmerNodePtr = std::shared_ptr<KmerNode>;
using NodeId = uint32_t;

class pangenome::Graph {
protected:
    std::unordered_map<std::string, SamplePtr>
        samples; // the samples this pangraph has information
    uint32_t next_id;

public:
    // TODO: move all attributes to private
    std::unordered_map<NodeId, NodePtr> nodes;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // declares all default constructors, destructors and assignment operators
    // explicitly constructor - builds a pangraph to hold information about the given
    // samples
    Graph(const std::vector<std::string>& sample_names = { "sample_1" });
    virtual ~Graph() = default; // destructor
    Graph(const Graph& other) = default; // copy default constructor
    Graph(Graph&& other) = default; // move default constructor
    Graph& operator=(const Graph& other) = default; // copy assignment operator
    Graph& operator=(Graph&& other) = default; // move assignment operator
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // getters
    const NodePtr& get_node(const NodeId& node_id) const { return nodes.at(node_id); }
    const NodePtr& get_node(const std::shared_ptr<LocalPRG>& prg) const
    {
        return nodes.at(prg->id);
    }
    const SamplePtr& get_sample(const std::string& sample_name) const
    {
        return samples.at(sample_name);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // adders
    /**
     * Adds a new node to represent the localPRG given by the parameter.
     * Used by pandora compare.
     * @param prg
     */
    void add_node(const std::shared_ptr<LocalPRG>& prg) { add_node(prg, prg->id); }
    /**
     * Adds a new node to represent the localPRG given by the parameter.
     * Used by pandora compare.
     * @param prg
     * @param node_id : sometimes we need more than 1 node representing a single PRG.
     * node_id provides uniqueness to each node, so that is why we have an over
     */
    void add_node(const std::shared_ptr<LocalPRG>& prg, uint32_t node_id);

    /**
     * Record we have hit the given PRG
     */
    void record_hit(const std::shared_ptr<LocalPRG>& prg);

    /**
     * Adds hits between the given PRG and sample described as a path of minimizer kmers
     * from the consensus path. This is just used in the global pangraph in pandora
     * compare. For each sample, we go to each PRG/node and get the consensus path (i.e.
     * the path in the PRG with the best read support from that sample). From this
     * consensus path, we get the minimizer kmers, and we add this minimizer kmer path
     * here. This is called when the pangraph represents a COLLECTION of samples - the
     * coverage of the PRG/node is the samples mapping to this PRG/node
     * @param node
     * @param sample
     * @param kmp : the minimizer kmers node path
     */
    void add_hits_between_PRG_and_sample(const NodePtr& node, const SamplePtr& sample,
        const std::vector<KmerNodePtr>& kmp);
    void add_hits_between_PRG_and_sample(uint32_t node_id,
        const std::string& sample_name, const std::vector<KmerNodePtr>& kmp)
    {
        add_hits_between_PRG_and_sample(
            get_node(node_id), get_sample(sample_name), kmp);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Remove nodes with covg <= thresh from graph
     * @param thresh
     */
    void remove_low_covg_nodes(const uint32_t& thresh);

    // TODO: possibly refactor the methods below
    std::unordered_map<uint32_t, NodePtr>::iterator remove_node(NodePtr);

    void copy_coverages_to_kmergraphs(const Graph&, const uint32_t&);
    std::vector<LocalNodePtr> infer_node_vcf_reference_path(const Node&,
        const std::shared_ptr<LocalPRG>&, const uint32_t&,
        const std::unordered_map<std::string, std::string>&, const uint32_t&) const;
    std::vector<LocalNodePtr> get_node_closest_vcf_reference(const Node&,
        const uint32_t&, const LocalPRG&,
        const uint32_t& max_num_kmers_to_average) const;
    // graph comparison
    bool operator==(const Graph& y) const;
    bool operator!=(const Graph& y) const;
    // graph read/write
    void save_matrix(
        const fs::path& filepath, const std::vector<std::string>& sample_names);
    friend std::ostream& operator<<(std::ostream& out, const Graph& m);
};

struct same_prg_id {
    uint32_t q;

    same_prg_id(const pangenome::NodePtr&);

    bool operator()(const std::pair<uint32_t, pangenome::NodePtr>&) const;
};

#endif
