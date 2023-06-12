#include <iostream>
#include <unordered_map>
#include <set>
#include <memory>
#include <vector>
#include <algorithm>
#include <boost/filesystem.hpp>

#include "utils.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "pangenome/pansample.h"
#include "fatal_error.h"

using namespace pangenome;

pangenome::Graph::Graph(const std::vector<std::string>& sample_names)
    : next_id { 0 }
{
    nodes.reserve(6000);

    // add the samples
    uint32_t sample_id = 0;
    for (const auto& sample_name : sample_names) {
        auto samplePointer = std::make_shared<Sample>(sample_name, sample_id++);
        this->samples[sample_name] = samplePointer;
    }
}

void pangenome::Graph::add_node(const std::shared_ptr<LocalPRG>& prg, uint32_t node_id)
{
    NodePtr node_ptr;
    auto it = nodes.find(node_id);
    const bool found_node = it != nodes.end();
    if (not found_node) {
        node_ptr = std::make_shared<Node>(
            prg, node_id, samples.size()); // TODO: refactor this - holding the
                                           // reference to PRG is enough
        nodes[node_id] = node_ptr;
    }
}

void pangenome::Graph::record_hit(const std::shared_ptr<LocalPRG>&prg)
{
    // add and get the new node
    add_node(prg);
    auto node_ptr = get_node(prg);

    // update the coverage
    node_ptr->covg += 1;
    node_ptr->kmer_prg_with_coverage.set_num_reads(node_ptr->covg);
}

// TODO: this should be a method of class Sample
void update_sample_info_with_this_node(
    const SamplePtr& sample, const NodePtr& node, const std::vector<KmerNodePtr>& kmp)
{
    sample->add_path(node->node_id, kmp);
}

// TODO: this should be a method of class Node
void update_node_info_with_this_sample(const NodePtr& node, const SamplePtr& sample)
{
    if (node->samples.find(sample) == node->samples.end())
        node->samples.insert(sample);
    node->covg += 1;
}

void pangenome::Graph::add_hits_between_PRG_and_sample(
    const NodePtr& node, const SamplePtr& sample, const std::vector<KmerNodePtr>& kmp)
{
    update_sample_info_with_this_node(sample, node, kmp);
    update_node_info_with_this_sample(node, sample);
}

// Remove the node n, and all references to it
std::unordered_map<uint32_t, NodePtr>::iterator pangenome::Graph::remove_node(NodePtr n)
{
    auto it = nodes.find(n->node_id);
    if (it != nodes.end()) {
        it = nodes.erase(it);
    }
    return it;
}

// remove the all instances of the pattern of nodes/orienations from graph
void pangenome::Graph::remove_low_covg_nodes(const uint32_t& thresh)
{
    BOOST_LOG_TRIVIAL(debug) << "Remove nodes with covg <= " << thresh << std::endl;
    for (auto it = nodes.begin(); it != nodes.end();) {
        if (it->second->covg <= thresh) {
            it = remove_node(it->second);
        } else {
            ++it;
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "Pangraph now has " << nodes.size() << " nodes";
}


// For each node in reference pangraph, copy the coverages over to sample_id in this
// pangraph
void pangenome::Graph::copy_coverages_to_kmergraphs(
    const Graph& ref_pangraph, const uint32_t& sample_id)
{
    const uint32_t ref_sample_id = 0;
    for (const auto& ref_node_entry : ref_pangraph.nodes) {
        const Node& ref_node = *ref_node_entry.second;

        const bool ref_node_is_in_this_pangraph
            = nodes.find(ref_node.node_id) != nodes.end();
        if (!ref_node_is_in_this_pangraph) {
            fatal_error(
                "Error copying coverages to kmer graphs: reference node does not "
                "exist in pangraph");
        }

        Node& pangraph_node = *nodes[ref_node.node_id];
        for (auto& kmergraph_node_ptr :
            pangraph_node.kmer_prg_with_coverage.kmer_prg->nodes) {
            const auto& knode_id = kmergraph_node_ptr->id;

            const bool kmer_graph_node_id_is_valid
                = knode_id < ref_node.kmer_prg_with_coverage.kmer_prg->nodes.size();
            if (!kmer_graph_node_id_is_valid) {
                fatal_error("Error copying coverages to kmer graphs: kmer graph node "
                            "id is not valid");
            }

            pangraph_node.kmer_prg_with_coverage.set_reverse_covg(knode_id,
                (uint16_t)(ref_node.kmer_prg_with_coverage.get_reverse_covg(
                    knode_id, ref_sample_id)),
                sample_id);
            pangraph_node.kmer_prg_with_coverage.set_forward_covg(knode_id,
                (uint16_t)(ref_node.kmer_prg_with_coverage.get_forward_covg(
                    knode_id, ref_sample_id)),
                sample_id);
        }
    }
}

std::vector<LocalNodePtr> pangenome::Graph::infer_node_vcf_reference_path(
    const Node& node, const std::shared_ptr<LocalPRG>& prg_ptr, const uint32_t& w,
    const std::unordered_map<std::string, std::string>& vcf_refs,
    const uint32_t& max_num_kmers_to_average) const
{
    BOOST_LOG_TRIVIAL(info) << "Infer VCF reference path";
    const auto& prg = *prg_ptr;
    if (vcf_refs.find(prg.name) != vcf_refs.end()) {
        const auto& vcf_reference_sequence = vcf_refs.at(prg.name);
        auto reference_path = prg.get_valid_vcf_reference(vcf_reference_sequence);
        if (!reference_path.empty())
            return reference_path;
    }
    return get_node_closest_vcf_reference(node, w, prg, max_num_kmers_to_average);
}

std::vector<LocalNodePtr> pangenome::Graph::get_node_closest_vcf_reference(
    const Node& node, const uint32_t& w, const LocalPRG& prg,
    const uint32_t& max_num_kmers_to_average) const
{
    // TODO: check if this is correct
    auto kmer_prg_with_coverage
        = node.kmer_prg_with_coverage; // TODO: is this indeed an assignment op?
    kmer_prg_with_coverage.zeroCoverages();

    for (const auto& sample_entry : this->samples) {
        const auto& sample = sample_entry.second;
        if (sample->paths.find(node.prg_id) == sample->paths.end()) {
            BOOST_LOG_TRIVIAL(debug) << "Could not find sample path for sample "
                                     << sample->name << " and prg " << node.prg_id;
            continue;
        }

        BOOST_LOG_TRIVIAL(debug) << "Increment coverages for path for sample "
                                 << sample->name << " and prg " << node.prg_id;
        const auto& sample_paths = sample->paths.at(node.prg_id);
        for (const auto& sample_path : sample_paths) {
            for (uint32_t i = 0; i != sample_path.size(); ++i) {
                const bool sample_path_node_is_valid
                    = (sample_path[i]->id
                          < kmer_prg_with_coverage.kmer_prg->nodes.size())
                    and (kmer_prg_with_coverage.kmer_prg->nodes[sample_path[i]->id]
                        != nullptr);
                if (!sample_path_node_is_valid) {
                    fatal_error("When getting the path closest to VCF reference, "
                                "a sample path node is not valid");
                }

                kmer_prg_with_coverage.increment_forward_covg(sample_path[i]->id, 0);
                kmer_prg_with_coverage.increment_reverse_covg(sample_path[i]->id, 0);
            }
        }
    }

    kmer_prg_with_coverage.kmer_prg->discover_k();
    kmer_prg_with_coverage.set_num_reads(node.covg);

    std::vector<KmerNodePtr> kmer_path;
    kmer_prg_with_coverage.find_max_path(kmer_path, "lin", max_num_kmers_to_average, 0);
    if (!kmer_path.empty()) {
        auto reference_path = prg.localnode_path_from_kmernode_path(kmer_path, w);
        BOOST_LOG_TRIVIAL(debug) << "Found reference path to return";
        return reference_path;
    } else {
        BOOST_LOG_TRIVIAL(debug) << "Could not infer reference path, so use top path";
        return prg.prg.top_path();
    }
}

same_prg_id::same_prg_id(const NodePtr& p)
    : q(p->prg_id) {};

bool same_prg_id::operator()(const std::pair<uint32_t, NodePtr>& n) const
{
    return (n.second->prg_id == q);
}

bool pangenome::Graph::operator==(const Graph& y) const
{
    // false if have different nodes
    for (const auto& c : nodes) {
        // if node id doesn't exist
        auto it = find_if(y.nodes.begin(), y.nodes.end(), same_prg_id(c.second));
        if (it == y.nodes.end()) {
            BOOST_LOG_TRIVIAL(warning) << "can't find node " << c.first;
            return false;
        }
    }
    for (const auto& c : y.nodes) {
        // if node id doesn't exist
        auto it = find_if(nodes.begin(), nodes.end(), same_prg_id(c.second));
        if (it == nodes.end()) {
            BOOST_LOG_TRIVIAL(debug) << "can't find node " << c.first;
            return false;
        }
    }

    // otherwise is true
    return true;
}

bool pangenome::Graph::operator!=(const Graph& y) const { return !(*this == y); }

// Saves a presence/absence/copynumber matrix for each node and each sample
void pangenome::Graph::save_matrix(
    const fs::path& filepath, const std::vector<std::string>& sample_names)
{
    // write a presence/absence matrix for samples and nodes
    fs::ofstream handle;
    handle.open(filepath);

    // save header line with sample names
    for (const auto& name : sample_names) {
        handle << "\t" << name;
    }
    handle << std::endl;

    // for each node, save number of each sample
    for (const auto& n : nodes) {
        handle << n.second->name;
        for (const auto& name : sample_names) {
            if (samples.find(name) == samples.end()
                or samples[name]->paths.find(n.second->node_id)
                    == samples[name]->paths.end()) {
                handle << "\t0";
            } else {
                handle << "\t" << samples[name]->paths[n.second->node_id].size();
            }
        }
        handle << std::endl;
    }
}

std::ostream& pangenome::operator<<(std::ostream& out, const pangenome::Graph& m)
{
    for (const auto& n : m.nodes) {
        out << n.second->prg_id << std::endl;
        out << *n.second << std::endl;
    }
    return out;
}
