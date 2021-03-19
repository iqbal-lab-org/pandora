#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <memory>
#include <vector>
#include <algorithm>
#include <boost/filesystem.hpp>

#include "utils.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "pangenome/pansample.h"
#include "fastaq_handler.h"
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

void pangenome::Graph::add_read(const uint32_t& read_id)
{
    auto it = reads.find(read_id);
    const bool found = it != reads.end();
    if (not found) {
        auto read_ptr = std::make_shared<Read>(read_id);
        reads[read_id] = read_ptr;
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

/**
 * Updated the node information with this new read mapping to it
 * Increments node coverage by 1 and add the read to the set of reads the read cover
 * @param node_ptr
 * @param read_ptr
 */
// TODO: this should be a method of the Node class
void update_node_info_with_this_read(const NodePtr& node_ptr, const ReadPtr& read_ptr)
{
    node_ptr->covg += 1;
    node_ptr->reads.insert(read_ptr);

    const bool coverage_information_is_consistent_with_read_information
        = node_ptr->covg == node_ptr->reads.size();
    if (!coverage_information_is_consistent_with_read_information) {
        fatal_error("Error updating Pangraph node with read: coverage information "
                    "is not consistent with read information");
    }
}

// Checks that all hits in the cluster are from the given prg and read
void check_correct_hits(const uint32_t prg_id, const uint32_t read_id,
    const std::set<MinimizerHitPtr, pComp>& cluster)
{
    for (const auto& hit_ptr : cluster) {
        const bool hits_correspond_to_correct_read = read_id == hit_ptr->get_read_id();
        if (!hits_correspond_to_correct_read) {
            fatal_error("Minimizer hits error: hit should be on read id ", read_id,
                ", but it is on read id ", hit_ptr->get_read_id());
        }

        const bool hits_correspond_to_correct_prg = prg_id == hit_ptr->get_prg_id();
        if (!hits_correspond_to_correct_prg) {
            fatal_error("Minimizer hits error: hit should be on PRG id ", prg_id,
                ", but it is on PRG id ", hit_ptr->get_prg_id());
        }
    }
}

// TODO: this should be a method of the Read class
// Add the node to the vector of nodes along read
// as well as its orientation (unless the previous node along
// the read was the same node and orientation)
// Store the hits on the read
void update_read_info_with_node_and_cluster(ReadPtr& read_ptr, const NodePtr& node_ptr,
    const std::set<MinimizerHitPtr, pComp>& cluster)
{
    read_ptr->add_hits(node_ptr, cluster);
}

void pangenome::Graph::add_hits_between_PRG_and_read(
    const std::shared_ptr<LocalPRG>&
        prg, // the prg from where this cluster of hits come
    const uint32_t read_id, // the read id from where this cluster of reads come - TODO:
                            // this can be derived from the cluster
    std::set<MinimizerHitPtr, pComp>& cluster // the cluster itself
)
{
    check_correct_hits(prg->id, read_id,
        cluster); // assure this cluster corresponds to the given prg and read

    // add and get the new read
    add_read(read_id);
    auto read_ptr = get_read(read_id);

    // add and get the new node
    add_node(prg);
    auto node_ptr = get_node(prg);

    // update the info
    update_node_info_with_this_read(node_ptr, read_ptr);
    update_read_info_with_node_and_cluster(read_ptr, node_ptr, cluster);
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
    // removes all instances of node n and references to it in reads
    for (const auto& r : n->reads) {
        r->remove_all_nodes_with_this_id(n->node_id);
    }

    auto it = nodes.find(n->node_id);
    if (it != nodes.end()) {
        it = nodes.erase(it);
    }
    return it;
}

// Remove read from each node which contains it
// and from the graph
// Remove nodes which no longer have any reads
void pangenome::Graph::remove_read(const uint32_t read_id)
{
    for (const auto& n : reads[read_id]->nodes) {
        auto nSharedPtr = n.lock();
        nSharedPtr->covg -= 1;
        nSharedPtr->reads.erase(reads[read_id]);
        if (nSharedPtr->covg == 0) {
            remove_node(nSharedPtr);
        }
    }
    reads.erase(read_id);
}

// Remove a single instance of a node from a read while iterating through the nodes of
// the read
std::vector<WeakNodePtr>::iterator pangenome::Graph::remove_node_from_read(
    std::vector<WeakNodePtr>::iterator node_it, ReadPtr read_ptr)
{
    auto node_ptr = node_it->lock();

    BOOST_LOG_TRIVIAL(debug) << "remove node " << node_ptr->node_id << " from read "
                             << read_ptr->id;
    // remove node from read
    node_it = read_ptr->remove_node_with_iterator(node_it);

    // remove read from node
    auto read_it = node_ptr->reads.find(read_ptr);
    if (read_it != node_ptr->reads.end())
        node_ptr->reads.erase(read_it);

    if (node_ptr->reads.size() == 0)
        remove_node(node_ptr);

    return node_it;
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

// Create a copy of the node with node_id and replace the old copy with
// the new one in each of the reads in reads_along_tig (by looking for the context of
// node_id)
void pangenome::Graph::split_node_by_reads(std::unordered_set<ReadPtr>& reads_along_tig,
    std::vector<uint_least32_t>& node_ids, const std::vector<bool>& node_orients,
    const uint_least32_t node_id)
{
    if (reads_along_tig.empty()) {
        return;
    }

    // replace the first instance of node_id which it finds on the read
    // (in the context of node_ids) with a new node
    while (nodes.find(next_id) != nodes.end()) {
        next_id++;
    }

    // define new node
    NodePtr n = std::make_shared<Node>(nodes[node_id]->prg, next_id,
        nodes[node_id]->kmer_prg_with_coverage.get_total_number_samples());
    n->covg -= 1;
    nodes[next_id] = n;

    // switch old node to new node in reads
    std::unordered_multiset<ReadPtr>::iterator rit;
    std::pair<uint32_t, uint32_t> pos;
    for (const auto& r : reads_along_tig) {
        // ignore if this node does not contain this read
        rit = nodes[node_id]->reads.find(r);
        if (rit == nodes[node_id]->reads.end()) {
            continue;
        }

        // find iterator to the node in the read
        pos = r->find_position(node_ids, node_orients);
        auto it = r->find_node_by_id(node_id);

        // replace the node in the read
        if (it != r->nodes.end()) {
            r->replace_node_with_iterator(it, n);
            nodes[node_id]->reads.erase(rit);
            nodes[node_id]->covg -= 1;
            if (nodes[node_id]->covg == 0) {
                remove_node(nodes[node_id]);
            }
            n->reads.insert(r);
            n->covg += 1;
        }
    }

    // replace node in tig
    for (uint32_t i = 0; i < node_ids.size(); ++i) {
        if (node_ids[i] == node_id) {
            node_ids[i] = next_id;
            break;
        }
    }
}

// For each node in pangraph, make a copy of the kmergraph and use the hits
// stored on each read containing the node to add coverage to this graph
void pangenome::Graph::add_hits_to_kmergraphs(const uint32_t& sample_id)
{
    for (const auto& node_entries : nodes) {
        Node& pangraph_node = *node_entries.second;
        const bool pangraph_node_has_a_valid_kmer_prg_with_coverage
            = (pangraph_node.kmer_prg_with_coverage.kmer_prg != nullptr)
            and (not pangraph_node.kmer_prg_with_coverage.kmer_prg->nodes.empty());
        if (!pangraph_node_has_a_valid_kmer_prg_with_coverage) {
            fatal_error(
                "Error adding hits to kmer graph: pangraph node does not have a "
                "valid Kmer PRG with coverage");
        }
        uint32_t num_hits[2] = { 0, 0 };

        // add hits
        for (const auto& read_ptr : pangraph_node.reads) {
            const Read& read = *read_ptr;

            auto hits = read.get_hits_as_unordered_map();
            for (const auto& minimizer_hit_ptr : hits.at(pangraph_node.prg_id)) {
                const auto& minimizer_hit = *minimizer_hit_ptr;

                const bool minimizer_hit_kmer_node_id_is_valid
                    = (minimizer_hit.get_kmer_node_id()
                          < pangraph_node.kmer_prg_with_coverage.kmer_prg->nodes.size())
                    && (pangraph_node.kmer_prg_with_coverage.kmer_prg
                            ->nodes[minimizer_hit.get_kmer_node_id()]
                        != nullptr);
                if (!minimizer_hit_kmer_node_id_is_valid) {
                    fatal_error("Error adding hits to kmer graph: minimizer hit "
                                "kmer node is invalid");
                }

                if (minimizer_hit.is_forward()) {
                    pangraph_node.kmer_prg_with_coverage.increment_forward_covg(
                        minimizer_hit.get_kmer_node_id(), sample_id);
                } else {
                    pangraph_node.kmer_prg_with_coverage.increment_reverse_covg(
                        minimizer_hit.get_kmer_node_id(), sample_id);
                }

                const auto covg { minimizer_hit.is_forward()
                        ? pangraph_node.kmer_prg_with_coverage.get_forward_covg(
                            minimizer_hit.get_kmer_node_id(), sample_id)
                        : pangraph_node.kmer_prg_with_coverage.get_reverse_covg(
                            minimizer_hit.get_kmer_node_id(), sample_id) };
                if (covg == 1000) {
                    BOOST_LOG_TRIVIAL(debug)
                        << "Adding hit " << minimizer_hit
                        << " resulted in high coverage on node "
                        << *pangraph_node.kmer_prg_with_coverage.kmer_prg
                                ->nodes[minimizer_hit.get_kmer_node_id()];
                }
                num_hits[minimizer_hit.is_forward()] += 1;
            }
        }

        BOOST_LOG_TRIVIAL(debug)
            << "Added " << num_hits[1] << " hits in the forward direction and "
            << num_hits[0] << " hits in the reverse";
        pangraph_node.kmer_prg_with_coverage.set_num_reads(pangraph_node.covg);
    }
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
        const auto reference_path = prg.get_valid_vcf_reference(vcf_reference_sequence);
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

void pangenome::Graph::save_mapped_read_strings(
    const fs::path& readfilepath, const fs::path& outdir, const int32_t buff)
{
    BOOST_LOG_TRIVIAL(debug) << "Save mapped read strings and coordinates";
    fs::ofstream outhandle;
    FastaqHandler readfile(readfilepath.string());
    uint32_t start, end;

    // for each node in pangraph, find overlaps and write to a file
    std::vector<std::vector<uint32_t>> read_overlap_coordinates;
    for (const auto& node_ptr : nodes) {
        BOOST_LOG_TRIVIAL(debug)
            << "Find coordinates for node " << node_ptr.second->name;
        node_ptr.second->get_read_overlap_coordinates(read_overlap_coordinates);

        const auto node_outpath { outdir / node_ptr.second->get_name()
            / (node_ptr.second->get_name() + ".reads.fa") };
        fs::create_directories(node_outpath.parent_path());
        outhandle.open(node_outpath);

        for (const auto& coord : read_overlap_coordinates) {
            readfile.get_nth_read(coord[0]);
            start = (uint32_t)std::max((int32_t)coord[1] - buff, 0);
            end = std::min(coord[2] + (uint32_t)buff, (uint32_t)readfile.read.length());

            const bool read_coordinates_are_valid = (coord[1] < coord[2])
                && (start <= coord[1]) && (start <= readfile.read.length())
                && (coord[2] <= readfile.read.length()) && (end >= coord[2])
                && (start < end);
            if (!read_coordinates_are_valid) {
                fatal_error("When saving mapped reads, read coordinates are not valid");
            }

            outhandle << ">" << readfile.name << " pandora: " << coord[0] << " "
                      << start << ":" << end;
            if (coord[3]) {
                outhandle << " + " << std::endl;
            } else {
                outhandle << " - " << std::endl;
            }
            outhandle << readfile.read.substr(start, end - start) << std::endl;
        }
        outhandle.close();
        read_overlap_coordinates.clear();
    }

    readfile.close();
}

std::ostream& pangenome::operator<<(std::ostream& out, const pangenome::Graph& m)
{
    for (const auto& n : m.nodes) {
        out << n.second->prg_id << std::endl;
        out << *n.second << std::endl;
    }
    for (const auto& n : m.reads) {
        out << *(n.second) << std::endl;
    }

    return out;
}
