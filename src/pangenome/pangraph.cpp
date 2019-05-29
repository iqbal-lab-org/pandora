#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <memory>
#include <vector>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <boost/filesystem.hpp>

#include "utils.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "pangenome/pansample.h"
#include "minihit.h"
#include "fastaq_handler.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace pangenome;

pangenome::Graph::Graph() : next_id(0) {
    nodes.reserve(6000);
}

void pangenome::Graph::clear() {
    reads.clear();
    nodes.clear();
    samples.clear();
}

pangenome::Graph::~Graph() {
    clear();
}

// Finds or creates a read with id read_id
ReadPtr pangenome::Graph::get_read(const uint32_t &read_id) {
    auto it = reads.find(read_id);
    bool found = it != reads.end();
    if (not found) {
        auto read_ptr = std::make_shared<Read>(read_id);
        assert(read_ptr != nullptr);
        reads[read_id] = read_ptr;
        return read_ptr;
    }

    auto read_ptr = it->second;
    return read_ptr;
}

// Finds or creates a node with id node_id
NodePtr pangenome::Graph::get_node(const NodeId &node_id,
                        const uint32_t &prg_id,
                        const std::string &prg_name){
    NodePtr node_ptr;
    auto it = nodes.find(node_id);
    bool found_node = it != nodes.end();
    if (not found_node) {
        node_ptr = std::make_shared<Node>(prg_id, node_id, prg_name);
        assert(node_ptr != nullptr);
        nodes[node_id] = node_ptr;
    } else {
        node_ptr = it->second;
        node_ptr->covg += 1;
    }
    assert(node_id < std::numeric_limits<uint32_t>::max()
           or assert_msg("WARNING, prg_id reached max pangraph node size"));
    return node_ptr;
}

// increments the coverage (number of times seen in a read)
// and records read_id as a read containing this node,
NodePtr pangenome::Graph::add_coverage(ReadPtr &read_ptr,
                            const NodeId &node_id,
                            const uint32_t &prg_id,
                            const std::string &prg_name) {
    NodePtr node_ptr = get_node(node_id, prg_id, prg_name);
    node_ptr->reads.insert(read_ptr);
    assert(node_ptr->covg == node_ptr->reads.size());

    return node_ptr;
}

// Checks that all hits in the cluster are from the given prg and read
void check_correct_hits(const uint32_t prg_id,
                        const uint32_t read_id,
                        const std::set<MinimizerHitPtr, pComp> &cluster) {
    for (const auto &hit_ptr : cluster) {
        bool hits_correspond_to_correct_read = read_id == hit_ptr->get_read_id();
        assert(hits_correspond_to_correct_read);

        bool hits_correspond_to_correct_prg = prg_id == hit_ptr->get_prg_id();
        assert(hits_correspond_to_correct_prg);
    }
}

// Add the node to the vector of nodes along read
// as well as its orientation (unless the previous node along
// the read was the same node and orientation)
// Store the hits on the read
void record_read_info(ReadPtr &read_ptr,
                      const NodePtr &node_ptr,
                      std::set<MinimizerHitPtr, pComp> &cluster) {
    assert(read_ptr != nullptr);
    read_ptr->add_hits(node_ptr->node_id, cluster);
    bool orientation = !cluster.empty() and (*cluster.begin())->is_forward();
    if (read_ptr->get_nodes().empty()
        or node_ptr != read_ptr->get_nodes().back()
        or orientation != read_ptr->node_orientations.back()
        //or we think there really are 2 copies of gene
            ) {
        read_ptr->add_node(node_ptr);
        read_ptr->add_orientation(orientation);
    }
}


// Add a node corresponding to a cluster of hits against a given localPRG from a read
void pangenome::Graph::add_node(const uint32_t prg_id, //the prg from where this cluster of hits come //TODO: this can be inferred from cluster
                     const std::string &prg_name, //the prg name
                     const uint32_t read_id, //the read id from where this cluster of reads come //TODO: this can be inferred from cluster
                     std::set<MinimizerHitPtr, pComp> &cluster //the cluster itself
                     ) {
    check_correct_hits(prg_id, read_id, cluster); //assure this cluster corresponds to the given prg and read
    auto read_ptr = get_read(read_id);
    assert(read_ptr != nullptr);

    // add new node if it doesn't exist
    const auto &node_id = prg_id;
    auto node_ptr = add_coverage(read_ptr, node_id, prg_id, prg_name);
    assert(node_ptr != nullptr);

    record_read_info(read_ptr, node_ptr, cluster);
}

SamplePtr pangenome::Graph::get_sample(const std::string &sample_name, const uint32_t &sample_id) {
    SamplePtr s;
    auto sit = samples.find(sample_name);
    if (sit == samples.end()) {
        BOOST_LOG_TRIVIAL(debug) << sample_name << " is a new sample in pangraph";
        s = std::make_shared<Sample>(sample_name, sample_id);
        samples[sample_name] = s;
    } else {
        BOOST_LOG_TRIVIAL(debug) << "found sample " << sample_name;
        s = sit->second;
    }
    return s;
}

// Add a node corresponding to an instance of a localPRG found in a sample
void pangenome::Graph::add_node(const uint32_t prg_id, const std::string &prg_name, const std::string &sample_name,
                     const uint32_t &sample_id, const std::shared_ptr<LocalPRG> &prg,
                     const std::vector<KmerNodePtr> &kmp) {
    // add new node if it doesn't exist
    NodePtr n = get_node(prg_id, prg_id, prg_name);

    // add a new sample if it doesn't exist
    SamplePtr s = get_sample(sample_name, sample_id);
    s->add_path(prg_id, kmp);
    if (std::find(n->samples.begin(), n->samples.end(), s) == n->samples.end())
        n->samples.push_back(s);
}

// Remove the node n, and all references to it
std::unordered_map<uint32_t, NodePtr>::iterator pangenome::Graph::remove_node(NodePtr n) {
    //cout << "Remove graph node " << *n << endl;
    // removes all instances of node n and references to it in reads
    for (const auto &r : n->reads) {
        r->remove_node(n);
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
void pangenome::Graph::remove_read(const uint32_t read_id) {
    for (const auto &n : reads[read_id]->nodes) {
        //cout << "looking at read node " << n->node_id;
        n->covg -= 1;
        n->reads.erase(reads[read_id]);
        if (n->covg == 0) {
            remove_node(n);
        }
    }
    reads.erase(read_id);
}

// Remove a single instance of a node from a read while iterating through the nodes of the read
std::vector<NodePtr>::iterator pangenome::Graph::remove_node_from_read(std::vector<NodePtr>::iterator node_it, ReadPtr read_ptr) {

    BOOST_LOG_TRIVIAL(debug) << "remove node " << (*node_it)->node_id << " from read " << read_ptr->id;
    NodePtr node_ptr = *node_it;

    // remove node from read
    node_it = read_ptr->remove_node(node_it);

    // remove read from node
    auto read_it = node_ptr->reads.find(read_ptr);
    if (read_it != node_ptr->reads.end())
        node_ptr->reads.erase(read_it);

    if (node_ptr->reads.size() == 0)
        remove_node(node_ptr);

    return node_it;
}

// remove the all instances of the pattern of nodes/orienations from graph

// Remove nodes with covg <= thresh from graph
void pangenome::Graph::remove_low_covg_nodes(const uint32_t &thresh) {
    BOOST_LOG_TRIVIAL(debug) << "Remove nodes with covg <= " << thresh << std::endl;
    for (auto it = nodes.begin(); it != nodes.end();) {
        //cout << "look at node " << *(it->second) << endl;
        if (it->second->covg <= thresh) {
            //cout << "delete node " << it->second->name;
            it = remove_node(it->second);
            //cout << " so pangraph now has " << nodes.size() << " nodes" << endl;
        } else {
            ++it;
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "Pangraph now has " << nodes.size() << " nodes";
}

// Create a copy of the node with node_id and replace the old copy with
// the new one in each of the reads in reads_along_tig (by looking for the context of node_id)
void pangenome::Graph::split_node_by_reads(std::unordered_set<ReadPtr> &reads_along_tig, std::vector<uint_least32_t> &node_ids,
                                const std::vector<bool> &node_orients, const uint_least32_t node_id) {
    if (reads_along_tig.empty()) {
        return;
    }

    // replace the first instance of node_id which it finds on the read
    // (in the context of node_ids) with a new node
    while (nodes.find(next_id) != nodes.end()) {
        next_id++;
        assert(next_id < std::numeric_limits<uint32_t>::max() ||
               assert_msg("WARNING, next_id reached max pangraph node size"));
    }

    // define new node
    NodePtr n = std::make_shared<Node>(nodes[node_id]->prg_id, next_id, nodes[node_id]->name);
    n->covg -= 1;
    nodes[next_id] = n;

    // switch old node to new node in reads
    std::vector<NodePtr>::iterator it;
    std::unordered_multiset<ReadPtr>::iterator rit;
    std::pair<uint32_t, uint32_t> pos;
    for (const auto &r : reads_along_tig) {
        // ignore if this node does not contain this read
        rit = nodes[node_id]->reads.find(r);
        if (rit == nodes[node_id]->reads.end()) {
            continue;
        }

        //find iterator to the node in the read
        pos = r->find_position(node_ids, node_orients);
        it = std::find(r->nodes.begin() + pos.first, r->nodes.end(), nodes[node_id]);

        // replace the node in the read
        if (it != r->nodes.end()) {
            //cout << "replace node " << (*it)->node_id << " in this read" << endl;
            //cout << "read was " << *r << endl;
            r->replace_node(it, n);
            //cout << "read is now " << *r << endl;
            nodes[node_id]->reads.erase(rit);
            nodes[node_id]->covg -= 1;
            if (nodes[node_id]->covg == 0) {
                remove_node(nodes[node_id]);
            }
            n->reads.insert(r);
            n->covg += 1;
            //} else {
            //    cout  << "read does not contain this node in context of the tig" << endl;
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


/*Graph::unordered_set<ReadPtr> find_reads_on_node_path(const vector<uint16_t> node_path_ids,
                                                      const vector<bool> node_path_orients)
{
    unordered_set<ReadPtr> reads_on_node_path;

    // collect all reads on nodes
    for (const auto &i : node_path_ids)
    {
        for (const auto &r : nodes[n]->reads)
        {
            reads_on_node_path.insert(r);
        }
    }

    // remove reads which deviate from path
    for (auto rit = reads_on_node_path.begin(); rit != reads_on_node_path.end();)
    {
        if ((*rit)->find_position(node_path_ids, node_path_orients) == std::numeric_limits<uint>::max())
        {
            rit = reads_on_node_path.erase(rit);
        } else {
            rit++;
        }
    }
    return reads_on_node_path;
}*/


void pangenome::Graph::setup_kmergraphs(const std::vector<std::shared_ptr<LocalPRG>> &prgs,
                             const uint64_t &total_number_samples) {
    for (const auto &node_entries: this->nodes) {
        Node &pangraph_node = *node_entries.second;

        if (not pangraph_node.kmer_prg.nodes.empty())
            continue;

        BOOST_LOG_TRIVIAL(debug) << "setup kmergraphs for node " << pangraph_node.get_name();
        assert(pangraph_node.prg_id < prgs.size());
        pangraph_node.kmer_prg = prgs[pangraph_node.prg_id]->kmer_prg; //TODO: optimize this
        pangraph_node.kmer_prg.setup_coverages(total_number_samples);
    }
}

// For each node in pangraph, make a copy of the kmergraph and use the hits
// stored on each read containing the node to add coverage to this graph
void pangenome::Graph::add_hits_to_kmergraphs(const std::vector<std::shared_ptr<LocalPRG>> &prgs,
                                   const uint32_t &sample_id) {
    for (const auto &node_entries: nodes) {
        Node &pangraph_node = *node_entries.second;
        assert(not pangraph_node.kmer_prg.nodes.empty());

        uint32_t num_hits[2] = {0, 0};

        // add hits
        for (const auto &read_ptr: pangraph_node.reads) {
            const Read &read = *read_ptr;

            auto hits = read.getHits();
            for (const auto &minimizer_hit_ptr : hits.at(pangraph_node.prg_id)) {
                const auto &minimizer_hit = *minimizer_hit_ptr;

                assert(minimizer_hit.get_kmer_node_id() < pangraph_node.kmer_prg.nodes.size());
                assert(pangraph_node.kmer_prg.nodes[minimizer_hit.get_kmer_node_id()] != nullptr);

                auto &kmer_node = *pangraph_node.kmer_prg.nodes[minimizer_hit.get_kmer_node_id()];
                kmer_node.increment_covg(minimizer_hit.is_forward(), sample_id);

                if (pangraph_node.kmer_prg.nodes[minimizer_hit.get_kmer_node_id()]->get_covg(minimizer_hit.is_forward(), sample_id) ==
                    1000) {
                    BOOST_LOG_TRIVIAL(debug) << "Adding hit " << minimizer_hit
                                             << " resulted in high coverage on node "
                                             << *pangraph_node.kmer_prg.nodes[minimizer_hit.get_kmer_node_id()];
                }
                num_hits[minimizer_hit.is_forward()] += 1;
            }
        }

        BOOST_LOG_TRIVIAL(debug) << "Added " << num_hits[1] << " hits in the forward direction and "
                                 << num_hits[0]
                                 << " hits in the reverse";
        pangraph_node.kmer_prg.num_reads = pangraph_node.covg;
    }
}

// For each node in reference pangraph, copy the coverages over to sample_id in this pangraph
void pangenome::Graph::copy_coverages_to_kmergraphs(const Graph &ref_pangraph, const uint32_t &sample_id) {
    const uint32_t ref_sample_id = 0;
    for (const auto &ref_node_entry : ref_pangraph.nodes){
        const Node &ref_node = *ref_node_entry.second;
        assert(nodes.find(ref_node.node_id) != nodes.end());
        Node &pangraph_node = *nodes[ref_node.node_id];

        for (auto & kmergraph_node_ptr : pangraph_node.kmer_prg.nodes){
            const auto &knode_id = kmergraph_node_ptr->id;
            assert(knode_id < ref_node.kmer_prg.nodes.size());
            kmergraph_node_ptr->set_covg(ref_node.kmer_prg.nodes[knode_id]->get_covg(0, ref_sample_id), 0, sample_id);
            kmergraph_node_ptr->set_covg(ref_node.kmer_prg.nodes[knode_id]->get_covg(1, ref_sample_id), 1, sample_id);
        }
    }
}


std::vector<LocalNodePtr>
pangenome::Graph::infer_node_vcf_reference_path(const Node &node, const std::shared_ptr<LocalPRG> &prg_ptr, const uint32_t &w,
                                     const std::unordered_map<std::string, std::string> &vcf_refs) {
    BOOST_LOG_TRIVIAL(info) << "Infer VCF reference path";
    const auto &prg = *prg_ptr;
    if (vcf_refs.find(prg.name) != vcf_refs.end()){
        const auto &vcf_reference_sequence = vcf_refs.at(prg.name);
        const auto reference_path = prg.get_valid_vcf_reference(vcf_reference_sequence);
        if (!reference_path.empty())
            return reference_path;
    }
    return get_node_closest_vcf_reference(node, w, prg);
}

std::vector<LocalNodePtr>
pangenome::Graph::get_node_closest_vcf_reference(const Node &node, const uint32_t &w, const LocalPRG &prg) {
    auto kmer_graph = prg.kmer_prg;


    for (const auto &sample_entry: this->samples) {
        const auto &sample = sample_entry.second;
        if (sample->paths.find(node.prg_id) == sample->paths.end()){
            BOOST_LOG_TRIVIAL(debug) << "Could not find sample path for sample " << sample->name << " and prg " << node.prg_id;
            continue;
        }

        BOOST_LOG_TRIVIAL(debug) << "Increment coverages for path for sample " << sample->name << " and prg " << node.prg_id;
        const auto &sample_paths = sample->paths.at(node.prg_id);
        for (const auto &sample_path : sample_paths) {
            for (uint32_t i = 0; i != sample_path.size(); ++i) {
                assert(sample_path[i]->id < kmer_graph.nodes.size()
                       and kmer_graph.nodes[sample_path[i]->id] != nullptr);
                kmer_graph.nodes[sample_path[i]->id]->increment_covg(0, 0);
                kmer_graph.nodes[sample_path[i]->id]->increment_covg(1, 0);
            }
        }
    }

    kmer_graph.discover_k();
    kmer_graph.num_reads = node.covg;

    std::vector<KmerNodePtr> kmer_path;
    kmer_graph.find_lin_max_path(kmer_path, 0);
    if (!kmer_path.empty()) {
        auto reference_path = prg.localnode_path_from_kmernode_path(kmer_path, w);
        BOOST_LOG_TRIVIAL(debug) << "Found reference path to return";
        return reference_path;
    } else {
        BOOST_LOG_TRIVIAL(debug) << "Could not infer reference path, so use top path";
        return prg.prg.top_path();
    }
}

same_prg_id::same_prg_id(const NodePtr &p) : q(p->prg_id) {};

bool same_prg_id::operator()(const std::pair<uint32_t, NodePtr> &n) const { return (n.second->prg_id == q); }

bool pangenome::Graph::operator==(const Graph &y) const {
    // false if have different numbers of nodes
    /*if (y.nodes.size() != nodes.size()) {
        cout << "different num nodes " << nodes.size() << "!=" << y.nodes.size() << endl;
        return false;
    }*/

    // false if have different nodes
    for (const auto &c: nodes) {
        // if node id doesn't exist
        auto it = find_if(y.nodes.begin(), y.nodes.end(), same_prg_id(c.second));
        if (it == y.nodes.end()) {
            BOOST_LOG_TRIVIAL(warning) << "can't find node " << c.first;
            return false;
        }
    }
    for (const auto &c: y.nodes) {
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

bool pangenome::Graph::operator!=(const Graph &y) const {
    return !(*this == y);
}

// Saves a presence/absence/copynumber matrix for each node and each sample
void pangenome::Graph::save_matrix(const std::string &filepath, const std::vector<std::string> &sample_names) {
    // write a presence/absence matrix for samples and nodes
    std::ofstream handle;
    handle.open(filepath);

    // save header line with sample names
    for (const auto &name : sample_names) {
        handle << "\t" << name;
    }
    handle << std::endl;

    // for each node, save number of each sample
    for (const auto &n : nodes) {
        handle << n.second->name;
        for (const auto &name : sample_names) {
            if (samples.find(name) == samples.end()
                or samples[name]->paths.find(n.second->node_id) == samples[name]->paths.end()) {
                handle << "\t0";
            } else {
                handle << "\t" << samples[name]->paths[n.second->node_id].size();
            }
        }
        handle << std::endl;
    }
}

void pangenome::Graph::save_mapped_read_strings(const std::string &readfilepath, const std::string &outdir, const int32_t buff) {
    BOOST_LOG_TRIVIAL(debug) << "Save mapped read strings and coordinates";
    std::ofstream outhandle;
    FastaqHandler readfile(readfilepath);
    uint32_t start, end;

    // for each node in pangraph, find overlaps and write to a file
    std::vector<std::vector<uint32_t>> read_overlap_coordinates;
    for (const auto &node_ptr : nodes) {
        std::cout << "Find coordinates for node " << node_ptr.second->name;
        node_ptr.second->get_read_overlap_coordinates(read_overlap_coordinates);
        std::cout << "." << std::endl;
        fs::create_directories(outdir + "/" + node_ptr.second->get_name());
        outhandle.open(outdir + "/" + node_ptr.second->get_name() + "/" + node_ptr.second->get_name() + ".reads.fa");
        for (const auto &coord : read_overlap_coordinates) {
            readfile.get_id(coord[0]);
            start = (uint32_t) std::max((int32_t) coord[1] - buff, 0);
            end = std::min(coord[2] + (uint32_t) buff, (uint32_t) readfile.read.length());
            outhandle << ">" << readfile.name << " pandora: " << coord[0] << " " << start << ":" << end;
            if (coord[3])
                outhandle << " + " << std::endl;
            else
                outhandle << " - " << std::endl;
            assert(coord[1] < coord[2]);
            assert(start <= coord[1]);
            assert(start <= readfile.read.length());
            assert(coord[2] <= readfile.read.length());
            assert(end >= coord[2]);
            assert(start < end);
            outhandle << readfile.read.substr(start, end - start) << std::endl;
        }
        outhandle.close();
        read_overlap_coordinates.clear();
    }

    readfile.close();
}

std::ostream &pangenome::operator<<(std::ostream &out, const pangenome::Graph &m) {
    //cout << "printing pangraph" << endl;
    for (const auto &n : m.nodes) {
        std::cout << n.second->prg_id << std::endl;
        std::cout << *n.second << std::endl;
    }
    for (const auto &n : m.reads) {
        std::cout << *(n.second) << std::endl;
    }

    return out;
}
