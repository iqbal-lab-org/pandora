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
#include "localPRG.h"
#include "minihit.h"
#include "fastaq_handler.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace pangenome;

Graph::Graph() : next_id(0) {
    nodes.reserve(6000);
}

void Graph::reserve_num_reads(uint32_t &num_reads) {
    reads.reserve(num_reads);
}

void Graph::clear() {
    reads.clear();
    nodes.clear();
    samples.clear();
}

Graph::~Graph() {
    clear();
}

// Finds or creates a read with id read_id
ReadPtr Graph::get_read(const uint32_t &read_id) {
    auto it = reads.find(read_id);
    bool found = it != reads.end();
    if (not found) {
        auto read_ptr = make_shared<Read>(read_id);
        assert(read_ptr != nullptr);
        reads[read_id] = read_ptr;
        return read_ptr;
    }

    auto read_ptr = it->second;
    return read_ptr;
}

// Finds or creates a node with id node_id
// increments the coverage (number of times seen in a read)
// and records read_id as a read containing this node,
NodePtr Graph::add_coverage(ReadPtr &read_ptr,
                            const NodeId &node_id,
                            const uint32_t &prg_id,
                            const string &prg_name) {
    NodePtr node_ptr;
    auto it = nodes.find(node_id);
    bool found_node = it != nodes.end();
    if (not found_node) {
        node_ptr = make_shared<Node>(prg_id, node_id, prg_name);
        assert(node_ptr != nullptr);
        nodes[node_id] = node_ptr;
    } else {
        node_ptr = it->second;
        node_ptr->covg += 1;
    }
    node_ptr->reads.insert(read_ptr);
    assert(node_ptr->covg == node_ptr->reads.size());
    assert(node_id < numeric_limits<uint32_t>::max()
           or assert_msg("WARNING, prg_id reached max pangraph node size"));
    return node_ptr;
}

// Checks that all hits in the cluster are from the given prg and read
void check_correct_hits(const uint32_t prg_id,
                        const uint32_t read_id,
                        const set<MinimizerHitPtr, pComp> &cluster) {
    for (const auto &hit_ptr : cluster) {
        bool hits_correspond_to_correct_read = read_id == hit_ptr->read_id;
        assert(hits_correspond_to_correct_read);

        bool hits_correspond_to_correct_prg = prg_id == hit_ptr->prg_id;
        assert(hits_correspond_to_correct_prg);
    }
}

// Add the node to the vector of nodes along read
// as well as its orientation (unless the previous node along
// the read was the same node and orientation)
// Store the hits on the read
void record_read_info(ReadPtr &read_ptr,
                      const NodePtr &node_ptr,
                      set<MinimizerHitPtr, pComp> &cluster) {
    assert(read_ptr != nullptr);
    read_ptr->add_hits(node_ptr->node_id, cluster);
    bool orientation = !cluster.empty() and (*cluster.begin())->strand;
    if (read_ptr->nodes.empty()
	or node_ptr != read_ptr->nodes.back()
	or orientation != read_ptr->node_orientations.back()
	//or we think there really are 2 copies of gene
	) {
	    read_ptr->nodes.push_back(node_ptr);
        read_ptr->node_orientations.push_back(orientation);
    }
}


// Add a node corresponding to a cluster of hits against a given localPRG from a read
void Graph::add_node(const uint32_t prg_id,
                     const string &prg_name,
                     const uint32_t read_id,
                     set<MinimizerHitPtr, pComp> &cluster) {
    check_correct_hits(prg_id, read_id, cluster);
    auto read_ptr = get_read(read_id);
    assert(read_ptr != nullptr);

    // add new node if it doesn't exist
    const auto &node_id = prg_id;
    auto node_ptr = add_coverage(read_ptr, node_id, prg_id, prg_name);
    assert(node_ptr != nullptr);

    record_read_info(read_ptr, node_ptr, cluster);
}

// Add a node corresponding to an instance of a localPRG found in a sample
void Graph::add_node(const uint32_t prg_id, const string &prg_name, const string &sample_name,
                     const vector<KmerNodePtr> &kmp, const LocalPRG *prg) {
    // add new node if it doesn't exist
    NodePtr n;
    auto it = nodes.find(prg_id);
    if (it == nodes.end()) {
        //cout << "add node " << *n << endl;
        n = make_shared<Node>(prg_id, prg_id, prg_name);
        n->kmer_prg = prg->kmer_prg;
        nodes[prg_id] = n;
        //nodes[prg_id] = make_shared<Node>(prg_id, prg_id, prg_name);
        //it = nodes.find(prg_id);
    } else {
        n = it->second;
        n->covg += 1;
        //cout << "node " << *n << " already existed " << endl;
    }

    // add a new sample if it doesn't exist
    SamplePtr s;
    auto sit = samples.find(sample_name);
    if (sit == samples.end()) {
        //cout << "new sample " << sample_name << endl;
        s = make_shared<Sample>(sample_name);
        samples[sample_name] = s;
        //sit = samples.find(sample_name);
    } else {
        s = sit->second;
    }
    //cout << "sample " << sample_name  << " already existed " << endl;
    s->add_path(prg_id, kmp);
    n->samples.insert(s);
    n->add_path(kmp);
}

// Remove the node n, and all references to it
unordered_map<uint32_t, NodePtr>::iterator Graph::remove_node(NodePtr n) {
    //cout << "Remove graph node " << *n << endl;
    // removes all instances of node n and references to it in reads
    for (auto r : n->reads) {
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
void Graph::remove_read(const uint32_t read_id) {
    for (auto n : reads[read_id]->nodes) {
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
vector<NodePtr>::iterator Graph::remove_node_from_read(vector<NodePtr>::iterator node_it, ReadPtr read_ptr) {

    cout << "remove node " << (*node_it)->node_id << " from read " << read_ptr->id << endl;
    NodePtr node_ptr = *node_it;

    // remove node from read
    node_it = read_ptr->remove_node(node_it);

    // remove read from node
    auto read_it = node_ptr->reads.find(read_ptr);
    if (read_it != node_ptr->reads.end())
        node_ptr->reads.erase(read_it);

    if (node_ptr->reads.size()==0)
        remove_node(node_ptr);

    return node_it;
}

// remove the all instances of the pattern of nodes/orienations from graph

// Remove nodes with covg <= thresh from graph
void Graph::remove_low_covg_nodes(const uint32_t &thresh) {
    cout << now() << "Remove nodes with covg <= " << thresh << endl;
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
    cout << now() << "Pangraph now has " << nodes.size() << " nodes" << endl;
}

// Create a copy of the node with node_id and replace the old copy with
// the new one in each of the reads in reads_along_tig (by looking for the context of node_id)
void Graph::split_node_by_reads(unordered_set<ReadPtr> &reads_along_tig, vector<uint_least32_t> &node_ids,
                                const vector<bool> &node_orients, const uint_least32_t node_id) {
    if (reads_along_tig.empty()) {
        return;
    }

    // replace the first instance of node_id which it finds on the read
    // (in the context of node_ids) with a new node
    while (nodes.find(next_id) != nodes.end()) {
        next_id++;
        assert(next_id < numeric_limits<uint32_t>::max() ||
               assert_msg("WARNING, next_id reached max pangraph node size"));
    }

    // define new node
    NodePtr n = make_shared<Node>(nodes[node_id]->prg_id, next_id, nodes[node_id]->name);
    n->covg -= 1;
    nodes[next_id] = n;

    // switch old node to new node in reads
    vector<NodePtr>::iterator it;
    unordered_multiset<ReadPtr>::iterator rit;
    pair<uint32_t, uint32_t> pos;
    for (auto r : reads_along_tig) {
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
    for (const auto i : node_path_ids)
    {
        for (const auto r : nodes[n]->reads)
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

// For each node in pangraph, make a copy of the kmergraph and use the hits
// stored on each read containing the node to add coverage to this graph
void Graph::add_hits_to_kmergraphs(const vector<LocalPRG *> &prgs) {
    uint32_t num_hits[2];
    for (auto pnode : nodes) {
        // copy kmergraph
        pnode.second->kmer_prg = prgs[pnode.second->prg_id]->kmer_prg;
        assert(pnode.second->kmer_prg == prgs[pnode.second->prg_id]->kmer_prg);
        num_hits[0] = 0;
        num_hits[1] = 0;

        // add hits
        for (auto read : pnode.second->reads) {
            for (auto mh = read->hits[pnode.second->prg_id].begin();
                 mh != read->hits[pnode.second->prg_id].end(); ++mh) {
                //bool added = false;
                // update the covg in the kmer_prg
                //cout << "pnode " << pnode.second->prg_id << " knode " << (*mh)->knode_id << " strand " << (*mh)->strand << " updated from " << pnode.second->kmer_prg.nodes[(*mh)->knode_id]->covg[(*mh)->strand];
                assert((*mh)->knode_id < pnode.second->kmer_prg.nodes.size() and pnode.second->kmer_prg.nodes[(*mh)->knode_id]!=nullptr);
                pnode.second->kmer_prg.nodes[(*mh)->knode_id]->covg[(*mh)->strand] += 1;
                if (pnode.second->kmer_prg.nodes[(*mh)->knode_id]->covg[(*mh)->strand] == 1000){
                    cout << "Adding hit " << **mh << " resulted in high coverage on node " << *pnode.second->kmer_prg.nodes[(*mh)->knode_id] << endl;
                }
                //cout << " to " << pnode.second->kmer_prg.nodes[(*mh)->knode_id]->covg[(*mh)->strand] << endl;
                num_hits[(*mh)->strand] += 1;
            }
        }
        cout << now() << "Added " << num_hits[1] << " hits in the forward direction and " << num_hits[0]
             << " hits in the reverse" << endl;
        pnode.second->kmer_prg.num_reads = pnode.second->covg;
    }
}

same_prg_id::same_prg_id(const NodePtr &p) : q(p->prg_id) {};

bool same_prg_id::operator()(const pair<uint32_t, NodePtr> &n) const { return (n.second->prg_id == q); }

bool Graph::operator==(const Graph &y) const {
    // false if have different numbers of nodes
    /*if (y.nodes.size() != nodes.size()) {
        cout << "different num nodes " << nodes.size() << "!=" << y.nodes.size() << endl;
        return false;
    }*/

    // false if have different nodes
    for (const auto c: nodes) {
        // if node id doesn't exist
        auto it = find_if(y.nodes.begin(), y.nodes.end(), same_prg_id(c.second));
        if (it == y.nodes.end()) {
            cout << "can't find node " << c.first << endl;
            return false;
        }
    }
    for (const auto c: y.nodes) {
        // if node id doesn't exist
        auto it = find_if(nodes.begin(), nodes.end(), same_prg_id(c.second));
        if (it == nodes.end()) {
            cout << "can't find node " << c.first << endl;
            return false;
        }
    }

    // otherwise is true
    return true;
}

bool Graph::operator!=(const Graph &y) const {
    return !(*this == y);
}

// Saves a presence/absence/copynumber matrix for each node and each sample
void Graph::save_matrix(const string &filepath) {
    // write a presence/absence matrix for samples and nodes
    ofstream handle;
    handle.open(filepath);

    // save header line with sample names
    for (auto s : samples) {
        handle << "\t" << s.second->name;
    }
    handle << endl;

    // for each node, save number of each sample
    for (auto n : nodes) {
        handle << n.second->name;
        for (auto s : samples) {
            if (s.second->paths.find(n.second->node_id) == s.second->paths.end()) {
                handle << "\t0";
            } else {
                handle << "\t" << s.second->paths[n.second->node_id].size();
            }
        }
        handle << endl;
    }
}

void Graph::save_mapped_read_strings(const string& readfilepath, const string& outdir, const int32_t buff){
    cout << now() << "Save mapped read strings and coordinates" << endl;
    ofstream outhandle;
    FastaqHandler readfile(readfilepath);
    uint32_t start, end;

    // for each node in pangraph, find overlaps and write to a file
    vector<vector<uint32_t>> read_overlap_coordinates;
    for (auto node_ptr : nodes)
    {
	    cout << "Find coordinates for node " << node_ptr.second->name;
        node_ptr.second->get_read_overlap_coordinates(read_overlap_coordinates);
	    cout << "." << endl;
        make_dir(outdir + "/" + node_ptr.second->get_name());
        outhandle.open(outdir + "/" + node_ptr.second->get_name() + "/" + node_ptr.second->get_name() + ".reads.fa");
        for (auto coord : read_overlap_coordinates){
            readfile.get_id(coord[0]);
            start = (uint32_t) max((int32_t)coord[1]-buff, 0);
            end = min(coord[2]+(uint32_t)buff, (uint32_t)readfile.read.length());
            outhandle << ">" << readfile.name << " pandora: " << coord[0] << " " << start << ":" << end;
            if (coord[3] == true)
                outhandle << " + " << endl;
            else
                outhandle << " - " << endl;
	    assert(coord[1] < coord[2]);
	    assert(start <= coord[1]);
	    assert(start <= readfile.read.length());
	    assert(coord[2] <= readfile.read.length());
            assert(end >= coord[2]);
	    assert(start < end);
            outhandle << readfile.read.substr(start, end - start) << endl;
        }
        outhandle.close();
        read_overlap_coordinates.clear();
    }

    readfile.close();
}

void Graph::save_kmergraph_coverages(const string &outdir, const string &sample_name){
    cout << now() << "Save kmergraph coverages for sample " << sample_name << endl;
    make_dir(outdir + "/coverages");
    for (auto n : nodes){
        string node_file = outdir + "/coverages/" + n.second->name + ".csv";
        if ( !boost::filesystem::exists(node_file))
        {
            ofstream handle;
            handle.open(node_file);
            assert (!handle.fail() or assert_msg("Could not open file " << node_file));
            handle << "sample";
            for (auto m : n.second->kmer_prg.nodes)
                handle << "\t" << m->id;
            handle << endl;
            handle.close();
        }
        n.second->kmer_prg.append_coverages_to_file(node_file, sample_name);
    }
}

std::ostream &pangenome::operator<<(std::ostream &out, pangenome::Graph const &m) {
    //cout << "printing pangraph" << endl;
    /*for (const auto &n : m.nodes) {
        cout << n.second->prg_id << endl;
    }*/
    for (const auto &n : m.reads) {
        cout << *(n.second) << endl;
    }

    return out;
}
