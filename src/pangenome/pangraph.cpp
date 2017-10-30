#include <iostream>
#include <map>
#include <unordered_set>
#include <vector>
#include <fstream>
#include "utils.h"
#include "pangenome/pangraph.h"
#include "pangenome/panread.h"
#include "pangenome/pansample.h"
#include <cassert>
#include "minihit.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace pangenome;

Graph::Graph() : next_id(0) {}

void Graph::clear() {
    next_id = 0;
    reads.clear();
    nodes.clear();
    samples.clear();
}

Graph::~Graph() {
    clear();
}

// add a node corresponding to a cluster of hits against a given localPRG from a read
void Graph::add_node(const uint32_t prg_id, const string prg_name, const uint32_t read_id,
                        const set<MinimizerHitPtr, pComp> &cluster) {
    // check sensible things in new cluster
    // NB if cluster is empty, handle by adding in 0 orientation
    cout << now() << "Add node" << endl;
    for (const auto &it : cluster) {
        assert(read_id == it->read_id); // the hits should correspond to the read we are saying...
        assert(prg_id == it->prg_id); // and the prg we are saying
    }
    cout << now() << "Checked cluster was sensible" << endl;

    // add new node if it doesn't exist
    NodePtr n;
    auto it = nodes.find(prg_id);
    if (it == nodes.end()) {
        n = make_shared<Node>(prg_id, prg_id, prg_name);
        cout << "add node " << *n << endl;
        nodes[prg_id] = n;
    } else {
        n = it->second;
        it->second->covg += 1;
        cout << "node " << *n << " already existed " << endl;
    }

    // add a new read if it doesn't exist
    auto rit = reads.find(read_id);
    if (rit == reads.end()) {
        cout << "new read " << read_id << endl;
        ReadPtr r(make_shared<Read>(read_id));
        reads[read_id] = r;

        n->reads.insert(r);
        r->add_hits(prg_id, cluster);
        r->nodes.push_back(n);
        if(!cluster.empty())
        {
            r->node_orientations.push_back((*cluster.begin())->strand);
        } else {
            r->node_orientations.push_back(0);
        }
    } else {
        cout << "read " << read_id  << " already existed " << endl;
        n->reads.insert(rit->second);
        rit->second->add_hits(prg_id, cluster);
        rit->second->nodes.push_back(n);
        if(!cluster.empty())
        {
            rit->second->node_orientations.push_back((*cluster.begin())->strand);
        } else {
            rit->second->node_orientations.push_back(0);
        }
    }

    assert(n->covg == n->reads.size());
}

// add a node corresponding to an instance of a localPRG found in a sample
void Graph::add_node(const uint32_t prg_id, const string &prg_name, const string &sample_name,
                        const vector<KmerNodePtr> &kmp, const LocalPRG *prg) {
    // add new node if it doesn't exist
    NodePtr n;
    auto it = nodes.find(prg_id);
    if (it == nodes.end()) {
        n = make_shared<Node>(prg_id, prg_id, prg_name);
        n->kmer_prg = prg->kmer_prg;
        //cout << "add node " << *n << endl;
        nodes[prg_id] = n;
    } else {
        n = it->second;
        it->second->covg += 1;
        //cout << "node " << *n << " already existed " << endl;
    }

    // add a new sample if it doesn't exist
    auto sit = samples.find(sample_name);
    if (sit == samples.end()) {
        //cout << "new sample " << sample_name << endl;
        SamplePtr s(make_shared<Sample>(sample_name));
        s->add_path(prg_id, kmp);
        samples[sample_name] = s;
        n->samples.insert(s);
        n->add_path(kmp);
    } else {
        //cout << "sample " << sample_name  << " already existed " << endl;
        sit->second->add_path(prg_id, kmp);
        n->samples.insert(sit->second);
        n->add_path(kmp);
    }

}

// remove the node n, and all references to it
map<uint32_t, NodePtr>::iterator Graph::remove_node(NodePtr n)
{
    //cout << "Remove graph node " << *n << endl;
    // removes all instances of node n and references to it in reads
    for (auto r : n->reads)
    {
        r->remove_node(n);
    }

    auto it = nodes.find(n->node_id);
    if (it != nodes.end()) {
        it = nodes.erase(it);
    }
    return it;
}

// remove the all instances of the pattern of nodes/orienations from graph

// remove nodes with covg <= thresh from graph
void Graph::remove_low_covg_nodes(const uint &thresh) {
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


void Graph::add_hits_to_kmergraphs(const vector<LocalPRG *> &prgs) {
    uint num_hits[2];
    for (auto pnode : nodes) {
        // copy kmergraph
        pnode.second->kmer_prg = prgs[pnode.second->prg_id]->kmer_prg;
	    assert(pnode.second->kmer_prg == prgs[pnode.second->prg_id]->kmer_prg);
        num_hits[0] = 0;
        num_hits[1] = 0;

        // add hits
        for (auto read : pnode.second->reads)
        {
            for (auto mh = read->hits[pnode.second->node_id].begin();
                 mh != read->hits[pnode.second->node_id].end(); ++mh)
            {
                //bool added = false;
                // update the covg in the kmer_prg
                pnode.second->kmer_prg.nodes[(*mh)->knode_id]->covg[(*mh)->strand] += 1;
                num_hits[(*mh)->strand] += 1;
            }
        }
        //cout << now() << "Added " << num_hits[1] << " hits in the forward direction and " << num_hits[0]
        //     << " hits in the reverse" << endl;
        pnode.second->kmer_prg.num_reads = pnode.second->covg;
    }
}

bool Graph::operator==(const Graph &y) const {
    // false if have different numbers of nodes
    if (y.nodes.size() != nodes.size()) {
        cout << "different num nodes " << nodes.size() << "!=" << y.nodes.size() << endl;
        return false;
    }

    // false if have different nodes
    for (const auto c: nodes) {
        // if node id doesn't exist 
        auto it = y.nodes.find(c.first);
        if (it == y.nodes.end()) {
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


std::ostream & pangenome::operator<<(std::ostream &out, pangenome::Graph const &m) {
    //cout << "printing pangraph" << endl;
    for (const auto &n : m.nodes) {
        cout << n.second->prg_id << endl;
    }

    return out;
}
