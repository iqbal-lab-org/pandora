#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include "pangenome/pannode.h"
#include "pangenome/pansample.h"
#include "pangenome/panread.h"
#include "minihit.h"
#include "utils.h"
#include "localPRG.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace pangenome;

Node::Node(const uint32_t i, const uint32_t j, const string n) : prg_id(i), node_id(j), name(n), covg(1) {}

/*// copy constructor
Node::Node(const Node& other)
{
    prg_id = other.prg_id;
    node_id = other.node_id;
    name = other.name;
    covg = other.covg;
    kmer_prg = other.kmer_prg;
    reads = other.reads;
    samples = other.samples;
}

// Assignment operator
Node& Node::operator=(const Node& other)
{
    // check for self-assignment
    if (this == &other)
        return *this;

    prg_id = other.prg_id;
    node_id = other.node_id;
    name = other.name;
    covg = other.covg;
    kmer_prg = other.kmer_prg;
    reads = other.reads;
    samples = other.samples;

    return *this;
}*/

void Node::remove_read(ReadPtr r) {
    //removes single copy of read
    auto it = find(reads.begin(), reads.end(), r);
    if (it != reads.end()) {
        covg -= 1;
        reads.erase(it);
        //it = find(reads.begin(), reads.end(), r);
    }
}

string Node::get_name() const {
    if (prg_id != node_id) {
        return name + "." + to_string(node_id);
    } else {
        return name;
    }
}

void Node::add_path(const vector<KmerNodePtr> &kmp) {
    for (uint32_t i = 0; i != kmp.size(); ++i) {
        assert(kmp[i]->id < kmer_prg.nodes.size() and kmer_prg.nodes[kmp[i]->id] != nullptr);
        kmer_prg.nodes[kmp[i]->id]->covg[0] += 1;
        kmer_prg.nodes[kmp[i]->id]->covg[1] += 1;
    }
}

void Node::get_read_overlap_coordinates(vector<vector<uint32_t>> &read_overlap_coordinates) {
    read_overlap_coordinates.reserve(reads.size());
    vector<uint32_t> coordinate;

    auto read_count = 0;
    for (const auto &read_ptr : reads) {
        read_count++;
        if (read_ptr->hits.at(prg_id).size() < 2)
            continue;

        auto hit_ptr_iter = read_ptr->hits.at(prg_id).begin();
        uint32_t start = (*hit_ptr_iter)->read_start_position;
        uint32_t end = 0;
        for (const auto &hit_ptr : read_ptr->hits.at(prg_id)) {
            start = min(start, hit_ptr->read_start_position);
            end = max(end, hit_ptr->read_start_position + hit_ptr->prg_path.length());
        }

        assert(end > start or assert_msg(
                "Error finding the read overlap coordinates for node " << name << " and read " << read_ptr->id
                                                                       << " (the " << read_count << "th on this node)"
                                                                       << endl << "Found end " << end
                                                                       << " after found start " << start));
        coordinate = {read_ptr->id, start, end, (*hit_ptr_iter)->strand};
        read_overlap_coordinates.push_back(coordinate);
    }

    if (not read_overlap_coordinates.empty()) {
        sort(read_overlap_coordinates.begin(), read_overlap_coordinates.end(),
             [](const vector<uint32_t> &a, const vector<uint32_t> &b) {
                 for (uint32_t i = 0; i < a.size(); ++i) {
                     if (a[i] != b[i]) { return a[i] < b[i]; }
                 }
                 return false;
             });
    }

}

void Node::output_samples(const std::shared_ptr<LocalPRG> &prg, const string &outdir, const uint32_t w,
                          const string &vcf_ref) {
    vector<KmerNodePtr> kmp;
    kmp.reserve(800);
    vector<LocalNodePtr> refpath, sample_lmp;
    refpath.reserve(100);
    sample_lmp.reserve(100);

    // if we passed in a string to be the reference, find the corresponding path
    if (!vcf_ref.empty()) {
        refpath = prg->prg.nodes_along_string(vcf_ref);
        if (refpath.empty()) {
            refpath = prg->prg.nodes_along_string(rev_complement(vcf_ref));
        }
        if (refpath.empty()) {
            cout << now() << "Could not find reference sequence for " << name
                 << "in the PRG so using the closest path" << endl;
            kmer_prg.set_p(0.01);
            kmer_prg.num_reads = covg;
            kmer_prg.find_max_path(kmp);
            refpath = prg->localnode_path_from_kmernode_path(kmp, w);
        }
    } else {
        kmer_prg.set_p(0.01);
        kmer_prg.num_reads = covg;
        kmer_prg.find_max_path(kmp);
        refpath = prg->localnode_path_from_kmernode_path(kmp, w);
    }

    // create a with respect to this ref
    VCF vcf;
    prg->build_vcf(vcf, refpath);
    vcf.save(outdir + "/" + name + ".multisample.vcf", true, true, true, true, true, true, true);
    uint32_t count = 0;
    for (const auto &s : samples) {
        //cout << "new sample" << endl;
        count = 0;
        for (const auto &p : s->paths[prg_id]) {
            /*cout << s->name << " ";
            for (uint i=0; i!=p.size(); ++i)
            {
                cout << p[i]->id << " ";
            }
            cout << endl;*/
            sample_lmp = prg->localnode_path_from_kmernode_path(p, w);
            /*cout << "sample lmp:" << endl;
            for (uint i=0; i!=sample_lmp.size(); ++i)
            {
                cout << sample_lmp[i]->id << "->";
            }
            cout << endl;*/
            if (count == 0) {
                prg->add_sample_gt_to_vcf(vcf, refpath, sample_lmp, s->name);
                //prg->add_sample_covgs_to_vcf(vcf, kmer_prg, lmp, p, s->name);

            } else {
                prg->add_sample_gt_to_vcf(vcf, refpath, sample_lmp, s->name + to_string(count));
                //prg->add_sample_covgs_to_vcf(vcf, kmer_prg, lmp, p, s->name + to_string(count));
            }
            sample_lmp.clear();
            count++;
            //cout << "finished adding sample " << s->name << " path " << count << endl;
        }
    }
    vcf.save(outdir + "/" + name + ".multisample.vcf", true, true, true, true, true, true, true);
}

bool Node::operator==(const Node &y) const {
    return (node_id == y.node_id);
}

bool Node::operator!=(const Node &y) const {
    return (node_id != y.node_id);
}

bool Node::operator<(const Node &y) const {
    return (node_id < y.node_id);
}

std::ostream &pangenome::operator<<(std::ostream &out, pangenome::Node const &n) {
    out << n.node_id << "," << n.prg_id << " covg: " << n.covg;
    return out;
}
