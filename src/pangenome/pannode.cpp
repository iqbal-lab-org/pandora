#include <iostream>
#include <fstream>
#include <cassert>
#include "pangenome/pannode.h"
#include "pangenome/pansample.h"
#include "utils.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace pangenome;

Node::Node (const uint32_t i, const uint32_t j, const string n): prg_id(i), node_id(j), name(n), covg(1) {}

string Node::get_name()
{
    if (prg_id != node_id)
    {
	    return name + "." + to_string(node_id);
    } else {
	    return name;
    }
}

void Node::add_path(const vector<KmerNodePtr>& kmp)
{
    for (uint i=0; i!=kmp.size(); ++i)
    {
	    assert(kmer_prg.nodes.find(kmp[i]->id)!=kmer_prg.nodes.end() ||
                       assert_msg("Must have wrong kmergraph as has different nodes"));
        kmer_prg.nodes[kmp[i]->id]->covg[0] += 1;
	    kmer_prg.nodes[kmp[i]->id]->covg[1] += 1;
    }
}

void Node::output_samples(const LocalPRG* prg, const string& prefix, const uint w)
{
    vector<KmerNodePtr> kmp;
    kmp.reserve(800);
    vector<LocalNodePtr> lmp, sample_lmp;
    lmp.reserve(100);
    sample_lmp.reserve(100);

    // find best ref
    kmer_prg.set_p(0.01);
    kmer_prg.num_reads = covg;
    //cout << "have " << covg << " samples with this node" << endl;
    kmer_prg.find_max_path(kmp);
    lmp = prg->localnode_path_from_kmernode_path(kmp, w);
    /*cout << "best ref: " << endl;
    for (uint i=0; i!=lmp.size(); ++i)
    {
	cout << lmp[i]->id << "->";
    }
    cout << endl;*/

    // create a with respect to this ref
    VCF vcf;
    prg->build_vcf(vcf, lmp);
    vcf.save(prefix + "." + name + ".multisample.vcf", true, true, true, true, true, true, true);
    uint count = 0;
    for (auto s : samples)
    {
	//cout << "new sample" << endl;
	count = 0;
	for (const auto &p : s->paths[prg_id])
        {
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
	    if (count == 0)
	    {
                prg->add_sample_to_vcf(vcf, lmp, sample_lmp, s->name); 
	    } else {
		prg->add_sample_to_vcf(vcf, lmp, sample_lmp, s->name + to_string(count));
	    }
	    sample_lmp.clear();
	    count++;
	    //cout << "finished adding sample " << s->name << " path " << count << endl;
	}
    }
    vcf.save(prefix + "." + name + ".multisample.vcf", true, true, true, true, true, true, true);
    vcf.write_aligned_fasta(prefix + "." + name + ".multisample.fa", lmp);
}

/*// copy constructor
Node::Node(const Node& other)
{
    prg_id = other.prg_id;
    //node_id = other.node_id;
    name = other.name;
    covg = other.covg;
    kmer_prg = other.kmer_prg;
    //edges = other.edges; // shallow copies, so will point to same edges and reads
    //reads = other.reads;
}

// Assignment operatorNode& KmerNode::operator=(const KmerNode& other)
Node& Node::operator=(const Node& other)
{
    // check for self-assignment
    if (this == &other)
        return *this;

    prg_id = other.prg_id;
    //node_id = other.node_id;
    name = other.name;
    covg = other.covg;
    kmer_prg = other.kmer_prg;
    //edges = other.edges; // shallow copies, so will point to same edges and reads
    //reads = other.reads;

    return *this;
}*/

bool Node::operator == (const Node& y) const {
    return (node_id == y.node_id);
}

bool Node::operator != (const Node& y) const {
    return (node_id != y.node_id);
}

bool Node::operator < (const Node& y) const {
    return (node_id < y.node_id);
}

std::ostream& pangenome::operator<< (std::ostream & out, pangenome::Node const& n) {
    out << n.prg_id << " covg: " << n.covg;
    return out ;
}
