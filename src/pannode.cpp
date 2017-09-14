#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include "pannode.h"
#include "utils.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

PanNode::PanNode (const uint32_t i, const uint32_t j, const string n): prg_id(i), node_id(j), name(n), covg(1) {}

string PanNode::get_name()
{
    if (prg_id != node_id)
    {
	return name + "." + to_string(node_id);
    } else {
	return name;
    }
}

void PanNode::add_path(const vector<KmerNode*>& kmp)
{
    for (uint i=0; i!=kmp.size(); ++i)
    {
	assert(kmer_prg.nodes.find(kmp[i]->id)!=kmer_prg.nodes.end() || assert_msg("Must have wrong kmergraph as has different nodes"));
        kmer_prg.nodes[kmp[i]->id]->covg[0] += 1;
	kmer_prg.nodes[kmp[i]->id]->covg[1] += 1;
    }
}

/*void PanNode::output_samples_vcf(const LocalPRG* prg, const uint w)
{
    vector<KmerNode*> kmp;
    kmp.reserve(800);
    vector<LocalNode*> lmp, sample_lmp;
    lmp.reserve(100);
    sample_lmp.reserve(100);

    // find best ref
    kmer_prg.find_max_path(kmp);
    lmp = prg->localnode_path_from_kmernode_path(kmp, w);

    // create a with respect to this ref
    prg->build_vcf();
    for each sample
    {
	sample_lmp.clear();
	sample_lmp = prg->localnode_path_from_kmernode_path(sample_kmp, w);
        add_sample_to_vcf(sample_lmp); 
	sample_lmp.clear();
    }
    prg->vcf.save(prefix + "." + new_name + ".kmlp.vcf", true, true, true, true, true, true, true);
}*/
	
/*// copy constructor
PanNode::PanNode(const PanNode& other)
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
PanNode& PanNode::operator=(const PanNode& other)
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

bool PanNode::operator == (const PanNode& y) const {
    return (node_id == y.node_id);
}

bool PanNode::operator != (const PanNode& y) const {
    return (node_id != y.node_id);
}

bool PanNode::operator < (const PanNode& y) const {
    return (node_id < y.node_id);
}

std::ostream& operator<< (std::ostream & out, PanNode const& n) {
    out << n.prg_id << " covg: " << n.covg;
    return out ;
}
