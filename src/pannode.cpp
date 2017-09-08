#include <iostream>
#include <string>
#include <fstream>
#include "pannode.h"
#include "utils.h"

using namespace std;

PanNode::PanNode (const uint32_t i, const uint32_t j, const string n): prg_id(i), node_id(j), name(n), covg(1) {}

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
