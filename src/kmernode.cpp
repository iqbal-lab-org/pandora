#include <cassert>
#include <iostream>

#include <boost/log/trivial.hpp>

#include "kmernode.h"
#include "utils.h"

KmerNode::KmerNode(uint32_t i, const prg::Path& p)
    : id(i)
    , path(p)
    , khash(std::numeric_limits<uint64_t>::max())
    , num_AT(0)
{
}

// copy constructor
KmerNode::KmerNode(const KmerNode& other)
{
    id = other.id;
    path = other.path;
    khash = other.khash;
    num_AT = other.num_AT;
    // NB we don't do edges
}

// Assignment operator
KmerNode& KmerNode::operator=(const KmerNode& other)
{
    // check for self-assignment
    if (this == &other)
        return *this;

    id = other.id;
    path = other.path;
    khash = other.khash;
    num_AT = other.num_AT;
    // NB we don't do edges

    return *this;
}

std::ostream& operator<<(std::ostream& out, const KmerNode& kmer_node)
{
    out << kmer_node.id << " " << kmer_node.path << " ";
    return out;
}

bool KmerNode::operator==(const KmerNode& y) const { return path == y.path; }
