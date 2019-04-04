#include <iostream>
#include <cassert>

#include <boost/log/trivial.hpp>

#include "kmernode.h"
#include "utils.h"


KmerNode::KmerNode(uint32_t i, const prg::Path &p) : id(i), path(p),
                                                khash(std::numeric_limits<uint64_t>::max()), num_AT(0) {
    this->covg_new = {{0, 0}};
}

// copy constructor
KmerNode::KmerNode(const KmerNode &other) {
    id = other.id;
    path = other.path;
    khash = other.khash;
    num_AT = other.num_AT;
    // NB we don't do edges

    this->covg_new = other.covg_new;
}

// Assignment operator
KmerNode &KmerNode::operator=(const KmerNode &other) {
    // check for self-assignment
    if (this == &other)
        return *this;

    id = other.id;
    path = other.path;
    khash = other.khash;
    num_AT = other.num_AT;
    // NB we don't do edges
    this->covg_new = other.covg_new;

    return *this;
}

void KmerNode::increment_covg(const bool &strand, const uint32_t &sample_id) {
    assert(this->covg_new.size() > sample_id);

    if (strand)
        this->covg_new[sample_id].first++;
    else
        this->covg_new[sample_id].second++;
}

uint32_t KmerNode::get_covg(const bool &strand, const uint32_t &sample_id) {
    if (this->covg_new.size() <= sample_id)
        return 0;

    if (strand)
        return this->covg_new[sample_id].first;
    else
        return this->covg_new[sample_id].second;
}

void KmerNode::set_covg(const uint32_t &value, const bool &strand, const uint32_t &sample_id) {
    assert(this->covg_new.size() > sample_id);
    if (strand)
        this->covg_new[sample_id].first = value;
    else
        this->covg_new[sample_id].second = value;
}

std::ostream &operator<<(std::ostream &out, const KmerNode &kmer_node) {
    out << kmer_node.id << " " << kmer_node.path << " ";
    for (const auto &sample_coverage: kmer_node.covg_new) {
        out << "(" << sample_coverage.first << ", " << sample_coverage.second << ") ";
    }
    return out;
}

bool KmerNode::operator==(const KmerNode &y) const {
    return path == y.path;
}
