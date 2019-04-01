#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <algorithm>
#include "minirecord.h"
#include "minihit.h"
#include "prg/path.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


MinimizerHit::MinimizerHit(const uint32_t i, const Minimizer &m, const MiniRecord *r)
        : read_id(i), read_start_position(m.pos.start), prg_id(r->prg_id), prg_path(r->path), kmer_node_id(r->knode_id),
          is_forward((m.strand == r->strand)) {
    assert(m.pos.length == prg_path.length());
    assert(read_id < std::numeric_limits<uint32_t>::max() ||
           assert_msg("Variable sizes too small to handle this number of reads"));
    assert(prg_id < std::numeric_limits<uint32_t>::max() ||
           assert_msg("Variable sizes too small to handle this number of prgs"));
};


MinimizerHit::MinimizerHit(const uint32_t read_id, const Interval read_interval, const uint32_t prg_id,
                           const prg::Path prg_path, const uint32_t kmer_node_id, const bool is_forward)
        : read_id(read_id), read_start_position(read_interval.start), prg_id(prg_id), kmer_node_id(kmer_node_id),
          is_forward(is_forward) {
    this->prg_path.initialize(prg_path.path);
    assert(read_interval.length == this->prg_path.length());
};


bool MinimizerHit::operator==(const MinimizerHit &y) const {
    if (read_id != y.read_id) { return false; }
    if (!(read_start_position == y.read_start_position)) { return false; }
    if (prg_id != y.prg_id) { return false; }
    if (!(prg_path == y.prg_path)) { return false; }
    if (is_forward != y.is_forward) { return false; }
    return true;
}


bool MinimizerHit::operator<(const MinimizerHit &y) const {
    // first by the read they map too - should all be the same
    if (read_id < y.read_id) { return true; }
    if (y.read_id < read_id) { return false; }

    // then by the prg they map too
    if (prg_id < y.prg_id) { return true; }
    if (y.prg_id < prg_id) { return false; }

    // then by direction NB this bias is in favour of the forward direction
    if (is_forward < y.is_forward) { return false; }
    if (y.is_forward < is_forward) { return true; }

    // then by position on query string
    if (read_start_position < y.read_start_position) { return true; }
    if (y.read_start_position < read_start_position) { return false; }

    // then by position on target string
    if (prg_path < y.prg_path) { return true; }
    if (y.prg_path < prg_path) { return false; }

    return false;
}


std::ostream &operator<<(std::ostream &out, MinimizerHit const &m) {
    out << "(" << m.read_id << ", " << m.read_start_position << ", " << m.prg_id << ", " << m.prg_path << ", "
        << m.is_forward << ", " << m.kmer_node_id << ")";
    return out;
}

