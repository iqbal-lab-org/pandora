#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <algorithm>
#include "minirecord.h"
#include "minihit.h"
#include "prg/path.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

MinimizerHit::MinimizerHit(const uint32_t i, const Minimizer &minimizerFromRead, const MiniRecord &minimizerFromPRG) :
        read_id{i}, read_start_position{minimizerFromRead.pos.start}, read_strand{minimizerFromRead.strand}, minimizerFromPRG{minimizerFromPRG}{

    assert(minimizerFromRead.pos.length == minimizerFromPRG.path.length());
    assert(read_id < std::numeric_limits<uint32_t>::max() ||
           assert_msg("Variable sizes too small to handle this number of reads"));
    assert(minimizerFromPRG.prg_id < std::numeric_limits<uint32_t>::max() ||
           assert_msg("Variable sizes too small to handle this number of prgs"));
    assert(minimizerFromRead.pos.length == minimizerFromPRG.path.length());
}


bool MinimizerHit::operator==(const MinimizerHit &y) const {
    if (get_read_id() != y.get_read_id()) { return false; }
    if (!(get_read_start_position() == y.get_read_start_position())) { return false; }
    if (get_prg_id() != y.get_prg_id()) { return false; }
    if (!(get_prg_path() == y.get_prg_path())) { return false; }
    if (is_forward() != y.is_forward()) { return false; }
    return true;
}


bool MinimizerHit::operator<(const MinimizerHit &y) const {
    // first by the read they map too - should all be the same
    if (get_read_id() < y.get_read_id()) { return true; }
    if (y.get_read_id() < get_read_id()) { return false; }

    // then by the prg they map too
    if (get_prg_id() < y.get_prg_id()) { return true; }
    if (y.get_prg_id() < get_prg_id()) { return false; }

    // then by direction NB this bias is in favour of the forward direction
    if (is_forward() < y.is_forward()) { return false; }
    if (y.is_forward() < is_forward()) { return true; }

    // then by position on query string
    if (get_read_start_position() < y.get_read_start_position()) { return true; }
    if (y.get_read_start_position() < get_read_start_position()) { return false; }

    // then by position on target string
    if (get_prg_path() < y.get_prg_path()) { return true; }
    if (y.get_prg_path() < get_prg_path()) { return false; }

    return false;
}


std::ostream &operator<<(std::ostream &out, MinimizerHit const &m) {
    out << "(" << m.get_read_id() << ", " << m.get_read_start_position() << ", " << m.get_prg_id() << ", " << m.get_prg_path() << ", "
        << m.is_forward() << ", " << m.get_kmer_node_id() << ")";
    return out;
}

