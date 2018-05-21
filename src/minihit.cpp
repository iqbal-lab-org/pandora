#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <algorithm>
#include "minirecord.h"
#include "minihit.h"
#include "path.h"

using namespace std;

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

MinimizerHit::MinimizerHit(const uint32_t i, const Minimizer& m, const MiniRecord* r): read_id(i), read_start_position(m.pos.start), prg_id(r->prg_id), prg_path(r->path), knode_id(r->knode_id), strand((m.strand == r->strand))
{
    //cout << *m << " + " << *r << " = " << "(" << read_id << ", " << read_start_position << ", " << prg_id << ", " << prg_path << ", " << strand << ")" << endl;
    assert(m.pos.length==prg_path.length());
    assert(read_id < std::numeric_limits<uint32_t>::max() || assert_msg("Variable sizes too small to handle this number of reads"));
    assert(prg_id < std::numeric_limits<uint32_t>::max() || assert_msg("Variable sizes too small to handle this number of prgs"));
};

MinimizerHit::MinimizerHit(const uint32_t i, const Interval j, const uint32_t k, const Path p, const uint32_t n, const bool c): read_id(i), read_start_position(j.start), prg_id(k), knode_id(n), strand(c)
{
    prg_path.initialize(p.path);
    assert(j.length==prg_path.length());
};

bool MinimizerHit::operator == (const MinimizerHit& y) const {
    if (read_id != y.read_id) {return false;}
    if (!(read_start_position == y.read_start_position)) {return false;}
    if (prg_id != y.prg_id) {return false;}
    if (!(prg_path == y.prg_path)) {return false;}
    if (strand != y.strand) {return false;}
    return true;
}

bool MinimizerHit::operator < ( const MinimizerHit& y) const
{
    // first by the read they map too - should all be the same
    if (read_id < y.read_id){ return true; }
    if (y.read_id < read_id) { return false; }

    // then by the prg they map too
    if (prg_id < y.prg_id){ return true; }
    if (y.prg_id < prg_id) { return false; }

    // then by direction NB this bias is in favour of the forward direction
    if (strand < y.strand){ return false; }
    if (y.strand < strand) { return true; }

    // then by position on query string
    if (read_start_position < y.read_start_position) { return true; }
    if (y.read_start_position < read_start_position) { return false; }

    // then by position on target string
    if (prg_path < y.prg_path) { return true; }
    if (y.prg_path < prg_path) { return false; }

    return false;
}

std::ostream& operator<< (std::ostream & out, MinimizerHit const& m) {
    out << "(" << m.read_id << ", " << m.read_start_position << ", " << m.prg_id << ", " << m.prg_path << ", " << m.strand << ", " << m.knode_id << ")";
    return out ;
}

