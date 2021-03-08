#include <functional>
#include <iostream>
#include "minirecord.h"
#include "minihit.h"
#include "prg/path.h"

MinimizerHit::MinimizerHit(const uint32_t i, const Minimizer& minimizer_from_read,
    const MiniRecord& minimizer_from_PRG)
    : read_id { i }
    , read_start_position { minimizer_from_read.pos_of_kmer_in_read.start }
    , read_strand { minimizer_from_read.is_forward_strand }
    , minimizer_from_PRG { minimizer_from_PRG }
{
    const bool both_minimizers_have_same_length
        = minimizer_from_read.pos_of_kmer_in_read.length
        == minimizer_from_PRG.path.length();
    if (!both_minimizers_have_same_length) {
        fatal_error("Error when storing minimizers: minimizer from read/sequence "
                    "and from PRG have different lengths");
    }
}

bool MinimizerHit::operator==(const MinimizerHit& y) const
{
    if (get_read_id() != y.get_read_id()) {
        return false;
    }
    if (!(get_read_start_position() == y.get_read_start_position())) {
        return false;
    }
    if (get_prg_id() != y.get_prg_id()) {
        return false;
    }
    if (!(get_prg_path() == y.get_prg_path())) {
        return false;
    }
    if (is_forward() != y.is_forward()) {
        return false;
    }
    return true;
}

bool MinimizerHit::operator<(const MinimizerHit& y) const
{
    // first by the read they map too - should all be the same
    if (get_read_id() < y.get_read_id()) {
        return true;
    }
    if (y.get_read_id() < get_read_id()) {
        return false;
    }

    // then by the prg they map too
    if (get_prg_id() < y.get_prg_id()) {
        return true;
    }
    if (y.get_prg_id() < get_prg_id()) {
        return false;
    }

    // then by direction NB this bias is in favour of the forward direction
    if (is_forward() < y.is_forward()) {
        return false;
    }
    if (y.is_forward() < is_forward()) {
        return true;
    }

    // then by position on query string
    if (get_read_start_position() < y.get_read_start_position()) {
        return true;
    }
    if (y.get_read_start_position() < get_read_start_position()) {
        return false;
    }

    // then by position on target string
    if (get_prg_path() < y.get_prg_path()) {
        return true;
    }
    if (y.get_prg_path() < get_prg_path()) {
        return false;
    }

    return false;
}

std::ostream& operator<<(std::ostream& out, MinimizerHit const& m)
{
    out << "(" << m.get_read_id() << ", " << m.get_read_start_position() << ", "
        << m.get_prg_id() << ", " << m.get_prg_path() << ", " << m.is_forward() << ", "
        << m.get_kmer_node_id() << ")";
    return out;
}
