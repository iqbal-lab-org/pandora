#include <functional>
#include <iostream>
#include "minirecord.h"

MiniRecord::MiniRecord() {};

MiniRecord::MiniRecord(
    const uint32_t p, const prg::Path &q, const uint32_t n, const bool c)
    : prg_id(p)
    , knode_id(n)
    , strand(c) {};

MiniRecord::~MiniRecord() {};

bool MiniRecord::operator==(const MiniRecord& y) const
{
    if (prg_id != y.prg_id) {
        return false;
    }
    if (knode_id != y.knode_id) {
        return false;
    }
    if (strand != y.strand) {
        return false;
    }
    return true;
}

std::ostream& operator<<(std::ostream& out, MiniRecord const& m)
{
    out << "(" << m.prg_id << ", " << m.knode_id << ", " << m.strand
        << ")";
    return out;
}

std::istream& operator>>(std::istream& in, MiniRecord& m)
{
    in.ignore(1, '(');
    in >> m.prg_id;
    in.ignore(2, ' ');
    in >> m.knode_id;
    in.ignore(2, ' ');
    in >> m.strand;
    in.ignore(1, ')');
    return in;
}
