#include <cassert>
#include <functional>
#include <iostream>
#include <cstring>
#include "minirecord.h"
#include "path.h"

using namespace std;

MiniRecord::MiniRecord(uint32_t p, Path q): prg_id(p), path(q) {};

MiniRecord::~MiniRecord() {};

bool MiniRecord::operator == (const MiniRecord& y) const {
    if (prg_id != y.prg_id) {return false;}
    if (!(path == y.path)) {return false;}
    return true;
}

std::ostream& operator<< (std::ostream & out, MiniRecord const& m) {
    out << "(" << m.prg_id << ", " << m.path << ")";
    return out ;
}
