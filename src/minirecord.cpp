#include <cassert>
#include <functional>
#include <iostream>
#include <cstring>
#include "minirecord.h"
#include "path.h"

using namespace std;

MiniRecord::MiniRecord(uint32_t p, Path q): prg_id(p), path(q) {
    //path.initialize(q.path);
};

MiniRecord::~MiniRecord() {};

bool MiniRecord::operator == (const MiniRecord& y) const {
    //cout << prg_id << "," << y.prg_id << endl;
    //cout << path << "," << y.path << endl;
    if (prg_id != y.prg_id) {return false;}
    if (!(path == y.path)) {return false;}
    return true;
}

std::ostream& operator<< (std::ostream & out, MiniRecord const& m) {
    out << "(" << m.prg_id << ", " << m.path << ")";
    return out ;
}
