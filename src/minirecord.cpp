#include <cassert>
#include <functional>
#include <iostream>
#include <cstring>
#include "minirecord.h"
#include "path.h"

using namespace std;

MiniRecord::MiniRecord(){};

MiniRecord::MiniRecord(uint32_t p, Path q, bool c): prg_id(p), path(q), strand(c){
    //path.initialize(q.path);
};

MiniRecord::~MiniRecord() {};

bool MiniRecord::operator == (const MiniRecord& y) const {
    //cout << prg_id << "," << y.prg_id << endl;
    //cout << path << "," << y.path << endl;
    if (prg_id != y.prg_id) {return false;}
    if (!(path == y.path)) {return false;}
    if (strand != y.strand) {return false;}
    return true;
}

std::ostream& operator<< (std::ostream& out, MiniRecord const& m) {
    out << "(" << m.prg_id << ", " << m.path << ", " << m.strand << ")";
    return out ;
}

std::istream& operator>> (std::istream& in, MiniRecord& m) {
    in.ignore('(');
    in >> m.prg_id;
    in.ignore(' ');
    in >> m.path;
    in.ignore(' ');
    in >> m.strand;
    in.ignore(')');
    return in ;
}
