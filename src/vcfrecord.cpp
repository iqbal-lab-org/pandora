#include <cassert>
#include <functional>
#include <iostream>
#include <cstring>
#include "vcfrecord.h"

using namespace std;

VCFRecord::VCFRecord(std::string c, uint32_t p, std::string r, std::string a): chrom(c), pos(p), id("."), ref(r), alt(a), qual("."), filter("."), info("."){
};

VCFRecord::VCFRecord(): chrom("."), pos(0), id("."), ref("."), alt("."), qual("."), filter("."), info("."){
};

VCFRecord::~VCFRecord() {};

bool VCFRecord::operator == (const VCFRecord& y) const {
    if (chrom != y.chrom){return false;}
    if (pos != y.pos){return false;}
    if (ref != y.ref){return false;}
    if (alt != y.alt){return false;}
    return true;
}

bool VCFRecord::operator <  (const VCFRecord& y) const {
    if (chrom < y.chrom){return true;}
    if (chrom > y.chrom){return false;}
    if (pos < y.pos){return true;}
    if (pos > y.pos){return false;}
    if (ref < y.ref){return true;}
    if (ref > y.ref){return false;}
    if (alt < y.alt){return true;}
    if (alt > y.alt){return false;}
    return false;
}


std::ostream& operator<< (std::ostream& out, VCFRecord const& m) {
    out << m.chrom << "\t" << m.pos << "\t" << m.id << "\t" << m.ref << "\t" << m.alt << "\t" << m.qual << "\t" << m.filter << "\t" << m.info << endl;
    return out;
}

std::istream& operator>> (std::istream& in, VCFRecord& m) {
    in >> m.chrom;
    in.ignore(1,'\t');
    in >> m.pos;
    in.ignore(1,'\t');
    in >> m.id;
    in.ignore(1,'\t');
    in >> m.ref;
    in.ignore(1,'\t');
    in >> m.alt;
    in.ignore(1,'\t');
    in >> m.qual;
    in.ignore(1,'\t');
    in >> m.filter;
    in.ignore(1,'\t');
    in >> m.info;
    in.ignore(1,'\n');
    return in;
}


