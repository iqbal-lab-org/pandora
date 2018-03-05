#include <cassert>
#include <functional>
#include <iostream>
#include <cstring>
#include "vcfrecord.h"

using namespace std;

VCFRecord::VCFRecord(std::string c, uint32_t p, std::string r, std::string a, std::string i, std::string g): chrom(c), pos(p), id("."), ref(r), alt(a), qual("."), filter("."), info(i), format("GT:COV") {
    // fix so no empty strings
    if (ref == "") {ref = ".";}
    if (alt == "") {alt = ".";}

    // classify variants as SNPs, INDELs PH_SNPs or COMPLEX
    // need to think about how to handle cases where there are more than 2 options, not all of one type
    if (info == ".")
    {
        if (ref == "." and alt == ".") {}
        else if (ref == "." or alt == ".") {info = "SVTYPE=INDEL";}
        else if (ref.length() == 1 and alt.length() == 1) {info = "SVTYPE=SNP";}
	else if (alt.length() == ref.length()) {info = "SVTYPE=PH_SNPs";}
        else if (ref.length() < alt.length() and ref.compare(0, ref.length(), alt, 0, ref.length()) == 0) {info = "SVTYPE=INDEL";}
        else if (alt.length() < ref.length() and alt.compare(0, alt.length(), ref, 0, alt.length()) == 0) {info = "SVTYPE=INDEL";}
        else {info = "SVTYPE=COMPLEX";}
    }

    // add graph type info
    if (g != "")
    {
        info += ";";
        info += g;
    }
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
    out << m.chrom << "\t" << m.pos << "\t" << m.id << "\t" << m.ref << "\t" << m.alt << "\t" << m.qual << "\t" << m.filter << "\t" << m.info << "\t" << m.format;
    for (uint i=0; i!=m.samples.size(); ++i)
    {
	out << "\t" << m.samples[i];
    }
    out << endl;
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
    in.ignore(1,'\t');
    in >> m.format;
    string token;
    while (in >> token)
    {
	m.samples.push_back(token);
    }
    return in;
}


