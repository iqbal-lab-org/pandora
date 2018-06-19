#include <cassert>
#include <cmath>
#include <iostream>
#include "vcfrecord.h"
#include "utils.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

VCFRecord::VCFRecord(std::string c, uint32_t p, std::string r, std::string a, std::string i, std::string g) : chrom(c),
                                                                                                              pos(p),
                                                                                                              id("."),
                                                                                                              ref(r),
                                                                                                              alt(a),
                                                                                                              qual("."),
                                                                                                              filter("."),
                                                                                                              info(i),
                                                                                                              format({"GT"}) {
    // fix so no empty strings
    if (ref == "") { ref = "."; }
    if (alt == "") { alt = "."; }

    // classify variants as SNPs, INDELs PH_SNPs or COMPLEX
    // need to think about how to handle cases where there are more than 2 options, not all of one type
    if (info == ".") {
        if (ref == "." and alt == ".") {}
        else if (ref == "." or alt == ".") { info = "SVTYPE=INDEL"; }
        else if (ref.length() == 1 and alt.length() == 1) { info = "SVTYPE=SNP"; }
        else if (alt.length() == ref.length()) { info = "SVTYPE=PH_SNPs"; }
        else if (ref.length() < alt.length() and
                 ref.compare(0, ref.length(), alt, 0, ref.length()) == 0) { info = "SVTYPE=INDEL"; }
        else if (alt.length() < ref.length() and
                 alt.compare(0, alt.length(), ref, 0, alt.length()) == 0) { info = "SVTYPE=INDEL"; }
        else { info = "SVTYPE=COMPLEX"; }
    }

    // add graph type info
    if (g != "") {
        info += ";";
        info += g;
    }
};

VCFRecord::VCFRecord() : chrom("."), pos(0), id("."), ref("."), alt("."), qual("."), filter("."), info(".") {
};

VCFRecord::~VCFRecord() {};

bool VCFRecord::operator==(const VCFRecord &y) const {
    if (chrom != y.chrom) { return false; }
    if (pos != y.pos) { return false; }
    if (ref != y.ref) { return false; }
    if (alt != y.alt) { return false; }
    return true;
}

void VCFRecord::add_formats(std::vector<std::string> formats) {
    for (auto s : formats){
        if (find(format.begin(), format.end(),s) == format.end())
            format.push_back(s);
    }
}

float logfactorial(uint32_t n){
    float ret = 0;
    for (uint32_t i=1; i<=n; ++i){
        ret += log(i);
    }
    return ret;
}

void VCFRecord::likelihood(const uint32_t& expected_depth_covg, const float& error_rate) {
    for (auto sample : samples) {
        if (sample.find("REF_MEAN_FWD_COVG") != sample.end() and sample.find("REF_MEAN_REV_COVG") != sample.end()
            and sample.find("REF_MEAN_FWD_COVG") != sample.end() and sample.find("REF_MEAN_REV_COVG") != sample.end()) {
            auto c1 = sample["REF_MEAN_FWD_COVG"] + sample["REF_MEAN_REV_COVG"];
            auto c2 = sample["ALT_MEAN_FWD_COVG"] + sample["ALT_MEAN_REV_COVG"];
            if (c1 > 0)
                sample["REF_LIKELIHOOD"] = c1 * log(expected_depth_covg) - expected_depth_covg
                                       - logfactorial(c1) + c2 * log(error_rate);
            else
                sample["REF_LIKELIHOOD"] = numeric_limits<float>::lowest();
            if (c2 > 0)
                sample["ALT_LIKELIHOOD"] = c2 * log(expected_depth_covg) - expected_depth_covg
                                       - logfactorial(c2) + c1 * log(error_rate);
            else
                sample["REF_LIKELIHOOD"] = numeric_limits<float>::lowest();
        }
    }
    add_formats({"REF_LIKELIHOOD","ALT_LIKELIHOOD"});
}

void VCFRecord::confidence(){
    for (auto sample : samples) {
        if (sample.find("REF_LIKELIHOOD") != sample.end() and sample.find("ALT_LIKELIHOOD") != sample.end()){
            sample["CONFIDENCE"] = abs(sample["REF_LIKELIHOOD"]-sample["ALT_LIKELIHOOD"]);
        }
    }
    add_formats({"CONFIDENCE"});
}

void VCFRecord::swap_ref_and_alt_properties(unordered_map<std::string, uint8_t>& sample){
    swap(sample["REF_MEAN_FWD_COVG"],sample["ALT_MEAN_FWD_COVG"]);
    swap(sample["REF_MEAN_REV_COVG"],sample["ALT_MEAN_REV_COVG"]);
    swap(sample["REF_MED_FWD_COVG"],sample["ALT_MED_FWD_COVG"]);
    swap(sample["REF_MED_REV_COVG"],sample["ALT_MED_REV_COVG"]);
    swap(sample["REF_SUM_FWD_COVG"],sample["ALT_SUM_FWD_COVG"]);
    swap(sample["REF_SUM_REV_COVG"],sample["ALT_SUM_REV_COVG"]);
    swap(sample["REF_LIKELIHOOD"],sample["ALT_LIKELIHOOD"]);
}

void VCFRecord::regenotype(const uint8_t confidence_threshold){
    for (auto sample : samples) {
        if (sample.find("CONFIDENCE") != sample.end()){
            if (sample["CONFIDENCE"] > confidence_threshold){
                if (sample["GT"] == 0 and sample["ALT_LIKELIHOOD"] > sample["REF_LIKELIHOOD"]){
                    sample["GT"] = 1;
                    swap_ref_and_alt_properties(sample);
                } else if (sample["GT"] == 1 and sample["ALT_LIKELIHOOD"] < sample["REF_LIKELIHOOD"]){
                    sample["GT"] = 0;
                    swap_ref_and_alt_properties(sample);
                }
            } else {
                sample.erase("GT");
            }
        }
    }
}

bool VCFRecord::operator<(const VCFRecord &y) const {
    if (chrom < y.chrom) { return true; }
    if (chrom > y.chrom) { return false; }
    if (pos < y.pos) { return true; }
    if (pos > y.pos) { return false; }
    if (ref < y.ref) { return true; }
    if (ref > y.ref) { return false; }
    if (alt < y.alt) { return true; }
    if (alt > y.alt) { return false; }
    return false;
}


std::ostream &operator<<(std::ostream &out, VCFRecord const &m) {
    out << m.chrom << "\t" << m.pos << "\t" << m.id << "\t" << m.ref << "\t" << m.alt << "\t" << m.qual << "\t"
        << m.filter << "\t" << m.info << "\t";

    string last_format;
    if (!m.format.empty())
        last_format = m.format[m.format.size()-1];

    for (auto s : m.format){
        out << s;
        if (s != last_format)
            out << ":";
    }

    for(auto s : m.samples){
        out << "\t";
        for (auto f : m.format){
            if (s.find(f)!=s.end())
                out << +s[f];
            else
                out << ".";

            if (f != last_format)
                out << ":";
        }
    }

    out << endl;
    return out;
}

std::istream &operator>>(std::istream &in, VCFRecord &m) {
    string token;
    vector<string> sample_strings;
    unordered_map<string, uint8_t> sample_data;
    in >> m.chrom;
    in.ignore(1, '\t');
    in >> m.pos;
    in.ignore(1, '\t');
    in >> m.id;
    in.ignore(1, '\t');
    in >> m.ref;
    in.ignore(1, '\t');
    in >> m.alt;
    in.ignore(1, '\t');
    in >> m.qual;
    in.ignore(1, '\t');
    in >> m.filter;
    in.ignore(1, '\t');
    in >> m.info;
    in.ignore(1, '\t');
    in >> token;
    m.format = split(token, ":");
    sample_data.reserve(m.format.size());
    while (in >> token) {
        sample_strings = split(token, ":");
        assert(sample_strings.size()==m.format.size() or assert_msg("sample data does not fit format"));
        m.samples.push_back(sample_data);
        for (uint32_t i=0; i<m.format.size(); ++i){
            if (sample_strings[i] != ".")
                m.samples.back()[m.format[i]] = stoi(sample_strings[i]);
        }
    }
    return in;
}


