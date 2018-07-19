#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
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

void VCFRecord::add_formats(const vector<string>& formats) {
    for (const auto s : formats){
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

    unordered_map<string, float> m;
    m.reserve(3);
    if (regt_samples.size() == 0){
        for (auto sample : samples) {
            regt_samples.push_back(m);
        }
    }

    for (uint_least16_t i=0; i<samples.size(); ++i) {
        if (samples[i].find("REF_MEAN_FWD_COVG") != samples[i].end() and samples[i].find("REF_MEAN_REV_COVG") != samples[i].end()
            and samples[i].find("ALT_MEAN_FWD_COVG") != samples[i].end() and samples[i].find("ALT_MEAN_REV_COVG") != samples[i].end()) {
            auto c1 = samples[i]["REF_MEAN_FWD_COVG"] + samples[i]["REF_MEAN_REV_COVG"];
            auto c2 = samples[i]["ALT_MEAN_FWD_COVG"] + samples[i]["ALT_MEAN_REV_COVG"];
            if (c1 > 0)
                regt_samples[i]["REF_LIKELIHOOD"] = c1 * log(expected_depth_covg) - expected_depth_covg
                                       - logfactorial(c1) + c2 * log(error_rate);
            else
                regt_samples[i]["REF_LIKELIHOOD"] = numeric_limits<float>::lowest();
            if (c2 > 0)
                regt_samples[i]["ALT_LIKELIHOOD"] = c2 * log(expected_depth_covg) - expected_depth_covg
                                       - logfactorial(c2) + c1 * log(error_rate);
            else
                regt_samples[i]["ALT_LIKELIHOOD"] = numeric_limits<float>::lowest();
            regt_samples[i]["DP"] = expected_depth_covg;
        }
    }

    assert(regt_samples.size()==samples.size() or assert_msg(regt_samples.size()<< "!=" << samples.size()));
    add_formats({"REF_LIKELIHOOD","ALT_LIKELIHOOD", "DP"});
}

void VCFRecord::confidence(){
    for (auto&& sample : regt_samples) {
        if (sample.find("REF_LIKELIHOOD") != sample.end() and sample.find("ALT_LIKELIHOOD") != sample.end()) {
            sample["GT_CONF"] = abs(sample["REF_LIKELIHOOD"] - sample["ALT_LIKELIHOOD"]);
        }
    }
    add_formats({"GT_CONF"});
}

void VCFRecord::regenotype(const uint8_t confidence_threshold){
    for (uint_least16_t i=0; i<samples.size(); ++i) {
        if (regt_samples[i].find("GT_CONF") != regt_samples[i].end()){
            if (regt_samples[i]["GT_CONF"] > confidence_threshold){
                if (samples[i]["GT"] == 0 and regt_samples[i]["ALT_LIKELIHOOD"] > regt_samples[i]["REF_LIKELIHOOD"]){
                    samples[i]["GT"] = 1;
                } else if (samples[i]["GT"] == 1 and regt_samples[i]["ALT_LIKELIHOOD"] < regt_samples[i]["REF_LIKELIHOOD"]){
                    samples[i]["GT"] = 0;
                    //swap_ref_and_alt_properties(i);
                }
            } else {
                samples[i].erase("GT");
            }
        } else {
            samples[i].erase("GT");
        }
    }
}

bool VCFRecord::operator==(const VCFRecord &y) const {
    if (chrom != y.chrom) { return false; }
    if (pos != y.pos) { return false; }
    if (ref != y.ref) { return false; }
    if (alt != y.alt) { return false; }
    return true;
}

bool VCFRecord::operator!=(const VCFRecord &y) const {
    return !(*this==y);
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

    for(uint_least16_t i=0;i<m.samples.size(); ++i){
        out << "\t";
        for (auto f : m.format){
            if (m.samples[i].find(f)!=m.samples[i].end()) {
                out << +m.samples.at(i).at(f);
            } else if (m.regt_samples.size() > 0 and m.regt_samples[i].find(f)!=m.regt_samples[i].end()) {
                out << +m.regt_samples.at(i).at(f);
            } else {
                out << ".";
            }

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
    vector<string> float_strings = {"REF_LIKELIHOOD", "ALT_LIKELIHOOD","GT_CONF"};
    unordered_map<string, uint8_t> sample_data;
    unordered_map<string, float> regt_sample_data;
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
        m.regt_samples.push_back(regt_sample_data);
        for (uint32_t i=0; i<m.format.size(); ++i){
            if (sample_strings[i] != "."
                and find(float_strings.begin(),float_strings.end(),m.format[i])==float_strings.end())
                m.samples.back()[m.format[i]] = stoi(sample_strings[i]);
            else if (sample_strings[i] != "."
                     and find(float_strings.begin(),float_strings.end(),m.format[i])!=float_strings.end())
                m.regt_samples.back()[m.format[i]] = stof(sample_strings[i]);
        }
    }
    return in;
}


