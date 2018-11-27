#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include <boost/log/trivial.hpp>

#include "vcfrecord.h"
#include "utils.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


VCFRecord::VCFRecord(std::string c, uint32_t p, std::string r, std::string a, std::string i, std::string g) : chrom(c),
                                                                                                              pos(p),
                                                                                                              id("."),
                                                                                                              ref(r),
                                                                                                              qual("."),
                                                                                                              filter("."),
                                                                                                              info(i),
                                                                                                              format({"GT"}) {
    if (a == "")
        alt.push_back(".");
    else
        alt.push_back(a);

    // fix so no empty strings
    if (ref == "") { ref = "."; }

    // classify variants as SNPs, INDELs PH_SNPs or COMPLEX
    // need to think about how to handle cases where there are more than 2 options, not all of one type
    if (info == ".") {
        if (ref == "." and (alt.empty() or alt[0] == ".")) {}
        else if (ref == "." or alt.empty() or alt[0] == ".") { info = "SVTYPE=INDEL"; }
        else if (ref.length() == 1 and !alt.empty() and alt[0].length() == 1) { info = "SVTYPE=SNP"; }
        else if (!alt.empty() and alt[0].length() == ref.length()) { info = "SVTYPE=PH_SNPs"; }
        else if (!alt.empty() and ref.length() < alt[0].length() and
                 ref.compare(0, ref.length(), alt[0], 0, ref.length()) == 0) { info = "SVTYPE=INDEL"; }
        else if (!alt.empty() and alt[0].length() < ref.length() and
                 alt[0].compare(0, alt[0].length(), ref, 0, alt[0].length()) == 0) { info = "SVTYPE=INDEL"; }
        else { info = "SVTYPE=COMPLEX"; }
    }

    // add graph type info
    if (g != "") {
        info += ";";
        info += g;
    }
};

VCFRecord::VCFRecord() : chrom("."), pos(0), id("."), ref("."), qual("."), filter("."), info(".") {};

VCFRecord::VCFRecord(const VCFRecord &other) {
    chrom = other.chrom;
    pos = other.pos;
    id = other.id;
    ref = other.ref;
    alt = other.alt;
    qual = other.qual;
    filter = other.filter;
    info = other.info;
    format = other.format;
    samples = other.samples;
    regt_samples = other.regt_samples;
}

VCFRecord &VCFRecord::operator=(const VCFRecord &other) {
    // check for self-assignment
    if (this == &other)
        return *this;

    chrom = other.chrom;
    pos = other.pos;
    id = other.id;
    ref = other.ref;
    alt = other.alt;
    qual = other.qual;
    filter = other.filter;
    info = other.info;
    format = other.format;
    samples = other.samples;
    regt_samples = other.regt_samples;

    return *this;
}

VCFRecord::~VCFRecord() {};

void VCFRecord::clear() {
    chrom = ".";
    pos = 0;
    id = ".";
    ref = ".";
    alt.clear();
    qual = ".";
    filter = ".";
    info = ".";
    format.clear();
    for (auto s : samples)
        s.clear();
    samples.clear();
    for (auto r : regt_samples)
        r.clear();
    regt_samples.clear();
}

void VCFRecord::clear_sample(uint32_t i) {
    if (samples.size() > i) {
        samples.at(i).clear();
    }

    if (regt_samples.size() > i) {
        regt_samples.at(i).clear();
    }
    bool all_cleared(true);
    for (const auto &s : samples) {
        if (!s.empty()) {
            all_cleared = false;
            break;
        }
    }
    if (all_cleared) {
        clear();
    }
}

void VCFRecord::add_formats(const std::vector<std::string> &formats) {
    for (const auto &s : formats) {
        if (find(format.begin(), format.end(), s) == format.end())
            format.push_back(s);
    }
}

float logfactorial(uint32_t n) {
    float ret = 0;
    for (uint32_t i = 1; i <= n; ++i) {
        ret += log(i);
    }
    return ret;
}

void VCFRecord::likelihood(const uint32_t &expected_depth_covg, const float &error_rate) {
    std::unordered_map<std::string, std::vector<float>> m;
    m.reserve(2);

    //float p_non_zero = 1 - exp(-expected_depth_covg);
    if (regt_samples.empty()) {
        for (const auto &sample : samples) {
            regt_samples.push_back(m);
        }
    }

    for (uint_least16_t i = 0; i < samples.size(); ++i) {
        if (samples[i].find("MEAN_FWD_COVG") != samples[i].end()
            and samples[i].find("MEAN_REV_COVG") != samples[i].end()
            and samples[i]["MEAN_FWD_COVG"].size() == samples[i]["MEAN_REV_COVG"].size()
            and samples[i]["MEAN_FWD_COVG"].size() >= 2) {

            std::vector<uint16_t> covgs = {};
            for (uint j = 0; j < samples[i]["MEAN_FWD_COVG"].size(); ++j) {
                covgs.push_back(samples[i]["MEAN_FWD_COVG"][j] + samples[i]["MEAN_REV_COVG"][j]);
            }

            std::vector<float> likelihoods = {};
            float likelihood = 0;
            for (uint j = 0; j < covgs.size(); ++j) {
                auto other_covg = accumulate(covgs.begin(), covgs.end(), 0) - covgs[j];
                if (covgs[j] > 0)
                    likelihood = covgs[j] * log(expected_depth_covg) - expected_depth_covg
                                 - logfactorial(covgs[j]) + other_covg * log(error_rate);
                else
                    likelihood = other_covg * log(error_rate) - expected_depth_covg;
                likelihoods.push_back(likelihood);
            }
            regt_samples[i]["LIKELIHOOD"] = likelihoods;
        }
    }

    assert(regt_samples.size() == samples.size() or assert_msg(regt_samples.size() << "!=" << samples.size()));
    add_formats({"LIKELIHOOD"});
}

void VCFRecord::confidence() {
    for (auto &sample : regt_samples) {
        if (sample.find("LIKELIHOOD") != sample.end()) {
            assert(sample["LIKELIHOOD"].size() > 1);
            float max_lik = 0, max_lik2 = 0;
            for (const auto &likelihood : sample["LIKELIHOOD"]) {
                if (max_lik == 0 or likelihood > max_lik) {
                    max_lik2 = max_lik;
                    max_lik = likelihood;
                } else if (max_lik2 == 0 or likelihood > max_lik2) {
                    max_lik2 = likelihood;
                }
            }
            sample["GT_CONF"] = {std::abs(max_lik - max_lik2)};
        }
    }
    add_formats({"GT_CONF"});
}

void VCFRecord::genotype(const uint8_t confidence_threshold) {
    for (uint_least16_t i = 0; i < samples.size(); ++i) {
        if (regt_samples[i].find("GT_CONF") != regt_samples[i].end()) {
            if (regt_samples[i]["GT_CONF"][0] > confidence_threshold) {
                uint8_t allele = 0;
                float max_likelihood = 0;
                for (const auto &likelihood : regt_samples[i]["LIKELIHOOD"]) {
                    if (max_likelihood == 0 or likelihood > max_likelihood) {
                        samples[i]["GT"] = {allele};
                        max_likelihood = likelihood;
                    }
                    allele++;
                }
            } else {
                samples[i]["GT"].clear();
            }
        } else {
            samples[i]["GT"].clear();
        }
    }
}

bool VCFRecord::operator==(const VCFRecord &y) const {
    if (chrom != y.chrom) { return false; }
    if (pos != y.pos) { return false; }
    if (ref != y.ref) { return false; }
    if (alt.size() != y.alt.size()) { return false; }
    for (const auto &a : alt) {
        if (find(y.alt.begin(), y.alt.end(), a) == y.alt.end()) { return false; }
    }
    return true;
}

bool VCFRecord::operator!=(const VCFRecord &y) const {
    return !(*this == y);
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
    out << m.chrom << "\t" << m.pos + 1 << "\t" << m.id << "\t" << m.ref << "\t";

    if (m.alt.empty()) {
        out << ".";
    } else {
        std::string buffer = "";
        for (const auto &a : m.alt) {
            out << buffer << a;
            buffer = ",";
        }
    }
    out << "\t" << m.qual << "\t"
        << m.filter << "\t" << m.info << "\t";

    std::string last_format;
    if (!m.format.empty())
        last_format = m.format[m.format.size() - 1];

    for (const auto &s : m.format) {
        out << s;
        if (s != last_format)
            out << ":";
    }

    for (uint_least16_t i = 0; i < m.samples.size(); ++i) {
        out << "\t";
        for (const auto &f : m.format) {
            std::string buffer = "";
            if (m.samples[i].find(f) != m.samples[i].end() and not m.samples[i].at(f).empty()) {
                for (const auto &a : m.samples.at(i).at(f)) {
                    out << buffer << +a;
                    buffer = ",";
                }

            } else if (m.regt_samples.size() > i
                       and m.regt_samples[i].find(f) != m.regt_samples[i].end()
                       and not m.regt_samples[i].at(f).empty()) {
                for (const auto &a : m.regt_samples.at(i).at(f)) {
                    out << buffer << +a;
                    buffer = ",";
                }
            } else {
                out << ".";
            }

            if (f != last_format)
                out << ":";
        }
    }

    out << std::endl;
    return out;
}

std::istream &operator>>(std::istream &in, VCFRecord &m) {
    std::string token, alt_s;
    std::vector<std::string> sample_strings, sample_substrings;
    std::vector<std::string> float_strings = {"LIKELIHOOD", "GT_CONF"};
    std::unordered_map<std::string, std::vector<uint8_t>> sample_data;
    std::unordered_map<std::string, std::vector<float>> regt_sample_data;
    m.alt.clear();
    in >> m.chrom;
    in.ignore(1, '\t');
    in >> m.pos;
    m.pos -= 1;
    in.ignore(1, '\t');
    in >> m.id;
    in.ignore(1, '\t');
    in >> m.ref;
    in.ignore(1, '\t');
    in >> alt_s;
    m.alt.push_back(alt_s);
    int c = in.peek();
    while (c == ',') {
        in >> alt_s;
        m.alt.push_back(alt_s);
        c = in.peek();
    }
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
        assert(sample_strings.size() == m.format.size() or assert_msg("sample data does not fit format"));
        m.samples.push_back(sample_data);
        m.regt_samples.push_back(regt_sample_data);
        for (uint32_t i = 0; i < m.format.size(); ++i) {
            if (sample_strings[i] != "."
                and find(float_strings.begin(), float_strings.end(), m.format[i]) == float_strings.end()) {
                sample_substrings = split(sample_strings[i], ",");
                for (const auto &s : sample_substrings)
                    m.samples.back()[m.format[i]].push_back(stoi(s));
            } else if (sample_strings[i] != "."
                       and find(float_strings.begin(), float_strings.end(), m.format[i]) != float_strings.end()) {
                sample_substrings = split(sample_strings[i], ",");
                for (const auto &s : sample_substrings)
                    m.regt_samples.back()[m.format[i]].push_back(stof(s));
            }
        }
    }
    return in;
}


