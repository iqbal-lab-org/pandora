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

void VCFRecord::set_format(const uint32_t& sample_id, const std::string& format, const std::vector<uint16_t>& val){
    assert(samples.size() > sample_id);
    samples[sample_id][format] = val;
    add_formats({format});
}

void VCFRecord::set_format(const uint32_t& sample_id, const std::string& format, const std::vector<float>& val){
    std::unordered_map<std::string, std::vector<float>> m;
    m.reserve(3);
    for (uint i=regt_samples.size(); i<samples.size(); ++i) {
        regt_samples.push_back(m);
    }
    assert(regt_samples.size() > sample_id);
    regt_samples[sample_id][format] = val;
    add_formats({format});
}

void VCFRecord::set_format(const uint32_t& sample_id, const std::string& format, const uint16_t& val){
    std::vector<uint16_t> v = {};
    v.emplace_back(val);
    set_format(sample_id, format, v);
}

void VCFRecord::set_format(const uint32_t& sample_id, const std::string& format, const uint32_t& val){
    assert(val < std::numeric_limits<uint16_t>::max());
    uint16_t v = val;
    set_format(sample_id, format, v);
}

void VCFRecord::set_format(const uint32_t& sample_id, const std::string& format, const float& val){
    std::vector<float> v = {};
    v.emplace_back(val);
    set_format(sample_id, format, v);
}

void VCFRecord::append_format(const uint32_t& sample_id, const std::string& format, const uint16_t& val){
    assert(samples.size() > sample_id);
    if (samples[sample_id].find(format)!=samples[sample_id].end()){
        samples[sample_id][format].push_back(val);
    } else {
        set_format(sample_id, format, val);
    }
}

void VCFRecord::append_format(const uint32_t& sample_id, const std::string& format, const uint32_t& val){
    assert(val < std::numeric_limits<uint16_t>::max());
    uint16_t v = val;
    append_format(sample_id, format, v);
}

void VCFRecord::append_format(const uint32_t& sample_id, const std::string& format, const float& val){
    std::unordered_map<std::string, std::vector<float>> m;
    m.reserve(3);
    if (regt_samples.empty()) {
        for (const auto &sample : samples) {
            regt_samples.push_back(m);
        }
    }
    assert(regt_samples.size() > sample_id);
    if (regt_samples[sample_id].find(format)!=regt_samples[sample_id].end()){
        regt_samples[sample_id][format].push_back(val);
    } else {
        set_format(sample_id, format, val);
    }
}

std::vector<uint16_t> VCFRecord::get_format_u(const uint32_t& sample_id, const std::string& format){
    std::vector<uint16_t> empty_return;
    bool sample_exists = samples.size() > sample_id;
    if (!sample_exists)
        return empty_return;
    bool found_key_in_sample = samples[sample_id].find(format)!=samples[sample_id].end();
    if (!found_key_in_sample)
        return empty_return;
    return samples[sample_id][format];
}

std::vector<float> VCFRecord::get_format_f(const uint32_t& sample_id, const std::string& format){
    std::vector<float> empty_return;
    bool sample_exists = regt_samples.size() > sample_id;
    if (!sample_exists)
        return empty_return;
    bool found_key_in_sample = regt_samples[sample_id].find(format)!=regt_samples[sample_id].end();
    if (!found_key_in_sample)
        return empty_return;
    return regt_samples[sample_id][format];
}

float logfactorial(uint32_t n) {
    float ret = 0;
    for (uint32_t i = 1; i <= n; ++i) {
        ret += log(i);
    }
    return ret;
}

void VCFRecord::likelihood(const uint32_t &expected_depth_covg, const float &error_rate, const uint32_t &min_covg) {
    for (uint_least16_t i = 0; i < samples.size(); ++i) {
        const auto &fwd_covgs = get_format_u(i,"MEAN_FWD_COVG");
        const auto &rev_covgs = get_format_u(i,"MEAN_REV_COVG");
        const auto &gaps = get_format_f(i,"GAPS");
        if (!fwd_covgs.empty() and fwd_covgs.size() == rev_covgs.size() and fwd_covgs.size() == gaps.size()){

            std::vector<uint16_t> covgs = {};
            for (uint j = 0; j < fwd_covgs.size(); ++j) {
                uint32_t total_covg = fwd_covgs[j] + rev_covgs[j];
                if (total_covg >= min_covg)
                    covgs.push_back(total_covg);
                else
                    covgs.push_back(0);
            }

            float likelihood = 0;
            for (uint j = 0; j < covgs.size(); ++j) {
                auto other_covg = accumulate(covgs.begin(), covgs.end(), 0) - covgs[j];
                if (covgs[j] > 0) {
                    likelihood = covgs[j] * log(expected_depth_covg) - expected_depth_covg
                                 - logfactorial(covgs[j]) + other_covg * log(error_rate);
                    //std::cout << "likelihood before gaps " << likelihood << " gaps = " << gaps[j] << std::endl;
                    likelihood += (1 - gaps[j]) * log(1 - exp(-(float)expected_depth_covg)) - expected_depth_covg * gaps[j];
                    //std::cout << "likelihood after gaps " << likelihood << std::endl;
                } else {
                    likelihood = other_covg * log(error_rate) - expected_depth_covg;
                    //std::cout << "likelihood before gaps " << likelihood << " gaps = " << gaps[j] << std::endl;
                    likelihood += (1 - gaps[j]) * log(1 - exp(-(float)expected_depth_covg)) - expected_depth_covg * gaps[j];
                    //std::cout << "likelihood after gaps " << likelihood << std::endl;
                }
                append_format(i,"LIKELIHOOD",likelihood);
            }
        }
    }

    assert(regt_samples.size() == samples.size() or assert_msg(regt_samples.size() << "!=" << samples.size()));
}

void VCFRecord::confidence(const uint32_t &min_total_covg, const uint32_t &min_diff_covg) {
    for (uint i=0; i < regt_samples.size(); ++i) {
        auto& sample = regt_samples[i];
        if (sample.find("LIKELIHOOD") != sample.end()) {
            assert(sample["LIKELIHOOD"].size() > 1);
            float max_lik = 0, max_lik2 = 0;
            uint32_t max_coord = 0, max_coord2 = 0;
            for (uint j=0; j < sample["LIKELIHOOD"].size(); ++j ){
                const auto & likelihood = sample["LIKELIHOOD"][j];
                if (max_lik == 0 or likelihood > max_lik) {
                    max_coord2 = max_coord;
                    max_coord = j;
                    max_lik2 = max_lik;
                    max_lik = likelihood;
                } else if (max_lik2 == 0 or likelihood > max_lik2) {
                    max_lik2 = likelihood;
                    max_coord2 = j;
                }
            }

            assert(samples.size() > i);
            assert(samples[i].find("MEAN_FWD_COVG")!=samples[i].end());
            assert(samples[i]["MEAN_FWD_COVG"].size() > max_coord);

            const auto & max_covg = samples[i]["MEAN_FWD_COVG"][max_coord]+samples[i]["MEAN_REV_COVG"][max_coord];
            const auto & next_covg = samples[i]["MEAN_FWD_COVG"][max_coord2]+samples[i]["MEAN_REV_COVG"][max_coord2];
            bool enough_total_covg = (max_covg + next_covg >= min_total_covg);
            bool enough_difference_in_covg = (std::abs(max_covg-next_covg) >= min_diff_covg);
            if (enough_total_covg and enough_difference_in_covg)
                sample["GT_CONF"] = {std::abs(max_lik - max_lik2)};
            else
                sample["GT_CONF"] = {0};
        }
    }
    add_formats({"GT_CONF"});
}

void VCFRecord::genotype(const uint16_t confidence_threshold) {
    for (uint_least16_t i = 0; i < samples.size(); ++i) {
        if (regt_samples[i].find("GT_CONF") != regt_samples[i].end()) {
            if (regt_samples[i]["GT_CONF"][0] > confidence_threshold) {
                uint16_t allele = 0;
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

bool VCFRecord::contains_dot_allele() const {
    if (ref == "." or ref == "")
        return true;
    for (const auto &a : alt){
        if (a == "." or a == "")
            return true;
    }
    return false;
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
    std::unordered_map<std::string, std::vector<uint16_t>> sample_data;
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


