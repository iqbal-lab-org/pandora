#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <boost/log/trivial.hpp>

#include "vcfrecord.h"
#include "utils.h"
#include "sampleinfo.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


VCFRecord::VCFRecord(const std::string &chrom, uint32_t pos, const std::string &ref, const std::string &alt,
                     const std::string &info, const std::string &graph_type_info) : chrom(chrom),
                                                                                    pos(pos),
                                                                                    id("."),
                                                                                    ref(ref),
                                                                                    qual("."),
                                                                                    filter("."),
                                                                                    info(info),
                                                                                    format({"GT"}) {
    if (alt == "")
        this->alts.push_back(".");
    else
        this->alts.push_back(alt);

    if (this->ref == "") {
        this->ref = ".";
    }

    if (this->info == ".") {
        this->info = infer_SVTYPE();
    }

    if (graph_type_info != "") {
        this->info += ";";
        this->info += graph_type_info;
    }
}

VCFRecord::VCFRecord() : chrom("."), pos(0), id("."), ref("."), qual("."), filter("."), info(".") {}

std::string VCFRecord::infer_SVTYPE() const {
    // TODO: How to handle cases where there are more than 2 options, not all of one type
    if (ref == "." and (alts.empty() or alts[0] == "."))
        return ".";
    else if (ref == "." or alts.empty() or alts[0] == ".")
        return "SVTYPE=INDEL";
    else if (ref.length() == 1 and !alts.empty() and alts[0].length() == 1)
        return "SVTYPE=SNP";
    else if (!alts.empty() and alts[0].length() == ref.length())
        return "SVTYPE=PH_SNPs";
    else if (!alts.empty() and ref.length() < alts[0].length() and
             ref.compare(0, ref.length(), alts[0], 0, ref.length()) == 0)
        return "SVTYPE=INDEL";
    else if (!alts.empty() and alts[0].length() < ref.length() and
             alts[0].compare(0, alts[0].length(), ref, 0, alts[0].length()) == 0)
        return "SVTYPE=INDEL";
    else
        return "SVTYPE=COMPLEX";
}

void VCFRecord::clear_sample(uint32_t i) {
    if (sampleIndex_to_format_to_sampleInfo.size() > i) {
        sampleIndex_to_format_to_sampleInfo.at(i).clear();
    }

    if (sampleIndex_to_format_to_sampleGenotypedInfo.size() > i) {
        sampleIndex_to_format_to_sampleGenotypedInfo.at(i).clear();
    }
    bool all_cleared(true);
    for (const auto &s : sampleIndex_to_format_to_sampleInfo) {
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
    assert(sampleIndex_to_format_to_sampleInfo.size() > sample_id);
    sampleIndex_to_format_to_sampleInfo[sample_id][format] = val;
    add_formats({format});
}

void VCFRecord::set_format(const uint32_t& sample_id, const std::string& format, const std::vector<float>& val){
    size_t amount_to_push = 0;
    if (sampleIndex_to_format_to_sampleInfo.size() > sampleIndex_to_format_to_sampleGenotypedInfo.size())
        amount_to_push = sampleIndex_to_format_to_sampleInfo.size() - sampleIndex_to_format_to_sampleGenotypedInfo.size();
    sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(amount_to_push);
    assert(sampleIndex_to_format_to_sampleGenotypedInfo.size() > sample_id);
    sampleIndex_to_format_to_sampleGenotypedInfo[sample_id][format] = val;
    add_formats({format});
}

void VCFRecord::set_format(const uint32_t& sample_id, const std::string& format, const uint16_t& val){
    std::vector<uint16_t> v = {};
    v.emplace_back(val);
    set_format(sample_id, format, v);
}

void VCFRecord::set_format(const uint32_t& sample_id, const std::string& format, const uint32_t& val){
    if (val >= std::numeric_limits<uint16_t>::max()){
        BOOST_LOG_TRIVIAL(debug) << now() << "Value too large for sample id " << sample_id << " format " << format
                                 << " value " << val;
        uint16_t v = std::numeric_limits<uint16_t>::max() - 1;
        set_format(sample_id, format, v);
    } else {
        uint16_t v = val;
        set_format(sample_id, format, v);
    }
}

void VCFRecord::set_format(const uint32_t& sample_id, const std::string& format, const float& val){
    std::vector<float> v = {};
    v.emplace_back(val);
    set_format(sample_id, format, v);
}

void VCFRecord::append_format(const uint32_t& sample_id, const std::string& format, const uint16_t& val){
    assert(sampleIndex_to_format_to_sampleInfo.size() > sample_id);
    if (sampleIndex_to_format_to_sampleInfo[sample_id].find(format) != sampleIndex_to_format_to_sampleInfo[sample_id].end()){
        sampleIndex_to_format_to_sampleInfo[sample_id][format].push_back(val);
    } else {
        set_format(sample_id, format, val);
    }
}

void VCFRecord::append_format(const uint32_t& sample_id, const std::string& format, const uint32_t& val){
    if (val >= std::numeric_limits<uint16_t>::max()){
        BOOST_LOG_TRIVIAL(debug) << now() << "Value too large for sample id " << sample_id << " format " << format
                                 << " value " << val;
        uint16_t v = std::numeric_limits<uint16_t>::max() - 1;
        append_format(sample_id, format, v);
    } else {
        uint16_t v = val;
        append_format(sample_id, format, v);
    }
}

void VCFRecord::append_format(const uint32_t& sample_id, const std::string& format, const float& val){
    if (sampleIndex_to_format_to_sampleGenotypedInfo.empty()) {
        size_t amount_to_push = sampleIndex_to_format_to_sampleInfo.size();
        sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(amount_to_push);
    }

    assert(sampleIndex_to_format_to_sampleGenotypedInfo.size() > sample_id);
    if (sampleIndex_to_format_to_sampleGenotypedInfo[sample_id].find(format) != sampleIndex_to_format_to_sampleGenotypedInfo[sample_id].end()){
        sampleIndex_to_format_to_sampleGenotypedInfo[sample_id][format].push_back(val);
    } else {
        set_format(sample_id, format, val);
    }
}

std::vector<uint16_t> VCFRecord::get_format_u(const uint32_t& sample_id, const std::string& format){
    std::vector<uint16_t> empty_return;
    bool sample_exists = sampleIndex_to_format_to_sampleInfo.size() > sample_id;
    if (!sample_exists)
        return empty_return;
    bool found_key_in_sample = sampleIndex_to_format_to_sampleInfo[sample_id].find(format) != sampleIndex_to_format_to_sampleInfo[sample_id].end();
    if (!found_key_in_sample)
        return empty_return;
    return sampleIndex_to_format_to_sampleInfo[sample_id][format];
}

std::vector<float> VCFRecord::get_format_f(const uint32_t& sample_id, const std::string& format){
    std::vector<float> empty_return;
    bool sample_exists = sampleIndex_to_format_to_sampleGenotypedInfo.size() > sample_id;
    if (!sample_exists)
        return empty_return;
    bool found_key_in_sample = sampleIndex_to_format_to_sampleGenotypedInfo[sample_id].find(format) != sampleIndex_to_format_to_sampleGenotypedInfo[sample_id].end();
    if (!found_key_in_sample)
        return empty_return;
    return sampleIndex_to_format_to_sampleGenotypedInfo[sample_id][format];
}

float logfactorial(uint32_t n) {
    float ret = 0;
    for (uint32_t i = 1; i <= n; ++i) {
        ret += log(i);
    }
    return ret;
}

void VCFRecord::likelihood(const std::vector<uint32_t> &expected_depth_covg_v, const float &error_rate,
                           const uint32_t &min_allele_covg, const float &min_fraction_allele_covg) {
    for (uint_least16_t i = 0; i < sampleIndex_to_format_to_sampleInfo.size(); ++i) {
        assert(i < expected_depth_covg_v.size());
        const auto &expected_depth_covg = expected_depth_covg_v[i];
        const uint32_t min_covg = std::max(min_allele_covg, uint(min_fraction_allele_covg*expected_depth_covg));
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

    assert(sampleIndex_to_format_to_sampleGenotypedInfo.size() == sampleIndex_to_format_to_sampleInfo.size() or assert_msg(sampleIndex_to_format_to_sampleGenotypedInfo.size() << "!=" << sampleIndex_to_format_to_sampleInfo.size()));
}

void VCFRecord::confidence(const uint32_t &min_total_covg, const uint32_t &min_diff_covg) {
    for (uint i=0; i < sampleIndex_to_format_to_sampleGenotypedInfo.size(); ++i) {
        auto& sample = sampleIndex_to_format_to_sampleGenotypedInfo[i];
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

            assert(sampleIndex_to_format_to_sampleInfo.size() > i);
            assert(sampleIndex_to_format_to_sampleInfo[i].find("MEAN_FWD_COVG") != sampleIndex_to_format_to_sampleInfo[i].end());
            assert(sampleIndex_to_format_to_sampleInfo[i]["MEAN_FWD_COVG"].size() > max_coord);

            const auto & max_covg = sampleIndex_to_format_to_sampleInfo[i]["MEAN_FWD_COVG"][max_coord] + sampleIndex_to_format_to_sampleInfo[i]["MEAN_REV_COVG"][max_coord];
            const auto & next_covg = sampleIndex_to_format_to_sampleInfo[i]["MEAN_FWD_COVG"][max_coord2] + sampleIndex_to_format_to_sampleInfo[i]["MEAN_REV_COVG"][max_coord2];
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
    for (uint_least16_t i = 0; i < sampleIndex_to_format_to_sampleInfo.size(); ++i) {
        if (sampleIndex_to_format_to_sampleGenotypedInfo[i].find("GT_CONF") != sampleIndex_to_format_to_sampleGenotypedInfo[i].end()) {
            if (sampleIndex_to_format_to_sampleGenotypedInfo[i]["GT_CONF"][0] > confidence_threshold) {
                uint16_t allele = 0;
                float max_likelihood = 0;
                for (const auto &likelihood : sampleIndex_to_format_to_sampleGenotypedInfo[i]["LIKELIHOOD"]) {
                    if (max_likelihood == 0 or likelihood > max_likelihood) {
                        sampleIndex_to_format_to_sampleInfo[i]["GT"] = {allele};
                        max_likelihood = likelihood;
                    }
                    allele++;
                }
            } else {
                sampleIndex_to_format_to_sampleInfo[i]["GT"].clear();
            }
        } else {
            sampleIndex_to_format_to_sampleInfo[i]["GT"].clear();
        }
    }
}

bool VCFRecord::contains_dot_allele() const {
    if (ref == "." or ref == "")
        return true;
    for (const auto &a : alts){
        if (a == "." or a == "")
            return true;
    }
    return false;
}

bool VCFRecord::operator==(const VCFRecord &y) const {
    if (chrom != y.chrom) { return false; }
    if (pos != y.pos) { return false; }
    if (ref != y.ref) { return false; }
    if (alts.size() != y.alts.size()) { return false; }
    for (const auto &a : alts) {
        if (find(y.alts.begin(), y.alts.end(), a) == y.alts.end()) { return false; }
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
    if (alts < y.alts) { return true; }
    if (alts > y.alts) { return false; }
    return false;
}


std::string VCFRecord::to_string() const {
    std::stringstream out;
    out << this->chrom << "\t" << this->pos + 1 << "\t" << this->id << "\t" << this->ref << "\t";

    if (this->alts.empty()) {
        out << ".";
    } else {
        std::string buffer = "";
        for (const auto &a : this->alts) {
            out << buffer << a;
            buffer = ",";
        }
    }
    out << "\t" << this->qual << "\t"
        << this->filter << "\t" << this->info << "\t";

    std::string last_format;
    if (!this->format.empty())
        last_format = this->format[this->format.size() - 1];

    for (const auto &s : this->format) {
        out << s;
        if (s != last_format)
            out << ":";
    }

    for (uint_least16_t i = 0; i < this->sampleIndex_to_format_to_sampleInfo.size(); ++i) {
        out << "\t";
        for (const auto &f : this->format) {
            std::string buffer = "";
            if (this->sampleIndex_to_format_to_sampleInfo[i].find(f) != this->sampleIndex_to_format_to_sampleInfo[i].end() and not this->sampleIndex_to_format_to_sampleInfo[i].at(f).empty()) {
                for (const auto &a : this->sampleIndex_to_format_to_sampleInfo.at(i).at(f)) {
                    out << buffer << +a;
                    buffer = ",";
                }

            } else if (this->sampleIndex_to_format_to_sampleGenotypedInfo.size() > i
                       and this->sampleIndex_to_format_to_sampleGenotypedInfo[i].find(f) != this->sampleIndex_to_format_to_sampleGenotypedInfo[i].end()
                       and not this->sampleIndex_to_format_to_sampleGenotypedInfo[i].at(f).empty()) {
                for (const auto &a : this->sampleIndex_to_format_to_sampleGenotypedInfo.at(i).at(f)) {
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

    return out.str();
}


std::ostream &operator<<(std::ostream &out, VCFRecord const &vcf_record) {
    out << vcf_record.to_string();
    return out;
}

std::istream &operator>>(std::istream &in, VCFRecord &m) {
    std::string token, alt_s;
    std::vector<std::string> sample_strings, sample_substrings;
    std::vector<std::string> float_strings = {"LIKELIHOOD", "GT_CONF", "GAPS"};
    m.alts.clear();
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
    m.alts.push_back(alt_s);
    int c = in.peek();
    while (c == ',') {
        in >> alt_s;
        m.alts.push_back(alt_s);
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
    while (in >> token) {
        sample_strings = split(token, ":");
        assert(sample_strings.size() == m.format.size() or assert_msg("sample data does not fit format"));
        m.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
        m.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
        for (uint32_t i = 0; i < m.format.size(); ++i) {
            if (sample_strings[i] != "."
                and find(float_strings.begin(), float_strings.end(), m.format[i]) == float_strings.end()) {
                sample_substrings = split(sample_strings[i], ",");
                for (const auto &s : sample_substrings)
                    m.sampleIndex_to_format_to_sampleInfo.back()[m.format[i]].push_back(stoi(s));
            } else if (sample_strings[i] != "."
                       and find(float_strings.begin(), float_strings.end(), m.format[i]) != float_strings.end()) {
                sample_substrings = split(sample_strings[i], ",");
                for (const auto &s : sample_substrings)
                    m.sampleIndex_to_format_to_sampleGenotypedInfo.back()[m.format[i]].push_back(stof(s));
            }
        }
    }
    return in;
}


size_t VCFRecord::get_longest_allele_length() const {
    size_t longest_allele_length = this->ref.size();
    for (const std::string &alt : this->alts) {
        longest_allele_length = std::max(longest_allele_length, alt.size());
    }
    return longest_allele_length;

}


void VCFRecord::merge_sample_information(const VCFRecord &other) {
    for (size_t sample_index = 0; sample_index < this->sampleIndex_to_format_to_sampleInfo.size(); ++sample_index) {
        auto keys = {"MEAN_FWD_COVG", "MEAN_REV_COVG",
                     "MED_FWD_COVG", "MED_REV_COVG",
                     "SUM_FWD_COVG", "SUM_REV_COVG"};
        for (const auto &key: keys) {
            merge_sample_key(this->sampleIndex_to_format_to_sampleInfo[sample_index], other.sampleIndex_to_format_to_sampleInfo[sample_index], key);
        }

        keys = {"LIKELIHOOD", "GT_CONF", "GAPS"};
        if (!this->sampleIndex_to_format_to_sampleGenotypedInfo.empty() and !other.sampleIndex_to_format_to_sampleGenotypedInfo.empty()) {
            for (const auto &key: keys)
                merge_sample_key(this->sampleIndex_to_format_to_sampleGenotypedInfo[sample_index], other.sampleIndex_to_format_to_sampleGenotypedInfo[sample_index], key);
        }
    }
}


void VCFRecord::merge_gt(const VCFRecord &other) {
    for (size_t sample_index = 0; sample_index < this->sampleIndex_to_format_to_sampleInfo.size(); ++sample_index) {
        if (other.sampleIndex_to_format_to_sampleInfo[sample_index].find("GT") == other.sampleIndex_to_format_to_sampleInfo[sample_index].end()
            or other.sampleIndex_to_format_to_sampleInfo[sample_index].at("GT").empty()) {
            continue;
        } else if (this->sampleIndex_to_format_to_sampleInfo[sample_index].find("GT") == this->sampleIndex_to_format_to_sampleInfo[sample_index].end()
                   or this->sampleIndex_to_format_to_sampleInfo[sample_index]["GT"].empty()) {
            if (other.sampleIndex_to_format_to_sampleInfo[sample_index].at("GT")[0] == 0) {
                this->sampleIndex_to_format_to_sampleInfo[sample_index]["GT"] = {0};
            } else {
                uint32_t allele_offset = this->alts.size();
                uint16_t new_allele = other.sampleIndex_to_format_to_sampleInfo[sample_index].at("GT")[0] + allele_offset;
                this->sampleIndex_to_format_to_sampleInfo[sample_index]["GT"] = {new_allele};
            }
        } else if (this->sampleIndex_to_format_to_sampleInfo[sample_index]["GT"][0] != 0 or other.sampleIndex_to_format_to_sampleInfo[sample_index].at("GT")[0] != 0) {
            //conflict, try to resolve with likelihoods
            if (this->sampleIndex_to_format_to_sampleGenotypedInfo.size() > sample_index
                and this->sampleIndex_to_format_to_sampleGenotypedInfo[sample_index].find("LIKELIHOOD") != this->sampleIndex_to_format_to_sampleGenotypedInfo[sample_index].end()) {
                this->confidence();
                this->genotype(5); //TODO: solve this bug
            } else {
                this->sampleIndex_to_format_to_sampleInfo[sample_index]["GT"] = {};
            }
        }
    }
}