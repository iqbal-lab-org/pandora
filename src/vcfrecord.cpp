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