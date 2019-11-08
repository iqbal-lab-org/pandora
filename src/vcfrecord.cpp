#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/log/trivial.hpp>
#include "vcf.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


VCFRecord::VCFRecord(const std::string &chrom, uint32_t pos, const std::string &ref, const std::string &alt,
                     const std::string &info, const std::string &graph_type_info) : chrom(chrom),
                                                                                    pos(pos),
                                                                                    id("."),
                                                                                    ref(ref),
                                                                                    qual("."),
                                                                                    filter("."),
                                                                                    info(info) {
    if (alt != "")
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


std::string VCFRecord::get_format (bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const {
    bool only_one_flag_is_set = ((int)(genotyping_from_maximum_likelihood) + (int)(genotyping_from_coverage)) == 1;
    assert(only_one_flag_is_set);

    static std::vector<std::string> format_for_genotyping_from_maximum_likelihood = {"GT", "MEAN_FWD_COVG", "MEAN_REV_COVG", "MED_FWD_COVG", "MED_REV_COVG", "SUM_FWD_COVG", "SUM_REV_COVG", "GAPS"};
    static std::vector<std::string> format_for_genotyping_from_coverage = {"GT", "MEAN_FWD_COVG", "MEAN_REV_COVG", "MED_FWD_COVG", "MED_REV_COVG", "SUM_FWD_COVG", "SUM_REV_COVG", "GAPS", "LIKELIHOOD", "GT_CONF"};

    const std::vector<std::string> *format;
    if (genotyping_from_maximum_likelihood)
        format = &format_for_genotyping_from_maximum_likelihood;
    if (genotyping_from_coverage)
        format = &format_for_genotyping_from_coverage;

    std::stringstream out;
    for (const auto &field : *format) {
        out << field;
        bool is_not_last_field = field != format->back();
        if (is_not_last_field) {
            out << ":";
        }
    }

    return out.str();
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


std::string VCFRecord::to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const {
    std::stringstream out;

    out << this->chrom << "\t"
        << this->pos + 1 << "\t"
        << this->id << "\t"
        << this->ref << "\t"
        << this->alts_to_string() << "\t"
        << this->qual << "\t"
        << this->filter << "\t"
        << this->info << "\t"
        << this->get_format(genotyping_from_maximum_likelihood, genotyping_from_coverage) << "\t"
        << this->sample_infos_to_string(genotyping_from_maximum_likelihood, genotyping_from_coverage);

    return out.str();
}

std::string VCFRecord::alts_to_string() const {
    std::stringstream out;
    if (this->alts.empty()) {
        out << ".";
    } else {
        std::string buffer = "";
        for (const auto &alt : this->alts) {
            out << buffer << alt;
            buffer = ",";
        }
    }
    return out.str();
}


// TODO: check if we keep this, it is only used in tests - better to keep in a VCFMock class
/*
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
        m.sampleIndex_to_format_to_sampleGenotypedInfo.emplace_back_several_empty_sample_infos(1);
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
*/

size_t VCFRecord::get_longest_allele_length() const {
    size_t longest_allele_length = this->ref.size();
    for (const std::string &alt : this->alts) {
        longest_allele_length = std::max(longest_allele_length, alt.size());
    }
    return longest_allele_length;

}


bool VCFRecord::can_biallelic_record_be_merged_into_this (const VCFRecord &vcf_record_to_be_merged_in, uint32_t max_allele_length) const {
    bool ensure_we_are_merging_only_biallelic_records = vcf_record_to_be_merged_in.alts.size() == 1;
    assert(ensure_we_are_merging_only_biallelic_records);

    bool both_records_have_the_same_ref = this->has_non_null_reference()
                                          and vcf_record_to_be_merged_in.has_non_null_reference()
                                          and this->ref == vcf_record_to_be_merged_in.ref;

    bool all_alleles_have_at_most_max_allele_length =
            this->get_longest_allele_length() <= max_allele_length
            and vcf_record_to_be_merged_in.get_longest_allele_length() <= max_allele_length;

    bool vcf_record_should_be_merged_in = vcf_record_to_be_merged_in != (*this)
                                          and this->has_the_same_position(vcf_record_to_be_merged_in)
                                          and both_records_have_the_same_ref
                                          and all_alleles_have_at_most_max_allele_length;

    return vcf_record_should_be_merged_in;
}

