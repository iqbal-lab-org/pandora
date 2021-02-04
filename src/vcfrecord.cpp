#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/log/trivial.hpp>
#include <vcfrecord.h>

VCFRecord::VCFRecord(VCF const* parent_vcf, const std::string& chrom, uint32_t pos,
    const std::string& ref, const std::string& alt, const std::string& info,
    const std::string& graph_type_info)
    : parent_vcf(parent_vcf)
    , id(".")
    , qual(".")
    , filter(".")
    , info(info)
    , chrom(chrom)
    , pos(pos)
{
    add_new_samples(parent_vcf->samples.size());
    set_ref_and_clear_alts(ref);
    add_new_alt(alt);

    if (this->info == ".") {
        this->info = infer_SVTYPE();
    }

    if (graph_type_info != "") {
        this->info += ";";
        this->info += graph_type_info;
    }
}

VCFRecord::VCFRecord(VCF const* parent_vcf)
    : parent_vcf(parent_vcf)
    , id(".")
    , qual(".")
    , filter(".")
    , info(".")
    , chrom(".")
    , pos(0)
{
    add_new_samples(parent_vcf->samples.size());
    set_ref_and_clear_alts(".");
}

std::string VCFRecord::infer_SVTYPE() const
{
    // TODO: How to handle cases where there are more than 2 options, not all of one
    // type
    if (ref == "." and (alts.empty() or alts[0] == "."))
        return ".";
    else if (ref == "." or alts.empty() or alts[0] == ".")
        return "SVTYPE=INDEL";
    else if (ref.length() == 1 and !alts.empty() and alts[0].length() == 1)
        return "SVTYPE=SNP";
    else if (!alts.empty() and alts[0].length() == ref.length())
        return "SVTYPE=PH_SNPs";
    else if (!alts.empty() and ref.length() < alts[0].length()
        and ref.compare(0, ref.length(), alts[0], 0, ref.length()) == 0)
        return "SVTYPE=INDEL";
    else if (!alts.empty() and alts[0].length() < ref.length()
        and alts[0].compare(0, alts[0].length(), ref, 0, alts[0].length()) == 0)
        return "SVTYPE=INDEL";
    else
        return "SVTYPE=COMPLEX";
}

std::string VCFRecord::get_format(
    bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const
{
    bool only_one_flag_is_set
        = ((int)(genotyping_from_maximum_likelihood) + (int)(genotyping_from_coverage))
        == 1;
    // this will still remain an assert as it is responsibility of the dev to ensure
    // this method is not called with the two flags set
    assert(only_one_flag_is_set);

    static std::vector<std::string> format_for_genotyping_from_maximum_likelihood
        = { "GT", "MEAN_FWD_COVG", "MEAN_REV_COVG", "MED_FWD_COVG", "MED_REV_COVG",
              "SUM_FWD_COVG", "SUM_REV_COVG", "GAPS" };
    static std::vector<std::string> format_for_genotyping_from_coverage
        = { "GT", "MEAN_FWD_COVG", "MEAN_REV_COVG", "MED_FWD_COVG", "MED_REV_COVG",
              "SUM_FWD_COVG", "SUM_REV_COVG", "GAPS", "LIKELIHOOD", "GT_CONF" };

    const std::vector<std::string>* format;
    if (genotyping_from_maximum_likelihood)
        format = &format_for_genotyping_from_maximum_likelihood;
    if (genotyping_from_coverage)
        format = &format_for_genotyping_from_coverage;

    std::stringstream out;
    for (const auto& field : *format) {
        out << field;
        bool is_not_last_field = field != format->back();
        if (is_not_last_field) {
            out << ":";
        }
    }

    return out.str();
}

bool VCFRecord::contains_dot_allele() const
{
    if (ref == "." or ref == "")
        return true;
    for (const auto& a : alts) {
        if (a == "." or a == "")
            return true;
    }
    return false;
}

bool VCFRecord::operator==(const VCFRecord& y) const
{
    if (chrom != y.chrom) {
        return false;
    }
    if (pos != y.pos) {
        return false;
    }
    if (ref != y.ref) {
        return false;
    }
    if (alts.size() != y.alts.size()) {
        return false;
    }
    for (const auto& a : alts) {
        if (find(y.alts.begin(), y.alts.end(), a) == y.alts.end()) {
            return false;
        }
    }
    return true;
}

bool VCFRecord::operator!=(const VCFRecord& y) const { return !(*this == y); }

bool VCFRecord::operator<(const VCFRecord& y) const
{
    if (chrom < y.chrom) {
        return true;
    }
    if (chrom > y.chrom) {
        return false;
    }
    if (pos < y.pos) {
        return true;
    }
    if (pos > y.pos) {
        return false;
    }
    if (ref < y.ref) {
        return true;
    }
    if (ref > y.ref) {
        return false;
    }
    if (alts < y.alts) {
        return true;
    }
    if (alts > y.alts) {
        return false;
    }
    return false;
}

std::string VCFRecord::to_string(
    bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const
{
    std::stringstream out;

    out << this->chrom << "\t" << this->pos + 1 << "\t" << this->id << "\t" << this->ref
        << "\t" << this->alts_to_string() << "\t" << this->qual << "\t" << this->filter
        << "\t" << this->info << "\t"
        << this->get_format(
               genotyping_from_maximum_likelihood, genotyping_from_coverage)
        << "\t"
        << this->sample_infos_to_string(
               genotyping_from_maximum_likelihood, genotyping_from_coverage);

    return out.str();
}

std::string VCFRecord::alts_to_string() const
{
    std::stringstream out;
    if (this->alts.empty()) {
        out << ".";
    } else {
        std::string buffer = "";
        for (const auto& alt : this->alts) {
            out << buffer << alt;
            buffer = ",";
        }
    }
    return out.str();
}

size_t VCFRecord::get_longest_allele_length() const
{
    size_t longest_allele_length = this->ref.size();
    for (const std::string& alt : this->alts) {
        longest_allele_length = std::max(longest_allele_length, alt.size());
    }
    return longest_allele_length;
}

void VCFRecord::merge_record_into_this(const VCFRecord& other)
{
    // no need for merge
    bool other_record_has_no_alt = other.alts.size() == 0;
    if (other_record_has_no_alt)
        return;

    if(!there_are_no_common_alt_alleles_between_this_and_other(other)) {
        fatal_error("When merging two VCF records, they have common ALTs, this should not happen");
    }

    this->sampleIndex_to_sampleInfo.merge_other_samples_infos_into_this(
        other.sampleIndex_to_sampleInfo);
    add_new_alts(other.get_alts().begin(), other.get_alts().end());
}

bool VCFRecord::can_biallelic_record_be_merged_into_this(
    const VCFRecord& vcf_record_to_be_merged_in, uint32_t max_allele_length) const
{
    // TODO : when this function is called, we could have VCFRecords with no alts, but
    // they would still be biallelic because these are corrected afterwards
    // TODO : maybe fix this?
    // bool ensure_we_are_merging_only_biallelic_records =
    // vcf_record_to_be_merged_in.alts.size() == 1;
    bool we_are_merging_only_biallelic_records
        = vcf_record_to_be_merged_in.alts.size() <= 1;
    if(!we_are_merging_only_biallelic_records) {
        fatal_error("When merging two biallelic records, one of them is not biallelic");
    }

    bool both_records_have_the_same_ref = this->ref == vcf_record_to_be_merged_in.ref;

    bool all_alleles_have_at_most_max_allele_length
        = this->get_longest_allele_length() <= max_allele_length
        and vcf_record_to_be_merged_in.get_longest_allele_length() <= max_allele_length;

    bool vcf_record_should_be_merged_in = vcf_record_to_be_merged_in != (*this)
        and this->has_the_same_position(vcf_record_to_be_merged_in)
        and both_records_have_the_same_ref
        and all_alleles_have_at_most_max_allele_length
        and there_are_no_common_alt_alleles_between_this_and_other(
            vcf_record_to_be_merged_in);

    return vcf_record_should_be_merged_in;
}

void VCFRecord::correct_dot_alleles(
    char nucleotide, bool add_nucleotide_before_the_sequence)
{
    std::string prefix = "";
    std::string suffix = "";
    int pos_change = 0;

    if (add_nucleotide_before_the_sequence) {
        prefix += nucleotide;
        pos_change = -1;
    } else {
        suffix += nucleotide;
        pos_change = 0;
    }

    if (allele_is_dot(get_ref())) {
        ref = nucleotide;
    } else {
        ref = prefix + ref + suffix;
    }

    for (auto& alt : alts) {
        if (allele_is_dot(alt)) {
            alt = nucleotide;
        } else {
            alt = prefix + alt + suffix;
        }
    }

    pos += pos_change;
}

void VCFRecord::set_ref_and_clear_alts(std::string ref)
{
    if (ref == "") {
        ref = ".";
    }
    this->ref = ref;
    this->alts.clear();
    set_number_of_alleles_and_resize_coverage_information_for_all_samples(1);
}

void VCFRecord::add_new_alt(std::string alt)
{
    if (alt == "") {
        alt = ".";
    }

    bool alt_already_present = std::find(alts.begin(), alts.end(), alt) != alts.end();
    if (alt_already_present) {
        fatal_error("Error adding new ALT to a VCF record: ALT already exists");
    }

    alts.push_back(alt);
    set_number_of_alleles_and_resize_coverage_information_for_all_samples(
        alts.size() + 1);
}

void VCFRecord::add_new_samples(uint32_t number_of_samples)
{
    sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(
        number_of_samples, get_number_of_alleles(), parent_vcf->genotyping_options);
}

bool VCFRecord::there_are_no_common_alt_alleles_between_this_and_other(
    const VCFRecord& other) const
{
    std::set<std::string> all_unique_alt_alleles;
    all_unique_alt_alleles.insert(this->get_alts().begin(), this->get_alts().end());
    all_unique_alt_alleles.insert(other.get_alts().begin(), other.get_alts().end());
    size_t the_size_we_are_supposed_to_have_if_there_are_no_common_alt_alleles
        = this->get_alts().size() + other.get_alts().size();
    return all_unique_alt_alleles.size()
        == the_size_we_are_supposed_to_have_if_there_are_no_common_alt_alleles;
}

void VCFRecord::reset_sample_infos_to_contain_the_given_number_of_samples(
    uint32_t number_of_samples)
{
    sampleIndex_to_sampleInfo.clear();
    sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(
        number_of_samples, get_number_of_alleles(), parent_vcf->genotyping_options);
}