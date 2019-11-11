#ifndef __VCFRECORD_H_INCLUDED__   // if vcfrecord.h hasn't been included yet...
#define __VCFRECORD_H_INCLUDED__

#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include "sampleinfo.h"
#include "vcf.h"

class VCF;

class VCFRecord {
public:
    VCF const * parent_vcf;

    //#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
    std::string id; // not used
    std::string qual; // not used
    std::string filter; // not used
    std::string info;
    SampleIndexToSampleInfo sampleIndex_to_sampleInfo; //it is fine to leave this public

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // constructors, destructors, operator=, etc
    // TODO: make sure only consistent VCFs are built (e.g. at least two alleles: ref + 1 alt)?
    VCFRecord(VCF const * parent_vcf, const std::string &chrom, uint32_t pos, const std::string &ref, const std::string &alt,
              const std::string &info=".", const std::string &graph_type_info="");
    VCFRecord(VCF const * parent_vcf);
    VCFRecord(const VCFRecord &) = default;
    virtual std::shared_ptr<VCFRecord> make_copy_as_shared_ptr () const {
        return std::make_shared<VCFRecord>(*this);
    }
    virtual VCFRecord &operator=(const VCFRecord &) = default;
    virtual ~VCFRecord(){}
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // getters
    virtual inline const std::string &get_ref() const {
        return ref;
    }

    virtual inline const std::vector<std::string> &get_alts() const {
        return alts;
    }

    virtual inline uint32_t get_number_of_alleles() const {
        return 1 + alts.size();
    }

    virtual inline bool allele_is_valid (const std::string &allele) const {
        return allele != "" and allele != ".";
    }

    const std::string &get_chrom() const {
        return chrom;
    }

    uint32_t get_pos() const {
        return pos;
    }

    virtual inline uint32_t get_ref_end_pos () const {
        if (allele_is_valid(get_ref())) {
            return pos + get_ref().length();
        } else {
            return pos;
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // modifiers
    virtual inline void set_ref(std::string ref);
    virtual inline void add_new_alt(std::string alt);

    template <class ITERATOR_TYPE>
    inline void add_new_alts(ITERATOR_TYPE begin, ITERATOR_TYPE end) {
        for(;begin<end;++begin) {
            if (allele_is_valid(*begin)) {
                add_new_alt(*begin);
            }
        }
    }

    virtual inline void add_new_samples(uint32_t number_of_samples);

    virtual inline void reset_sample_infos_to_contain_the_given_number_of_samples (uint32_t number_of_samples);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // comparison operators
    virtual bool operator==(const VCFRecord &y) const;
    virtual bool operator!=(const VCFRecord &y) const;
    virtual bool operator<(const VCFRecord &y) const;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // to_string related methods
    virtual std::string alts_to_string() const;
    virtual std::string get_format (bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const;
    virtual std::string sample_infos_to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const {
        return this->sampleIndex_to_sampleInfo.to_string(genotyping_from_maximum_likelihood, genotyping_from_coverage);
    }
    virtual std::string to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // methods querying the INFO
    virtual inline bool graph_type_is_simple() const {
        return this->info.find("GRAPHTYPE=SIMPLE") != std::string::npos;
    }
    virtual inline bool graph_type_is_nested() const {
        return this->info.find("GRAPHTYPE=NESTED") != std::string::npos;
    }
    virtual inline bool graph_type_has_too_many_alts() const {
        return this->info.find("GRAPHTYPE=TOO_MANY_ALTS") != std::string::npos;
    }
    virtual inline bool svtype_is_SNP() const {
        return this->info.find("SVTYPE=SNP") != std::string::npos;
    }
    virtual inline bool svtype_is_indel() const {
        return this->info.find("SVTYPE=INDEL") != std::string::npos;
    }
    virtual inline bool svtype_is_PH_SNPs() const {
        return this->info.find("SVTYPE=PH_SNPs") != std::string::npos;
    }
    virtual inline bool svtype_is_complex() const {
        return this->info.find("SVTYPE=COMPLEX") != std::string::npos;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // general querying
    virtual inline bool ref_allele_is_inside_given_interval(const std::string &chrom, uint32_t pos_from, uint32_t pos_to) const {
        return this->chrom == chrom
               and pos_from <= this->pos
               and this->pos + this->ref.length() <= pos_to;
    }
    virtual inline bool is_SNP () const {
        return ref.length() == 1 and alts.size()==1 and alts[0].length() == 1;
    }
    virtual inline bool has_the_same_position(const VCFRecord &other) const {
        return this->chrom==other.chrom and this->pos==other.pos;
    }
    virtual inline bool has_non_null_reference () const {
        return this->ref != "." and this->ref != "";
    }
    virtual size_t get_longest_allele_length() const;
    virtual bool contains_dot_allele() const;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // genotyping related methods
    virtual inline void genotype_from_coverage() {
        sampleIndex_to_sampleInfo.genotype_from_coverage();
    }

    virtual inline void solve_incompatible_gt_conflict_with (VCFRecord &other) {
        this->sampleIndex_to_sampleInfo.solve_incompatible_gt_conflict_with(other.sampleIndex_to_sampleInfo);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // merging related methods
    virtual inline void merge_record_into_this(const VCFRecord &other) {
        assert(there_are_no_common_alt_alleles_between_this_and_other(other));
        this->sampleIndex_to_sampleInfo.merge_other_samples_infos_into_this(other.sampleIndex_to_sampleInfo);
        add_new_alts(other.get_alts().begin(), other.get_alts().end());
    }

    virtual bool can_biallelic_record_be_merged_into_this (const VCFRecord &vcf_record_to_be_merged_in, uint32_t max_allele_length = 10000) const;



    virtual inline void clear() {
        *this = VCFRecord(parent_vcf);
    }

    virtual inline void correct_dot_alleles_adding_nucleotide_before (char nucleotide) {
        correct_dot_alleles(nucleotide, true);
    }

    virtual inline void correct_dot_alleles_adding_nucleotide_after (char nucleotide) {
        correct_dot_alleles(nucleotide, false);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // TODO: check if we keep this, it is only used in tests - better to keep in a VCFMock class
    // friend std::istream &operator>>(std::istream &in, VCFRecord &m);

protected:
    std::string ref;
    std::vector<std::string> alts;
    std::string chrom;
    uint32_t pos;


    std::string infer_SVTYPE() const;

    virtual void correct_dot_alleles (char nucleotide, bool add_nucleotide_before_the_sequence);

    virtual inline void set_number_of_alleles_and_resize_coverage_information_for_all_samples (uint32_t number_of_alleles) {
        sampleIndex_to_sampleInfo.set_number_of_alleles_and_resize_coverage_information_for_all_samples(number_of_alleles);
    }

    virtual inline bool there_are_no_common_alt_alleles_between_this_and_other(const VCFRecord &other) const;
};

#endif
