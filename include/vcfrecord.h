#ifndef __VCFRECORD_H_INCLUDED__   // if vcfrecord.h hasn't been included yet...
#define __VCFRECORD_H_INCLUDED__

#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include "sampleinfo.h"

struct VCFRecord {
private:
    std::string infer_SVTYPE() const;


public:
    //#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
    std::string chrom;
    uint32_t pos;
    std::string id; // not used
    std::string ref;
    std::vector<std::string> alts;
    std::string qual; // not used
    std::string filter; // not used
    std::string info;

    SamplesInfos sampleIndex_to_sampleInfo;

    // TODO: make sure only consistent VCFs are built (e.g. at least two alleles: ref + 1 alt)?
    VCFRecord(const std::string &chrom, uint32_t pos, const std::string &ref, const std::string &alt,
              const std::string &info=".", const std::string &graph_type_info="");
    VCFRecord();
    VCFRecord(const VCFRecord &) = default;
    VCFRecord &operator=(const VCFRecord &) = default;
    ~VCFRecord() = default;

    std::string get_format (bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const;

    inline bool ref_allele_is_inside_given_interval(const std::string &chrom, uint32_t pos_from, uint32_t pos_to) const {
        return this->chrom == chrom
               and pos_from <= this->pos
               and this->pos + this->ref.length() <= pos_to;
    }

    inline bool is_SNP () const {
        return ref.length() == 1 and alts.size()==1 and alts[0].length() == 1;

    }


    inline void clear() {
        *this = VCFRecord();
    }


    inline void genotype() {
        sampleIndex_to_sampleInfo.genotype_from_coverage();
    }

    bool contains_dot_allele() const;

    bool operator==(const VCFRecord &y) const;

    bool operator!=(const VCFRecord &y) const;

    bool operator<(const VCFRecord &y) const;

    std::string to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const;

    // TODO: check if we keep this, it is only used in tests - better to keep in a VCFMock class
    // friend std::istream &operator>>(std::istream &in, VCFRecord &m);


    inline bool graph_type_is_simple() const {
        return this->info.find("GRAPHTYPE=SIMPLE") != std::string::npos;
    }
    inline bool graph_type_is_nested() const {
        return this->info.find("GRAPHTYPE=NESTED") != std::string::npos;
    }
    inline bool graph_type_has_too_many_alts() const {
        return this->info.find("GRAPHTYPE=TOO_MANY_ALTS") != std::string::npos;
    }
    inline bool svtype_is_SNP() const {
        return this->info.find("SVTYPE=SNP") != std::string::npos;
    }
    inline bool svtype_is_indel() const {
        return this->info.find("SVTYPE=INDEL") != std::string::npos;
    }
    inline bool svtype_is_PH_SNPs() const {
        return this->info.find("SVTYPE=PH_SNPs") != std::string::npos;
    }
    inline bool svtype_is_complex() const {
        return this->info.find("SVTYPE=COMPLEX") != std::string::npos;
    }
    inline bool has_the_same_position(const VCFRecord &other) const {
        return this->chrom==other.chrom and this->pos==other.pos;
    }
    inline bool has_non_null_reference () const {
        return this->ref != "." and this->ref != "";
    }

    size_t get_longest_allele_length() const;

    // MERGING-RELATED METHODS
    inline void merge_record_into_this(const VCFRecord &other) {
        this->sampleIndex_to_sampleInfo.merge_other_samples_infos_into_this(other.sampleIndex_to_sampleInfo, alts.size());
        add_alts(other);
    }

private:
    inline void add_alts(const VCFRecord &other) {
        this->alts.insert(this->alts.end(), other.alts.begin(), other.alts.end());
    }
    // END OF MERGING-RELATED METHODS
};

#endif
