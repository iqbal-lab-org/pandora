#ifndef __VCF_H_INCLUDED__ // if vcf.h hasn't been included yet...
#define __VCF_H_INCLUDED__

#include <ostream>
#include <vector>
#include <memory>
#include <string>
#include <cstdint>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <map>

#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include "vcfrecord.h"
#include "IITree.h"
#include "OptionsAggregator.h"

#include "vcfrecord.h"
#include "utils.h"

class LocalNode;
class VCFRecord;

namespace fs = boost::filesystem;

typedef std::shared_ptr<LocalNode> LocalNodePtr;

class VCF {
public:
    // TODO : protect this member
    GenotypingOptions const* genotyping_options;

    // TODO : protect this member
    std::vector<std::string> samples;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // constructor/destructors
    // TODO: make VCF not sample updatable - we have to specify the samples upfront -
    // this can be done in pandora and
    //  will render the code less complicated and less error-prone
    VCF(GenotypingOptions const* genotyping_options)
        : genotyping_options(genotyping_options)
    {
    }
    virtual ~VCF() { }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // comparison operators
    virtual bool operator==(const VCF& y) const;
    virtual bool operator!=(const VCF& y) const;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // static constants for common VCF fields
    static const std::string VARIANT_CLASS_ID;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // adders
    virtual void add_record(const std::string& chrom, uint32_t position,
        const std::string& ref, const std::string& alt, const std::string& info = ".",
        const std::string& graph_type_info = "");
    virtual void add_record(const VCFRecord& vcf_record);
    virtual VCFRecord& add_or_update_record_restricted_to_the_given_samples(
        VCFRecord& vr, const std::vector<std::string>& sample_names);
    virtual void add_samples(const std::vector<std::string>& sample_names);
    virtual void append_vcf(const VCF&);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // getters
    const std::vector<std::shared_ptr<VCFRecord>>& get_records() const
    {
        return records;
    }

    virtual inline size_t get_VCF_size() const { return records.size(); }
    virtual ptrdiff_t get_sample_index(const std::string&);
    virtual std::vector<VCFRecord*> get_all_records_overlapping_the_given_record(
        const VCFRecord& vcf_record);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // genotyping related methods
    virtual void add_a_new_record_discovered_in_a_sample_and_genotype_it(
        const std::string& sample_name, const std::string& chrom, const uint32_t pos,
        const std::string& ref, const std::string& alt);
    virtual void set_sample_gt_to_ref_allele_for_records_in_the_interval(
        const std::string& sample_name, const std::string& chrom,
        const uint32_t& pos_from, const uint32_t& pos_to);
    virtual void genotype(const bool do_local_genotyping);
    virtual void make_gt_compatible();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // merging related methods
    virtual inline VCF merge_multi_allelic(uint32_t max_allele_length = 100000) const
    {
        VCF merged_VCF(this->genotyping_options);
        merge_multi_allelic_core(merged_VCF, max_allele_length);
        return merged_VCF;
    }
    virtual VCF correct_dot_alleles(
        const std::string& vcf_ref, const std::string& chrom) const;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // to_string methods
    virtual std::string header() const;
    virtual std::string to_string(bool genotyping_from_maximum_likelihood,
                                  bool genotyping_from_coverage, bool output_dot_allele = false,
                                  bool graph_is_simple = true, bool graph_is_nested = true,
                                  bool graph_has_too_many_alts = true, bool variant_class_is_snp = true,
                                  bool variant_class_is_indel = true, bool variant_class_is_ph_snps = true,
                                  bool variant_class_is_complex = true);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // misc methods
    virtual void
    sort_records(); // TODO: remove this method and store the records always sorted
    virtual bool pos_in_range(const uint32_t, const uint32_t, const std::string&) const;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // serialization operations
    // TODO : give a file handler instead of filepath to allow for mocking
    // TODO : the file handler must be a wrapper on ostream, to allow for mocking in
    // fact
    virtual void save(const fs::path& filepath, bool genotyping_from_maximum_likelihood,
        bool genotyping_from_coverage, bool output_dot_allele = false,
        bool graph_is_simple = true, bool graph_is_nested = true,
        bool graph_has_too_many_alts = true, bool variant_class_is_snp = true,
        bool variant_class_is_indel = true, bool variant_class_is_ph_snps = true,
        bool variant_class_is_complex = true);

    // concatenate several VCF files that were previously written to disk as .vcfs into
    // a single VCF file
    static void concatenate_VCFs(
        const std::vector<fs::path>& VCF_paths_to_be_concatenated,
        const fs::path& final_VCF_file);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

protected:
    std::vector<std::shared_ptr<VCFRecord>> records;
    /* will contain, for each chromosome, an interval tree containing VCF records
       interval and a pointer to the VCF Record itself to allow
       VCF::make_gt_compatible() to execute a lot faster than serial search */
    std::map<std::string, IITree<uint32_t, VCFRecord*>> chrom_to_record_interval_tree;

    // add a VCF record to this VCF
    virtual void add_record_core(const VCFRecord& vr);

    // find a VCRRecord in records
    virtual inline std::vector<std::shared_ptr<VCFRecord>>::iterator
    find_record_in_records(const VCFRecord& vr);
    virtual inline std::vector<std::shared_ptr<VCFRecord>>::const_iterator
    find_record_in_records(const VCFRecord& vr) const;

    virtual void update_other_samples_of_this_record(VCFRecord* reference_record);

    virtual void merge_multi_allelic_core(
        VCF& merged_vcf, uint32_t max_allele_length) const;

    virtual inline std::string get_current_date() const;
};

#endif
