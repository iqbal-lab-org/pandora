#ifndef __VCF_H_INCLUDED__   // if vcf.h hasn't been included yet...
#define __VCF_H_INCLUDED__

#include <ostream>
#include <vector>
#include <memory>
#include <string>
#include <cstdint>
#include "vcfrecord.h"
#include "IITree.h"
#include <map>
#include "OptionsAggregator.h"


class LocalNode;

typedef std::shared_ptr<LocalNode> LocalNodePtr;

class VCF {
public:
    std::vector<std::shared_ptr<VCFRecord>> records;
    GenotypingOptions const * genotyping_options;
    std::vector<std::string> samples;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //constructor/destructors
    // TODO: make VCF not sample updatable - we have to specify the samples upfront - this can be done in pandora
    // and it makes a lot of things easier and less error prone
    VCF(GenotypingOptions const * genotyping_options) : genotyping_options(genotyping_options){}
    virtual ~VCF() = default;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // comparison operators
    virtual bool operator==(const VCF &y) const;
    virtual bool operator!=(const VCF &y) const;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // adders
    virtual void add_record(std::string c, uint32_t p, std::string r, std::string a, std::string i = ".", std::string g = "");
    virtual VCFRecord &add_record(VCFRecord &, const std::vector<std::string> &sample_names);
    virtual void add_samples(const std::vector<std::string>);
    virtual void append_vcf(const VCF &);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // getters
    virtual inline size_t get_VCF_size() const {
        return records.size();
    }
    virtual ptrdiff_t get_sample_index(const std::string &);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // genotyping related methods
    virtual void add_a_new_record_discovered_in_a_sample_and_genotype_it(const std::string &sample_name, const std::string &chrom, const uint32_t pos, const std::string &ref,
                                                                         const std::string &alt);
    virtual void set_sample_gt_to_ref_allele_for_records_in_the_interval(const std::string &sample_name, const std::string &chrom, const uint32_t &pos_from, const uint32_t &pos_to);
    virtual void genotype();
    virtual void make_gt_compatible();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // merging related methods
    virtual VCF merge_multi_allelic(uint32_t max_allele_length = 10000) const;
    virtual void correct_dot_alleles(const std::string &, const std::string &);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // to_string methods
    virtual std::string header() const;
    virtual std::string to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage,
                                  bool output_dot_allele = false, bool graph_is_simple = true, bool graph_is_nested = true, bool graph_has_too_many_alts = true, bool sv_type_is_snp = true, bool sv_type_is_indel = true,
                                  bool sv_type_is_ph_snps = true, bool sv_type_is_complex = true);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // misc methods
    virtual void sort_records();
    virtual void clean();
    virtual bool pos_in_range(const uint32_t, const uint32_t, const std::string &) const;

    /**
     * Concatenate several VCF files that were previously written to disk as .vcfs into a single VCF file
     * @param VCFPathsToBeConcatenated : vector containing paths to the .vcfs to be concatenated
     * @param sink : where to put the concatenated VCFs
     */
    static void concatenateVCFs(const std::vector<std::string> &VCFPathsToBeConcatenated, const std::string &sink);

    // serialization operations
    virtual void save(const std::string &filepath, bool output_dot_allele = false, bool graph_is_simple = true, bool graph_is_nested = true, bool graph_has_too_many_alts = true, bool sv_type_is_snp = true, bool sv_type_is_indel = true,
                      bool sv_type_is_ph_snps = true, bool sv_type_is_complex = true);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // TODO: check if we keep this, it is only used in tests - better to keep in a VCFMock class
    // void load(const std::string &filepath);
protected:
    /* will contain, for each chromosome, an interval tree containing VCF records interval and a pointer to the VCF Record itself to allow
       VCF::make_gt_compatible() to execute a lot faster than serial search */
    std::map<std::string, IITree<uint32_t, VCFRecord*>> chrom_to_record_interval_tree;

    //add a VCF record to this VCF
    virtual void add_record_core(const VCFRecord &vr);

    //find a VCRRecord in records
    virtual inline std::vector<std::shared_ptr<VCFRecord>>::iterator find_record_in_records(const VCFRecord &vr) {
        return find_if(records.begin(), records.end(), [&vr](const std::shared_ptr<VCFRecord> &record) { return *record==vr; });
    }
    virtual inline std::vector<std::shared_ptr<VCFRecord>>::const_iterator find_record_in_records(const VCFRecord &vr) const {
        return find_if(records.begin(), records.end(), [&vr](const std::shared_ptr<VCFRecord> &record) { return *record==vr; });
    }

    virtual std::vector<VCFRecord*> get_all_records_overlapping_the_given_record (const VCFRecord &vcf_record);

    virtual void update_other_samples_of_this_record(VCFRecord *reference_record);

};

#endif
