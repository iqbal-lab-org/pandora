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


class LocalNode;

typedef std::shared_ptr<LocalNode> LocalNodePtr;

class VCF {
private:
    /* will contain, for each chromosome, an interval tree containing VCF records interval and a pointer to the VCF Record itself to allow
       VCF::make_gt_compatible() to execute a lot faster than serial search */
    std::map<std::string, IITree<uint32_t, VCFRecord*>> chrom2recordIntervalTree;

    //add a VCF record to this VCF
    void add_record_core(const VCFRecord &vr);

    //find a VCRRecord in records
    std::vector<std::shared_ptr<VCFRecord>>::iterator find_record_in_records(const VCFRecord &vr) {
        return find_if(records.begin(), records.end(), [&vr](const std::shared_ptr<VCFRecord> &record) { return *record==vr; });
    }
    std::vector<std::shared_ptr<VCFRecord>>::const_iterator find_record_in_records(const VCFRecord &vr) const {
        return find_if(records.begin(), records.end(), [&vr](const std::shared_ptr<VCFRecord> &record) { return *record==vr; });
    }


public:
    std::vector<std::shared_ptr<VCFRecord>> records;
    std::vector<std::string> samples;

    //constructor/destructors
    VCF() = default;
    virtual ~VCF() = default;

    inline size_t get_VCF_size() const {
        return records.size();
    }

    void add_record(std::string c, uint32_t p, std::string r, std::string a, std::string i = ".", std::string g = "");

    VCFRecord &add_record(VCFRecord &, const std::vector<std::string> &sample_names);

    void add_samples(const std::vector<std::string>);

    void add_formats(const std::vector<std::string> &);

    ptrdiff_t get_sample_index(const std::string &);

    void add_sample_gt(const std::string &name, const std::string &c, const uint32_t p, const std::string &r,
                       const std::string &a);

    void add_sample_ref_alleles(const std::string &, const std::string &, const uint32_t &, const uint32_t &);

    void append_vcf(const VCF &);

    void sort_records();

    bool pos_in_range(const uint32_t, const uint32_t, const std::string &) const;

    void genotype(const std::vector<uint32_t> &, const float &, const uint16_t, const uint32_t &min_allele_covg,
                  const float &min_fraction_allele_covg, const uint32_t &min_site_total_covg,
                  const uint32_t &min_site_diff_covg, bool snps_only);

    void clean();

    VCF merge_multi_allelic(uint32_t max_allele_length = 10000) const;

    void correct_dot_alleles(const std::string &, const std::string &);

    void make_gt_compatible();

    bool operator==(const VCF &y) const;

    bool operator!=(const VCF &y) const;

    /**
     * Concatenate several VCF files that were previously written to disk as .vcfs into a single VCF file
     * @param VCFPathsToBeConcatenated : vector containing paths to the .vcfs to be concatenated
     * @param sink : where to put the concatenated VCFs
     */
    static void concatenateVCFs(const std::vector<std::string> &VCFPathsToBeConcatenated, const std::string &sink);



    // serialization operations
    void save(const std::string &filepath, bool output_dot_allele = false, bool graph_is_simple = true, bool graph_is_nested = true, bool graph_has_too_many_alts = true, bool sv_type_is_snp = true, bool sv_type_is_indel = true,
              bool sv_type_is_ph_snps = true, bool sv_type_is_complex = true);
    virtual std::string header() const;
    std::string to_string(bool output_dot_allele = false, bool graph_is_simple = true, bool graph_is_nested = true, bool graph_has_too_many_alts = true, bool sv_type_is_snp = true, bool sv_type_is_indel = true,
                          bool sv_type_is_ph_snps = true, bool sv_type_is_complex = true);
    void load(const std::string &filepath);
};

#endif
