#ifndef __VCF_H_INCLUDED__   // if vcf.h hasn't been included yet...
#define __VCF_H_INCLUDED__

#include <ostream>
#include <vector>
#include <memory>
#include <string>
#include <cstdint>
#include "vcfrecord.h"


class LocalNode;

typedef std::shared_ptr<LocalNode> LocalNodePtr;

class VCF {
public:
    std::vector<VCFRecord> records;
    std::vector<std::string> samples;

    VCF();

    ~VCF();

    void add_record(std::string c, uint32_t p, std::string r, std::string a, std::string i = ".", std::string g = "");

    VCFRecord &add_record(VCFRecord &, const std::vector<std::string> &sample_names);

    void add_samples(const std::vector<std::string>);

    void add_formats(const std::vector<std::string> &);

    ptrdiff_t get_sample_index(const std::string &);

    void add_sample_gt(const std::string &name, const std::string &c, const uint32_t p, const std::string &r,
                       const std::string &a);

    void add_sample_ref_alleles(const std::string &, const std::string &, const uint32_t &, const uint32_t &);

    void clear();

    void append_vcf(const VCF &);

    void sort_records();

    bool pos_in_range(const uint32_t, const uint32_t, const std::string &) const;

    void genotype(const uint32_t &, const float &, const uint16_t, const uint32_t &min_allele_covg,
                      const float &min_fraction_allele_covg, const uint32_t &min_site_total_covg,
                      const uint32_t &min_site_diff_covg, bool snps_only);

    void clean();

    void merge_multi_allelic(uint32_t max_allele_length = 10000);

    void correct_dot_alleles(const std::string &, const std::string &);

    void make_gt_compatible();

    std::string header();

    void save(const std::string &, bool simple = false, bool complexgraph = false, bool toomanyalts = false,
              bool snp = false, bool indel = false, bool phsnps = false, bool complexvar = false);

    void load(const std::string &);

    bool operator==(const VCF &y) const;

    bool operator!=(const VCF &y) const;

    friend std::ostream &operator<<(std::ostream &out, const VCF &m);
};

#endif
