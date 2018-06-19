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

    void add_record(VCFRecord &);

    ptrdiff_t get_sample_index(const std::string&);

    void add_sample_gt(const std::string &name, const std::string &c, const uint32_t p, const std::string &r,
                       const std::string &a);

    void add_sample_ref_alleles(const std::string &, const std::string &, const uint32_t &, const uint32_t &);

    void clear();

    void sort_records();

    bool pos_in_range(const uint32_t , const uint32_t);

    void regenotype(const uint32_t&, const float&,const uint8_t);

    void save(const std::string &, bool simple = false, bool complexgraph = false, bool toomanyalts = false,
              bool snp = false, bool indel = false, bool phsnps = false, bool complexvar = false);

    void load(const std::string &);

    void write_aligned_fasta(const std::string &, const std::vector<LocalNodePtr> &);

    bool operator==(const VCF &y) const;
};

#endif
