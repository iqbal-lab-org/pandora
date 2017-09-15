#ifndef __VCF_H_INCLUDED__   // if vcf.h hasn't been included yet...
#define __VCF_H_INCLUDED__

#include <stdint.h>
#include <ostream>
#include <vector>
#include "vcfrecord.h"

class VCF {
  public:
    std::vector<VCFRecord> records;
    std::vector<std::string> samples;
    VCF();
    ~VCF();
    void add_record(std::string c, uint32_t p, std::string r, std::string a, std::string i=".", std::string g="");
    void add_record(VCFRecord&);
    void add_sample_gt(const std::string& name, const std::string& c, const uint32_t p, const std::string& r, const std::string& a);
    void clear();
    void save(const std::string&, bool simple=false, bool complexgraph=false, bool toomanyalts=false, bool snp=false, bool indel=false, bool phsnps=false, bool complexvar=false);
    void load(const std::string&);
    bool operator == (const VCF& y) const;
};

#endif
