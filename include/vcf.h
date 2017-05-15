#ifndef __VCF_H_INCLUDED__   // if vcf.h hasn't been included yet...
#define __VCF_H_INCLUDED__

#include <stdint.h>
#include <ostream>
#include <vector>
#include "vcfrecord.h"

class VCF {
  public:
    std::vector<VCFRecord> records;
    VCF();
    ~VCF();
    void add_record(std::string c, uint32_t p, std::string r, std::string a);
    void add_record(VCFRecord&);
    void clear();
    void save(const std::string&);
    void load(const std::string&);
};

#endif
