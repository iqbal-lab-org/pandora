#ifndef __VCFRECORD_H_INCLUDED__   // if vcfrecord.h hasn't been included yet...
#define __VCFRECORD_H_INCLUDED__

#include <iostream>
#include <algorithm>

struct VCFRecord
{
    //#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
    std::string chrom;
    uint32_t pos;
    std::string id;
    std::string ref;
    std::string alt;
    std::string qual;
    std::string filter;
    std::string info;
    std::vector<std::string> samples;

    VCFRecord(std::string, uint32_t, std::string, std::string, std::string i = ".", std::string g = "");
    VCFRecord();
    ~VCFRecord();
    bool operator == (const VCFRecord& y) const;
    bool operator <  (const VCFRecord& y) const;
    friend std::ostream& operator<< (std::ostream& out, const VCFRecord& m);
    friend std::istream& operator>> (std::istream& in, VCFRecord& m);
};

#endif
