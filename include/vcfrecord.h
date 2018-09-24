#ifndef __VCFRECORD_H_INCLUDED__   // if vcfrecord.h hasn't been included yet...
#define __VCFRECORD_H_INCLUDED__

#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>

struct VCFRecord {
    //#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
    std::string chrom;
    uint32_t pos;
    std::string id; // not used
    std::string ref;
    std::vector<std::string> alt;
    std::string qual; // not used
    std::string filter; // not used
    std::string info;
    std::vector<std::string> format; //e.g. "GT"
    std::vector<std::unordered_map<std::string, std::vector<uint8_t>>> samples;      // should have an entry for each sample in vcf,
    std::vector<std::unordered_map<std::string, std::vector<float>>> regt_samples;   // in the same order

    VCFRecord(std::string, uint32_t, std::string, std::string, std::string i = ".", std::string g = "");

    VCFRecord();

    VCFRecord(const VCFRecord &);

    VCFRecord &operator=(const VCFRecord &);

    ~VCFRecord();

    void clear();

    void clear_sample(uint32_t);

    void add_formats(const std::vector<std::string>&);

    void likelihood(const uint32_t&, const float&);

    void confidence();

    void genotype(const uint8_t);

    bool operator==(const VCFRecord &y) const;

    bool operator!=(const VCFRecord &y) const;

    bool operator<(const VCFRecord &y) const;

    friend std::ostream &operator<<(std::ostream &out, const VCFRecord &m);

    friend std::istream &operator>>(std::istream &in, VCFRecord &m);
};

#endif
