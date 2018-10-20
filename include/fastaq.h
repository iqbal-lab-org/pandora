#ifndef __FASTAQ_H_INCLUDED__   // if fastaq.h hasn't been included yet...
#define __FASTAQ_H_INCLUDED__

#include <vector>
#include <unordered_map>
#include <iostream>
#include <cmath>


struct Fastaq {
    bool gzipped;
    bool fastq;
    std::vector<std::string> names;
    std::unordered_map<std::string, std::string> headers;
    std::unordered_map<std::string, std::string> sequences;
    std::unordered_map<std::string, std::string> scores;

    Fastaq(bool gz = false, bool fq = false);

    char covg_to_score(const uint_least16_t &, const uint_least16_t &, const bool &alt = false);

    char alt_covg_to_score(const uint_least16_t &covg);

    void add_entry(const std::string &, const std::string &, const std::vector<uint32_t> &,
                   const uint_least16_t, const std::string header = "");

    void add_entry(const std::string &, const std::string &, const std::string header = "");

    void clear();

    void save(const std::string &);

    bool operator==(const Fastaq &y) const;

    bool operator!=(const Fastaq &y) const;

    friend std::ostream &operator<<(std::ostream &out, const Fastaq &m);

    friend std::istream &operator>>(std::istream &in, Fastaq &m);
};

#endif
