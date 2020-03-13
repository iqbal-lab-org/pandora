#ifndef __SEQ_H_INCLUDED__ // if seq.h hasn't been included yet...
#define __SEQ_H_INCLUDED__

#include <string>
#include <cstdint>
#include <set>
#include <ostream>
#include "minimizer.h"

class Seq {
public:
    uint32_t id;
    std::string name;
    std::string seq;
    std::set<Minimizer> sketch;

    Seq(uint32_t, const std::string&, const std::string&, uint32_t, uint32_t);

    ~Seq();

    void initialize(
        uint32_t, const std::string&, const std::string&, uint32_t, uint32_t);

    bool add_letter_to_get_next_kmer(const char&, const uint64_t&, const uint64_t&,
        uint32_t&, uint64_t (&)[2], uint64_t (&)[2]);

    void add_minimizing_kmers_to_sketch(const std::vector<Minimizer>&, const uint64_t&);

    void minimize_window(std::vector<Minimizer>&, uint64_t&);

    void add_new_smallest_minimizer(std::vector<Minimizer>&, uint64_t&);

    void minimizer_sketch(const uint32_t w, const uint32_t k);

    friend std::ostream& operator<<(std::ostream& out, const Seq& data);
};

#endif
