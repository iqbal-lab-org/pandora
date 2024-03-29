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
    std::set<Minimizer> sketch;

    // the original sequence
    std::string full_seq;

    // the original sequence is split into several valid subsequences (composed of ACGT only)
    // TODO: now that we are storing the original sequence, these valid subsequences
    // TODO: could be more efficiently stored with std::string_view instead,
    // TODO: although I don't think the RAM improvement will be significant, so
    // TODO: just adding this as a future note.
    std::vector<std::string> subseqs;  // these are the subsequences themselves
    std::vector<size_t> offsets;  // these are the subsequences offsets on the original string


    Seq(uint32_t, const std::string&, const std::string&, uint32_t, uint32_t);

    ~Seq();

    void initialize(
        uint32_t, const std::string&, const std::string&, uint32_t, uint32_t);

    void add_letter_to_get_next_kmer(const char&, const uint64_t&, const uint64_t&,
        uint32_t&, uint64_t (&)[2], uint64_t (&)[2]);

    void add_minimizing_kmers_to_sketch(const std::vector<Minimizer>&, const uint64_t&);

    void minimize_window(std::vector<Minimizer>&, uint64_t&);

    void add_new_smallest_minimizer(std::vector<Minimizer>&, uint64_t&);

    uint64_t length() const;

    friend std::ostream& operator<<(std::ostream& out, const Seq& data);

    void minimizer_sketch(const uint32_t w, const uint32_t k);

private:
    void minimizer_sketch(const std::string &s, const size_t seq_offset,
        const uint32_t w, const uint32_t k);

};

#endif
