#ifndef __SEQ_H_INCLUDED__ // if seq.h hasn't been included yet...
#define __SEQ_H_INCLUDED__

#include <string>
#include <cstdint>
#include <set>
#include <ostream>
#include "minimizer.h"
#include "inthash.h"

class Seq {
public:
    uint32_t id;
    std::string name;
    std::set<Minimizer> sketch;

    // the original sequence is split into several valid subsequences (composed of ACGT only)
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

    /// Take a substring of the Seq object. This extracts the corresponding substrings
    /// from the `subseqs` based on the `offsets` positions.
    /// If the requested substring extends past the end of the string, i.e. the count
    /// is greater than size() - pos (e.g. if count == npos), the returned substring is [pos, size()).
    /// Where the interval to be extract overlaps a "gap" (i.e., where the original
    /// read had an ambiguous base(s)), TODO: what do we do in this?
    /// \param pos position of the first character to include
    /// \param count length of the substring
    /// \return String containing the substring [pos, pos+count) or [pos, size()).
    std::string substr(size_t pos, size_t count = std::string::npos) const;

private:
    void minimizer_sketch(const std::string &s, const size_t seq_offset,
        const uint32_t w, const uint32_t k);

};

#endif
