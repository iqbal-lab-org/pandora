#include <iostream>
#include <vector>
#include <zconf.h>

#include <boost/log/trivial.hpp>

#include "inthash.h"
#include "minimizer.h"
#include "seq.h"
#include "utils.h"

using std::vector;

Seq::Seq(uint32_t i, const std::string& n, const std::string& p, uint32_t w, uint32_t k)
    : id(i)
    , name(n)
{
    seq = split_ambiguous(p);
    minimizer_sketch(w, k);
}

Seq::~Seq() { sketch.clear(); }

void Seq::initialize(
    uint32_t i, const std::string& n, const std::string& p, uint32_t w, uint32_t k)
{
    id = i;
    name = n;
    seq = split_ambiguous(p);
    sketch.clear();
    minimizer_sketch(w, k);
}

bool Seq::add_letter_to_get_next_kmer(const char& letter, const uint64_t& shift1,
    const uint64_t& mask, uint32_t& buff, uint64_t (&kmer)[2], uint64_t (&kh)[2])
{
    uint32_t c = nt4((uint8_t)letter);
    if (c < 4) { // not an ambiguous base
        kmer[0] = (kmer[0] << 2 | c) & mask; // forward k-mer
        kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
        kh[0] = hash64(kmer[0], mask);
        kh[1] = hash64(kmer[1], mask);
        buff++;
        return true;
    } else {
        BOOST_LOG_TRIVIAL(debug)
            << now() << "bad letter - found a non AGCT base in read so skipping read "
            << name;
        sketch.clear();
        return false;
    }
}

uint64_t find_smallest_kmer_value(
    const vector<Minimizer>& window, uint& pos_of_smallest)
{
    uint64_t smallest = std::numeric_limits<uint64_t>::max();
    uint i = 0;
    for (const auto& minimizer : window) {
        if (minimizer.canonical_kmer_hash <= smallest) {
            smallest = minimizer.canonical_kmer_hash;
            pos_of_smallest = i;
        }
        i++;
    }
    return smallest;
}

void Seq::add_minimizing_kmers_to_sketch(
    const vector<Minimizer>& window, const uint64_t& smallest)
{
    for (const auto& minimizer : window) {
        if (minimizer.canonical_kmer_hash == smallest) {
            sketch.insert(minimizer);
        }
    }
}

// finds the minimizer in the window, add the minimizer to the sketch set and erase
// everything until the minimizer
void Seq::minimize_window(vector<Minimizer>& window, uint64_t& smallest)
{
    uint pos_of_smallest;
    smallest = find_smallest_kmer_value(window, pos_of_smallest);
    add_minimizing_kmers_to_sketch(window, smallest);
    window.erase(window.begin(), window.begin() + pos_of_smallest + 1);
}

// add the last element of the window (a Minimizer) to the sketch, update the smallest
// and clear the window
void Seq::add_new_smallest_minimizer(vector<Minimizer>& window, uint64_t& smallest)
{
    sketch.insert(window.back());
    smallest = window.back().canonical_kmer_hash;
    window.clear();
}

void Seq::minimizer_sketch(const uint32_t w, const uint32_t k)
{
    // initializations
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1,
        smallest = std::numeric_limits<uint64_t>::max(), kmer[2] = { 0, 0 },
        kh[2] = { 0, 0 };
    uint32_t buff = 0;
    vector<Minimizer> window; // will store all k-mers as Minimizer in the window
    window.reserve(w);

    for (auto &s : seq) {
        const bool sequence_too_short_to_sketch = s.length() + 1 < w + k;
        if (sequence_too_short_to_sketch)
            continue;

        for (const char letter : s) {
            const bool added = add_letter_to_get_next_kmer(letter, shift1, mask, buff,
                kmer,
                kh); // add the next base and remove the first one to get the next kmer
            if (not added)
                return;

            if (buff >= k) {
                window.push_back(Minimizer(
                    std::min(kh[0], kh[1]), buff - k, buff, (kh[0] <= kh[1])));
            }

            if (window.size() == w) {
                minimize_window(window,
                    smallest); // finds the minimizer in the window, add the minimizer to the sketch set and erase everything until the minimizer
            } else if (buff >= w + k
                and window.back().canonical_kmer_hash <= smallest) {
                add_new_smallest_minimizer(window,
                    smallest); // add the last element of the window (a Minimizer) to the sketch, update the smallest and clear the window
            }

            const bool window_has_shortened = window.size() < w;
            if (!window_has_shortened) {
                fatal_error(
                    "Error when sketching sequence: a minimizer should have been added "
                    "and windows should have size < ",
                    w, " (is ", window.size(), ")");
            }
        }
    }
}

std::ostream& operator<<(std::ostream& out, Seq const& data)
{
    out << data.name;
    return out;
}

uint64_t Seq::length() const
{
    uint64_t l{0};
    for (auto &s: seq) {
        l += s.length();
    }
    return l;
}
