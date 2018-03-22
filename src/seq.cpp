#include <iostream>
#include <vector>
#include <cassert>
#include "inthash.h"
#include "minimizer.h"
#include "seq.h"
#include "utils.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using std::vector;
using namespace std;

Seq::Seq(uint32_t i, string n, string p, uint32_t w, uint32_t k) : id(i), name(n), seq(p) {
    minimizer_sketch(w, k);
}

Seq::~Seq() {
    sketch.clear();
}

void Seq::initialize(uint32_t i, string n, string p, uint32_t w, uint32_t k) {
    id = i;
    name = n;
    seq = p;
    sketch.clear();
    minimizer_sketch(w, k);
}

bool Seq::add_letter_to_get_next_kmer(const char& letter,
                                      const uint64_t& shift1,
                                      const uint64_t& mask,
                                      uint32_t&  buff,
                                      uint64_t (&kmer)[2],
                                      uint64_t (&kh)[2]){
    uint c = nt4((uint8_t) letter);
    if (c < 4) { // not an ambiguous base
        kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
        kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
        kh[0] = hash64(kmer[0], mask);
        kh[1] = hash64(kmer[1], mask);
        buff++;
        return true;
    } else {
        cout << now() << "bad letter - found a non AGCT base in read so skipping read " << name << endl;
        sketch.clear();
        return false;
    }
}

uint64_t find_smallest_kmer_value(const vector<Minimizer>& window, uint & pos_of_smallest){
    uint64_t smallest = std::numeric_limits<uint64_t>::max();
    uint i = 0;
    for (const auto minimizer : window) {
        if (minimizer.kmer <= smallest) {
            smallest = minimizer.kmer;
            pos_of_smallest = i;
        }
        i++;
    }
    return smallest;
}

void Seq::add_minimizing_kmers_to_sketch(const vector<Minimizer>& window, const uint64_t& smallest){
    for (const auto minimizer : window) {
        if (minimizer.kmer == smallest) {
            sketch.insert(minimizer);
            //num_minis_found += 1;
        }
    }
}

void Seq::minimize_window(vector<Minimizer>& window, uint64_t& smallest) {
    uint pos_of_smallest;
    smallest = find_smallest_kmer_value(window, pos_of_smallest);
    add_minimizing_kmers_to_sketch(window, smallest);
    window.erase(window.begin(), window.begin() + pos_of_smallest + 1);
}

void Seq::add_new_smallest_minimizer(vector<Minimizer>& window, uint64_t& smallest) {
    sketch.insert(window.back());
    smallest = window.back().kmer;
    window.clear();
}

void Seq::minimizer_sketch(const uint32_t w, const uint32_t k) {
    bool sequence_too_short_to_sketch = seq.length() + 1 < w + k;
    if (sequence_too_short_to_sketch)
        return;

    // initializations
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1,
            smallest = std::numeric_limits<uint64_t>::max(),
            kmer[2] = {0, 0}, kh[2] = {0, 0};
    uint32_t buff = 0;
    vector<Minimizer> window;
    window.reserve(w);

    for (const char letter : seq) {
        bool added = add_letter_to_get_next_kmer(letter, shift1, mask, buff, kmer, kh);
        if (not added)
            return;

        if (buff >= k) {
            window.push_back(Minimizer(min(kh[0], kh[1]), buff - k, buff, (kh[0] <= kh[1])));
        }

        if (window.size() == w) {
            minimize_window(window, smallest);
        } else if (buff >= w + k and window.back().kmer <= smallest) {
            add_new_smallest_minimizer(window,smallest);
        }
        assert(window.size() < w ||
               assert_msg("we can't have added a smallest kmer correctly as window still has size " << window.size()));
    }
    //cout << now() << "Sketch size " << sketch.size() << " for read " << name << endl;
    return;
}

std::ostream &operator<<(std::ostream &out, Seq const &data) {
    out << data.name;
    return out;
}

