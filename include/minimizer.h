#ifndef __MINIMIZER_H_INCLUDED__   // if minimizer.h hasn't been included yet...
#define __MINIMIZER_H_INCLUDED__

#include <ostream>
#include <cstdint>
#include "interval.h"

//Represent a sequence minimizer of a READ!
struct Minimizer {
    uint64_t kmer; //this is the minimum canonical kmer hashed value in fact
    Interval pos; //position of the kmer in the read
    bool strand; //strand of the kmer

    Minimizer() {};

    Minimizer(uint64_t, uint32_t, uint32_t, bool);

    ~Minimizer();

    bool operator<(const Minimizer &y) const;

    bool operator==(const Minimizer &y) const;

    friend std::ostream &operator<<(std::ostream &out, const Minimizer &m);
};

#endif
