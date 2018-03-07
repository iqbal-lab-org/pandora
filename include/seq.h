#ifndef __SEQ_H_INCLUDED__   // if seq.h hasn't been included yet...
#define __SEQ_H_INCLUDED__

#include <string>
#include <set>
#include <ostream>
#include "minimizer.h"

class Seq {
public:
    uint32_t id;
    std::string name;
    std::string seq;
    std::set<Minimizer *, pMiniComp> sketch;

    Seq(uint32_t, std::string, std::string, uint32_t, uint32_t);

    ~Seq();

    void initialize(uint32_t, std::string, std::string, uint32_t, uint32_t);

    void minimizer_sketch(const uint32_t w, const uint32_t k);

    friend std::ostream &operator<<(std::ostream &out, const Seq &data);
};

#endif
