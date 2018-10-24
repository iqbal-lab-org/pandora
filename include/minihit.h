#ifndef __MINIHIT_H_INCLUDED__   // if minihit.h hasn't been included yet...
#define __MINIHIT_H_INCLUDED__

#include <ostream>
#include <cstdint>
#include "minimizer.h"
#include "minirecord.h"


struct MinimizerHit {
    uint32_t read_id;
    uint32_t read_start_position;
    uint32_t prg_id;
    prg::Path prg_path;
    uint32_t knode_id;
    bool strand; // forward or reverse complement

    MinimizerHit(const uint32_t i, const Minimizer &m, const MiniRecord *r);

    MinimizerHit(const uint32_t i, const Interval j, const uint32_t k, const prg::Path p, const uint32_t n,
                 const bool c); // second allowed constructor
    bool operator<(const MinimizerHit &y) const;

    bool operator==(const MinimizerHit &y) const;

    friend std::ostream &operator<<(std::ostream &out, const MinimizerHit &m);
};

#endif
