#ifndef __MINIHIT_H_INCLUDED__   // if minihit.h hasn't been included yet...
#define __MINIHIT_H_INCLUDED__

#include <functional>
#include <stdint.h>
#include <ostream>
#include "minimizer.h"
#include "minirecord.h"

using namespace std;

struct MinimizerHit
{
    uint32_t read_id;
    Interval read_interval;
    uint32_t prg_id;
    Path prg_path;
    bool strand; // forward or reverse complement

    MinimizerHit(const uint32_t i, const Minimizer* m, const MiniRecord* r);
    MinimizerHit(const uint32_t i, const Interval j, const uint32_t k, const Path p, const bool c); // second allowed constructor
    bool operator < ( const MinimizerHit& y) const;
    bool operator == (const MinimizerHit& y) const;
    friend ostream& operator<< (ostream& out, const MinimizerHit& m);
};

#endif
