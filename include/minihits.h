#ifndef __MINIHITS_H_INCLUDED__   // if minihits.h hasn't been included yet...
#define __MINIHITS_H_INCLUDED__

struct MinimizerHit;

#include <set>
#include <stdint.h>
#include "minimizer.h"
#include "minirecord.h"

using namespace std;

struct pComp
{
  bool operator()(MinimizerHit* lhs, MinimizerHit* rhs);
};

struct clusterComp
{
  bool operator()(set<MinimizerHit*> lhs, set<MinimizerHit*> rhs);
};

class MinimizerHits {
  public:
    MinimizerHits() {}
    ~MinimizerHits();
    set<MinimizerHit*, pComp> hits;
    void add_hit(const uint32_t i, const Minimizer* m, const MiniRecord r, const uint8_t c);
    friend ostream& operator<< (ostream& out, const MinimizerHits& m);
};

#endif
