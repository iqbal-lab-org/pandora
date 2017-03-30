#ifndef __MINIHITS_H_INCLUDED__   // if minihits.h hasn't been included yet...
#define __MINIHITS_H_INCLUDED__

struct MinimizerHit;

#include <set>
#include <unordered_set>
#include <stdint.h>
#include "minimizer.h"
#include "minirecord.h"

struct pComp
{
  bool operator()(MinimizerHit* lhs, MinimizerHit* rhs);
};

struct pComp_path
{
  bool operator()(MinimizerHit* lhs, MinimizerHit* rhs);
};

struct clusterComp
{
  bool operator()(std::set<MinimizerHit*, pComp> lhs, std::set<MinimizerHit*, pComp> rhs);
};

class MinimizerHits {
  public:
    MinimizerHits();
    ~MinimizerHits();
    std::unordered_set<MinimizerHit*> uhits;
    std::set<MinimizerHit*, pComp> hits;
    void add_hit(const uint32_t i, const Minimizer* m, const MiniRecord* r);
    void sort();
    friend std::ostream& operator<< (std::ostream& out, const MinimizerHits& m);
};

#endif
