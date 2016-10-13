#ifndef __MINIMIZER_H_INCLUDED__   // if minimizer.h hasn't been included yet...
#define __MINIMIZER_H_INCLUDED__

#include <string>
#include <functional>
#include <ostream>
#include "interval.h"

using namespace std;

struct Minimizer
{
    string kmer;
    Interval pos;
    bool operator < ( const Minimizer& str) const;
    Minimizer(string, uint32_t, uint32_t);
    ~Minimizer();
    friend ostream& operator<< (ostream& out, const Minimizer& m); 
};

struct pMiniComp
{
  bool operator()(Minimizer* lhs, Minimizer* rhs);
};


#endif
