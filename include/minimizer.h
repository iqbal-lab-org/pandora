#ifndef __MINIMIZER_H_INCLUDED__   // if minimizer.h hasn't been included yet...
#define __MINIMIZER_H_INCLUDED__

#include <string>
#include <functional>
#include <ostream>
#include "path.h"
#include "interval.h"

using namespace std;

struct Minimizer
{
    string miniKmer;
    Path path;
    uint32_t startPosOnString; // actual position along prg, including numbers
    uint32_t endPosOnString;
    bool operator < ( const Minimizer& str) const;
    Minimizer(string, deque<Interval>);
    ~Minimizer();
    friend ostream& operator<< (ostream& out, const Minimizer& m); 
};

struct pMiniComp
{
  bool operator()(Minimizer* lhs, Minimizer* rhs);
};


#endif
