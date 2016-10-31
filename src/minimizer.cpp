#include <iostream>
#include <cstring>
#include <functional>
#include <cassert>
#include <cmath>
#include <stdint.h>
#include "minimizer.h"
#include "path.h"
#include "interval.h"

using namespace std;

Minimizer::Minimizer(uint64_t s, uint32_t a, uint32_t b): kmer(s)
{
    pos = Interval(a,b);
    assert(s<=pow(4,pos.length)); // used to check kmer length same as interval length. 
					  // Can't any more but at least know if s is too big for a kmer of interval size to have generated it.
}

Minimizer::~Minimizer()
{
    //if (path != NULL) {delete path;}
}

bool Minimizer::operator < ( const Minimizer& y) const
{
    if (kmer < y.kmer) { return true; }
    if ( y.kmer < kmer ) { return false; }

    if (pos.start < y.pos.start) { return true; }
    if ( y.pos.start < pos.start ) { return false; }

    if (pos.end < y.pos.end) { return true; }
    if ( y.pos.end < pos.end ) { return false; }

    // if both are completely equal (based on strict weak ordering)
    // then just return false since equality doesn't yield less than
    return false;
}

bool Minimizer::operator == (const Minimizer& y) const {
    if (kmer != y.kmer) {return false;}
    if (!(pos == y.pos)) {return false;}
    return true;
}

bool pMiniComp::operator()(Minimizer* lhs, Minimizer* rhs) {
        return (*lhs)<(*rhs);
}

std::ostream& operator<< (std::ostream & out, Minimizer const& m) {
    out << "(" << m.kmer << ", " << m.pos << ")";
    return out ;
}
