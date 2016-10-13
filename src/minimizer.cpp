#include <iostream>
#include <cstring>
#include <functional>
#include <cassert>
#include "minimizer.h"
#include "path.h"
#include "interval.h"

using namespace std;

Minimizer::Minimizer(string s, uint32_t a, uint32_t b): kmer(s)
{
    pos = Interval(a,b);
    assert(kmer.length()==pos.length);
}

Minimizer::~Minimizer()
{
    //if (path != NULL) {delete path;}
}

bool Minimizer::operator < ( const Minimizer& str) const
{
    if (kmer < str.kmer) { return true; }
    if ( str.kmer < kmer ) { return false; }

    if (pos.start < str.pos.start) { return true; }
    if ( str.pos.start < pos.start ) { return false; }

    if (pos.end < str.pos.end) { return true; }
    if ( str.pos.end < pos.end ) { return false; }

    // if both are completely equal (based on strict weak ordering)
    // then just return false since equality doesn't yield less than
    return false;
}

bool pMiniComp::operator()(Minimizer* lhs, Minimizer* rhs) {
        return (*lhs)<(*rhs);
}

std::ostream& operator<< (std::ostream & out, Minimizer const& m) {
    out << "(" << m.kmer << ", " << m.pos << ")";
    return out ;
}
