#include <iostream>
#include <cstring>
#include <functional>
#include "minimizer.h"
#include "path.h"
#include "interval.h"

using namespace std;

Minimizer::Minimizer(string s, deque<interval> l)
{
    miniKmer = s;
    path = Path();
    path.initialize(l);
    startPosOnString = path.start;
    endPosOnString = path.end;
}

Minimizer::~Minimizer()
{
    //if (path != NULL) {delete path;}
}

bool Minimizer::operator < ( const Minimizer& str) const
{
    if (miniKmer < str.miniKmer) { return true; }
    if ( str.miniKmer < miniKmer ) { return false; }

    if (startPosOnString < str.startPosOnString) { return true; }
    if ( str.startPosOnString < startPosOnString ) { return false; }

    if (endPosOnString < str.endPosOnString) { return true; }
    if ( str.endPosOnString < endPosOnString ) { return false; }

    // if both are completely equal (based on strict weak ordering)
    // then just return false since equality doesn't yield less than
    return false;
}

bool pMiniComp::operator()(Minimizer* lhs, Minimizer* rhs) {
        return (*lhs)<(*rhs);
}

void Minimizer::print() const
{
    cout << "(" << miniKmer << ", ";
    path.print();
    cout << ") ";
}
