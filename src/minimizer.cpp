#include <iostream>
#include <cstring>
#include <functional>
#include "minimizer.h"
#include "path.h"
#include "interval.h"

using namespace std;

Minimizer::Minimizer(string s, list<interval> l)
{
    miniWord = s;
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
    if (miniWord < str.miniWord) { return true; }
    if ( str.miniWord < miniWord ) { return false; }

    if (startPosOnString < str.startPosOnString) { return true; }
    if ( str.startPosOnString < startPosOnString ) { return false; }

    if (endPosOnString < str.endPosOnString) { return true; }
    if ( str.endPosOnString < endPosOnString ) { return false; }

    // if both are completely equal (based on strict weak ordering)
    // then just return false since equality doesn't yield less than
    return false;
}

void Minimizer::print() const
{
    cout << "(" << miniWord << ", ";
    path.print();
    cout << ") ";
}
