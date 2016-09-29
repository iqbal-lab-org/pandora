#include <iostream>
#include "interval.h"

using namespace std;

interval::interval(uint32_t s, uint32_t e) 
{
    start = s;
    end = e;
    length = end - start; // intervals need to be exclusive of end point so that empty strings can be represented
}

void interval::print() const
{
    cout << "[" << start << ", " << end << ") ";
}
