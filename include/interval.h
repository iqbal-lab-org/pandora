#ifndef __INTERVAL_H_INCLUDED__   // if interval.h hasn't been included yet...
#define __INTERVAL_H_INCLUDED__

#include <cstdint> //or <stdint.h>

using namespace std;

struct interval {
    uint32_t start;
    uint32_t end; //in pilot, longest prg was 208,562 characters long
    uint32_t length;

    interval(uint32_t=0, uint32_t=0);
    void print() const;
};
#endif
