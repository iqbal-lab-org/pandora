#ifndef __PATH_H_INCLUDED__   // if path.h hasn't been included yet...
#define __PATH_H_INCLUDED__

#include <deque>
#include <cstdint> //or <stdint.h>
#include "interval.h"

using namespace std;

struct interval;

class Path {
  public:
    deque<interval> path;
    uint32_t length;
    uint32_t start;
    uint32_t end;

    Path();
    ~Path();
    void add_start_interval(interval);
    void add_end_interval(interval);
    void initialize(deque<interval>);
    void print() const;
};

#endif
