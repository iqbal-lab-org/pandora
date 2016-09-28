#ifndef __PATH_H_INCLUDED__   // if path.h hasn't been included yet...
#define __PATH_H_INCLUDED__

#include <list>
#include <cstdint> //or <stdint.h>
#include "interval.h"

using namespace std;

struct interval;

class Path {
  public:
    list<interval> path;
    uint32_t length;
    uint32_t start;
    uint32_t end;

    Path(list<interval>);
    ~Path();
    void add_start_interval(interval);
    void add_end_interval(interval);
    void print() const;
};

#endif
