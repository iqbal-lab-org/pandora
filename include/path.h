#ifndef __PATH_H_INCLUDED__   // if path.h hasn't been included yet...
#define __PATH_H_INCLUDED__

#include <deque>
#include <cstdint> //or <stdint.h>
#include <iostream>
#include <functional>
#include "interval.h"

class Path {
  public:
    std::deque<Interval> path;
    uint32_t start;
    uint32_t end;

    Path();
    ~Path();
    Path(const Path&);
    Path& operator=(const Path&);

    void initialize(const std::deque<Interval>&);
    uint32_t length() const;
    void add_start_interval(const Interval&);
    void add_end_interval(const Interval&);
    Path subpath(const uint32_t, const uint32_t) const;
    bool is_branching(const Path& y) const;
    bool operator < (const Path& y) const;
    bool operator == (const Path& y) const;
    bool operator != (const Path& y) const;
  friend std::ostream& operator<< (std::ostream& out, const Path& p); 
  friend std::istream& operator>> (std::istream& in, Path& p);
  friend bool equal_except_null_nodes (const Path&, const Path&);
};

bool equal_except_null_nodes (const Path&, const Path&);
#endif
