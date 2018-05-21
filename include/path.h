#ifndef __PATH_H_INCLUDED__   // if path.h hasn't been included yet...
#define __PATH_H_INCLUDED__

#include <deque>
#include <vector>
#include <iostream>
#include <cstdint>
#include "interval.h"


class Path {
public:
    std::vector<Interval> path;

    void initialize(const std::deque<Interval> &);
    void initialize(const std::vector<Interval> &);
    void initialize(const Interval &);

    uint32_t get_start() const;
    uint32_t get_end() const;
    uint32_t length() const;

    void add_end_interval(const Interval &);

    Path subpath(const uint32_t, const uint32_t) const;

    bool is_branching(const Path &) const;

    bool is_subpath(const Path &) const;

    bool operator<(const Path &y) const;

    bool operator==(const Path &y) const;

    bool operator!=(const Path &y) const;

    friend std::ostream &operator<<(std::ostream &out, const Path &p);

    friend std::istream &operator>>(std::istream &in, Path &p);

    friend Path get_union(const Path&, const Path&);

    friend bool equal_except_null_nodes(const Path &, const Path &);
};

bool equal_except_null_nodes(const Path &, const Path &);

#endif
