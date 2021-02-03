#ifndef __INTERVAL_H_INCLUDED__ // if interval.h hasn't been included yet...
#define __INTERVAL_H_INCLUDED__

#include <cstdint>
#include <iostream>
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include "fatal_error.h"

struct Interval {
    uint32_t start;
    uint32_t length; // in pilot, longest prg was 208,562 characters long

    Interval(uint32_t s = 0, uint32_t e = 0);

    // non-inclusive
    uint32_t get_end() const;

    friend std::ostream& operator<<(std::ostream& out, const Interval& i);

    friend std::istream& operator>>(std::istream& in, Interval& i);

    bool operator==(const Interval& y) const;

    bool operator!=(const Interval& y) const;

    bool operator<(const Interval& y) const;

    bool empty() const;

    bool is_close(const Interval& other, uint32_t dist = 0) const;

    static bool sorted_interval_vector_has_overlapping_intervals (const std::vector<Interval> &intervals);
};

// Merge intervals within dist of each other. Changes the vector inplace.
// Time complexity: O(N log(N)) Sort the vector, then a single linear pass through
void merge_intervals_within(std::vector<Interval>&, uint32_t);

#endif
