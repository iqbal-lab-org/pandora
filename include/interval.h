#ifndef __INTERVAL_H_INCLUDED__   // if interval.h hasn't been included yet...
#define __INTERVAL_H_INCLUDED__

#include <cstdint>
#include <iostream>


struct Interval {
    uint32_t start;
    uint32_t length;
    uint32_t end;
    // in pilot, longest prg was 208,562 characters long

    Interval(uint32_t= 0, uint32_t= 0);

    uint32_t get_end() const;

    friend std::ostream &operator<<(std::ostream &out, const Interval &i);

    friend std::istream &operator>>(std::istream &in, Interval &i);

    bool operator==(const Interval &y) const;

    bool operator!=(const Interval &y) const;

    bool operator<(const Interval &y) const;
};

#endif
