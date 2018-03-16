#include <iostream>
#include <cassert>
#include "interval.h"

using namespace std;

Interval::Interval(uint32_t s, uint32_t e) : start(s) {
    assert(e >= start); // not a real interval ;
    assert (e - start < std::numeric_limits<uint16_t>::max());
    length = (uint16_t) (e - start); // intervals need to be exclusive of end point so that empty strings can be represented
}

uint32_t Interval::get_end() const {
    return start + (uint32_t) length;
}

ostream &operator<<(ostream &out, Interval const &i) {
    out << "[" << i.start << ", " << i.get_end() << ")";
    return out;
}

istream &operator>>(istream &in, Interval &i) {
    uint32_t end;
    in.ignore(1, '[');
    in >> i.start;
    in.ignore(2, ' ');
    in >> end;
    in.ignore(1, ')');
    i.length = (uint16_t) (end - i.start);
    return in;
}

bool Interval::operator==(const Interval &y) const {
    return (start == y.start and length == y.length);
}

bool Interval::operator!=(const Interval &y) const {
    return !(start == y.start and length == y.length);
}

bool Interval::operator<(const Interval &y) const {
    if (start < y.start) { return true; }
    if (start > y.start) { return false; }
    if (length < y.length) { return true; }
    if (length > y.length) { return false; }
    return false;
}
