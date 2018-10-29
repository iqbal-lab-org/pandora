#include <limits>
#include <cassert>
#include "interval.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


Interval::Interval(uint32_t s, uint32_t e) : start(s), end(e) {
    assert(this->end >= this->start);
    assert(this->end - this->start <= std::numeric_limits<uint32_t>::max());

    // intervals must be end exclusive so that empty strings can be represented
    this->length = this->end - this->start;
}

uint32_t Interval::get_end() const {
    return start + length;
}

std::ostream &operator<<(std::ostream &out, Interval const &i) {
    out << "[" << i.start << ", " << i.get_end() << ")";
    return out;
}

std::istream &operator>>(std::istream &in, Interval &i) {
    uint32_t end;
    in.ignore(1, '[');
    in >> i.start;
    in.ignore(2, ' ');
    in >> i.end;
    in.ignore(1, ')');

    i.length = i.end - i.start;
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
