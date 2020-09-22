#include "interval.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

Interval::Interval(uint32_t s, uint32_t e)
    : start(s)
{
    assert(e >= start); // not a real interval ;
    assert(e - start < std::numeric_limits<uint32_t>::max()
        or assert_msg("Interval [" << s << "," << e << ") was too long"));
    length = e - start; // intervals need to be exclusive of end point so that empty
                        // strings can be represented
}

uint32_t Interval::get_end() const { return start + length; }

std::ostream& operator<<(std::ostream& out, Interval const& i)
{
    out << "[" << i.start << ", " << i.get_end() << ")";
    return out;
}

std::istream& operator>>(std::istream& in, Interval& i)
{
    uint32_t end;
    in.ignore(1, '[');
    in >> i.start;
    in.ignore(2, ' ');
    in >> end;
    in.ignore(1, ')');
    assert(end >= i.start);
    assert(end - i.start < std::numeric_limits<uint32_t>::max()
        or assert_msg("Interval [" << i.start << "," << end << ") was too long"));
    i.length = end - i.start;
    return in;
}

bool Interval::operator==(const Interval& y) const
{
    return (start == y.start and length == y.length);
}

bool Interval::operator!=(const Interval& y) const
{
    return !(start == y.start and length == y.length);
}

bool Interval::operator<(const Interval& y) const
{
    if (start < y.start) {
        return true;
    }
    if (start > y.start) {
        return false;
    }
    if (length < y.length) {
        return true;
    }
    if (length > y.length) {
        return false;
    }
    return false;
}

bool Interval::empty() const { return length == 0; }

void merge_intervals_within(std::vector<Interval>& intervals, const uint32_t dist)
{
    if (intervals.size() < 2) { // nothing to merge
        return;
    }
    std::sort(intervals.begin(), intervals.end());

    // assumption that iv1 < iv2
    auto close { [dist](const Interval& iv1, const Interval& iv2) -> bool {
        const auto leading_edge { iv1.get_end() + dist };
        return leading_edge > iv2.start;
    } };

    unsigned long prev_idx { 0 };
    for (std::size_t i = 1; i < intervals.size(); ++i) {
        const auto& current_iv { intervals[i] };
        auto& prev_iv { intervals[prev_idx] };

        if (close(prev_iv, current_iv)) {
            const auto new_end { std::max(prev_iv.get_end(), current_iv.get_end()) };
            intervals[prev_idx] = Interval(prev_iv.start, new_end);

            // remove the current interval as it was merged into the previous interval
            intervals.erase(intervals.begin() + i);
            --i; // vector elements get relocated after erase so stay at current idx
        } else {
            ++prev_idx;
        }
    }

    intervals.resize(prev_idx + 1);
}
