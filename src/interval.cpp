#include "interval.h"

Interval::Interval(uint32_t s, uint32_t e)
    : start(s)
{
    if (e < start) {
        fatal_error("Error when building interval: interval end cannot be less than the interval start");
    }
    // intervals need to be exclusive of end so that empty strings can be represented
    length = e - start;
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

bool Interval::is_close(const Interval& other, const uint32_t dist) const
{
    const auto iv1 = std::min(*this, other);
    const auto iv2 = std::max(other, *this);

    const auto leading_edge { iv1.get_end() + dist };
    return leading_edge > iv2.start;
}

void merge_intervals_within(std::vector<Interval>& intervals, const uint32_t dist)
{
    if (intervals.size() < 2) { // nothing to merge
        return;
    }
    std::sort(intervals.begin(), intervals.end());

    unsigned long prev_idx { 0 };
    for (std::size_t i = 1; i < intervals.size(); ++i) {
        const auto& current_iv { intervals[i] };
        auto& prev_iv { intervals[prev_idx] };

        if (prev_iv.is_close(current_iv, dist)) {
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

bool Interval::sorted_interval_vector_has_overlapping_intervals (const std::vector<Interval> &intervals) {
    for (uint32_t index = 1; index < intervals.size(); ++index) {
        bool there_is_overlap = intervals[index - 1].get_end() > intervals[index].start;
        if (there_is_overlap) {
            return true;
        }
    }
    return false;
}