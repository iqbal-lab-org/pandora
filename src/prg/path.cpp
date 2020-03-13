#include <cassert>
#include <localPRG.h>
#include "prg/path.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

uint32_t prg::Path::get_start() const
{
    if (path.size() < 1)
        return 0;
    return path.front().start;
}

uint32_t prg::Path::get_end() const
{
    if (path.size() < 1)
        return 0;
    return path.back().start + (uint32_t)path.back().length;
}

uint32_t prg::Path::length() const
{
    uint32_t length = 0;
    for (const auto& interval : path)
        length += interval.length;
    return length;
}

void prg::Path::add_end_interval(const Interval& i)
{
    memoizedDirty = true;
    assert(i.start >= get_end()
        || assert_msg("tried to add interval starting at "
            << i.start << " to end of path finishing at " << get_end()));
    path.push_back(i);
}

std::vector<LocalNodePtr> prg::Path::nodes_along_path(const LocalPRG& localPrg)
{
    // sanity check
    assert(
        (isMemoized == false
            || (isMemoized == true && localPRGIdOfMemoizedLocalNodePath == localPrg.id))
        || assert_msg("Memoization bug: memoized a local node path for PRG with id"
            << localPRGIdOfMemoizedLocalNodePath << " but PRG id " << localPrg.id
            << " is also trying to use this memoized path"));

    if (isMemoized == false
        || memoizedDirty == true) { // checks if we must do memoization
        // yes
        // not memoized, memoize
        memoizedLocalNodePath = localPrg.nodes_along_path_core(*this);

        // update the memoization control variables
        localPRGIdOfMemoizedLocalNodePath = localPrg.id;
        isMemoized = true;
        memoizedDirty = false;

        // return the local node path
        return memoizedLocalNodePath;
    } else if (memoizedDirty == false) {
        // redudant call, return the memoized local node path
        return memoizedLocalNodePath;
    } else {
        fatal_error("Bug on prg::Path::nodes_along_path()");
    }
}

prg::Path prg::Path::subpath(const uint32_t start, const uint32_t len) const
{
    // function now returns the path starting at position start along the path, rather
    // than at position start on linear PRG, and for length len cout << "find subpath of
    // " << *this << " from " << start << " for length " << len << endl;
    assert(start + len <= length());
    prg::Path p;
    std::deque<Interval> d;
    uint32_t covered_length = 0;
    uint32_t added_len = 0;
    for (const auto& interval : path) {
        if ((covered_length <= start and covered_length + interval.length > start
                and p.path.empty())
            or (covered_length == start and interval.length == 0 and p.path.empty())) {
            assert(added_len == 0);
            d = { Interval(interval.start + start - covered_length,
                std::min(interval.get_end(),
                    interval.start + start - covered_length + len - added_len)) };
            p.initialize(d);
            added_len
                += std::min(len - added_len, interval.length - start + covered_length);
            /// cout << "added first interval " << p.path.back() << " and added length
            /// is now " << added_len << endl;
        } else if (covered_length >= start and added_len <= len) {
            p.add_end_interval(Interval(interval.start,
                std::min(interval.get_end(), interval.start + len - added_len)));
            added_len += std::min(len - added_len, interval.length);
            // cout << "added interval " << p.path.back() << " and added length is now "
            // << added_len << endl;
        }
        covered_length += interval.length;
        // cout << "covered length is now " << covered_length << endl;
        if (added_len >= len) {
            break;
        }
    }
    assert(added_len == len);
    return p;
}

bool prg::Path::is_branching(const prg::Path& y)
    const // returns true if the two paths branch together or coalesce apart
{
    // simple case, one ends before the other starts -> return false
    if (get_end() < y.get_start() or y.get_end() < get_start()) {
        return false;
    }

    // otherwise, check it out
    bool overlap = false;
    std::vector<Interval>::const_iterator it, it2;
    for (it = path.begin(); it != path.end(); ++it) {
        if (overlap) {
            if (it->start != it2->start) {
                // had paths which overlapped and now don't
                // cout << "branch" << endl;
                return true; // represent branching paths at this point
            }
            ++it2;
            if (it2 == y.path.end()) {
                return false;
                break;
            }
        } else {
            for (it2 = y.path.begin(); it2 != y.path.end(); ++it2) {
                // cout << *it << " " << *it2 << endl;
                if ((it->get_end() > it2->start and it->start < it2->get_end())
                    or (*it == *it2)) {
                    // then the paths overlap
                    overlap = true;
                    // cout << "overlap" << endl;
                    // first query the previous intervals if they exist, to see if they
                    // overlap
                    if (it != path.begin() && it2 != y.path.begin()
                        && (--it)->get_end() != (--it2)->get_end()) {
                        // cout << "coal" << endl;
                        return true; // represent coalescent paths at this point
                    } else {
                        ++it2;
                        if (it2 == y.path.end()) {
                            return false;
                        }
                        break; // we will step through intervals from here comparing to
                               // path
                    }
                }
            }
        }
    }

    return false;
}

bool prg::Path::is_subpath(const prg::Path& big_path) const
{
    if (big_path.length() < length() or big_path.get_start() > get_start()
        or big_path.get_end() < get_end() or is_branching(big_path)) {
        // cout << "fell at first hurdle " << (big_path.length() < length());
        // cout << (big_path.start > start) << (big_path.end < end);
        // cout << (is_branching(big_path)) << endl;
        return false;
    }

    uint32_t offset = 0;
    for (const auto& big_i : big_path.path) {
        if (big_i.get_end() >= get_start()) {
            offset += get_start() - big_i.start;
            if (offset + length() > big_path.length()) {
                return false;
            } else if (big_path.subpath(offset, length()) == *this) {
                return true;
            }
            break;
        }
        offset += big_i.length;
    }

    return false;
}

bool prg::Path::operator<(const prg::Path& y) const
{
    auto it2 = y.path.begin();
    auto it = path.begin();
    while (it != path.end() and it2 != y.path.end()) {
        if (!(*it
                == *it2)) // for the first interval which is not the same in both paths
        {
            return (*it < *it2); // return interval comparison
        }
        it++;
        it2++;
    }
    if (it == path.end() and it2 != y.path.end()) {
        // if path is shorter than comparison path, but equal otherwise, return that it
        // is smaller
        return true;
    }

    return false; // shouldn't ever call this
}

bool prg::Path::operator==(const prg::Path& y) const
{
    if (path.size() != y.path.size()) {
        return false;
    }
    auto it2 = y.path.begin();
    for (auto it = path.begin(); it != path.end();) {
        if (!(*it == *it2)) {
            return false;
        }
        it++;
        it2++;
    }
    return true;
}

// tests if the paths are equal except for null nodes at the start or
// ends of the paths
bool equal_except_null_nodes(const prg::Path& x, const prg::Path& y)
{
    auto it2 = y.begin();
    for (auto it = x.begin(); it != x.end();) {
        while (it != x.end() and it->length == 0) {
            it++;
        }
        while (it2 != y.end() and it2->length == 0) {
            it2++;
        }

        if (it == x.end() and it2 == y.end()) {
            break;
        } else if (it == x.end() or it2 == y.end()) {
            return false;
        }

        if (it->length > 0 and it2->length > 0 and !(*it == *it2)) {
            return false;
        } else {
            it++;
            it2++;
        }
    }
    return true;
}

bool prg::Path::operator!=(const prg::Path& y) const { return (!(path == y.path)); }

std::ostream& prg::operator<<(std::ostream& out, const prg::Path& p)
{
    uint32_t num_intervals = p.path.size();
    out << num_intervals << "{";
    for (auto it = p.path.begin(); it != p.path.end(); ++it) {
        out << *it;
    }
    out << "}";
    return out;
}

std::istream& prg::operator>>(std::istream& in, prg::Path& p)
{
    uint32_t num_intervals;
    in >> num_intervals;
    std::deque<Interval> d(num_intervals, Interval());
    in.ignore(1, '{');
    for (uint32_t i = 0; i != num_intervals; ++i) {
        in >> d[i];
    }
    in.ignore(1, '{');
    p.initialize(d);
    return in;
}

prg::Path prg::get_union(const prg::Path& x, const prg::Path& y)
{
    auto xit = x.path.begin();
    auto yit = y.path.begin();

    prg::Path p;
    assert(x < y);

    if (x.get_end() < y.get_start() or x.is_branching(y)) {
        return p;
    } else if (x.path.empty()) {
        return y;
    }

    while (
        xit != x.path.end() and yit != y.path.end() and xit->get_end() < yit->start) {
        if (p.path.empty()) {
            p.initialize(*xit);
        } else {
            p.add_end_interval(*xit);
        }
        xit++;
    }
    if (xit != x.path.end() and yit != y.path.end() and xit->start <= yit->get_end()) {
        // then we have overlap
        if (p.path.empty()) {
            p.initialize(
                Interval(xit->start, std::max(yit->get_end(), xit->get_end())));
        } else {
            p.add_end_interval(
                Interval(xit->start, std::max(yit->get_end(), xit->get_end())));
        }
        while (yit != --y.path.end()) {
            yit++;
            p.add_end_interval(*yit);
        }
    }
    return p;
}