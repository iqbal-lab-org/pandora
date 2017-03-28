#include <iostream>
#include <functional>
#include "path.h"
#include "interval.h"
#include "errormessages.h"
#include <cassert>
#include <deque>

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

Path::Path(){}

uint32_t Path::length() const
{
    uint32_t i = 0;
    for (deque<Interval>::const_iterator it = path.begin(); it!=path.end(); ++it)
    {
	i+=it->length;
    }
    return i;
}

void Path::initialize(const deque<Interval>& q)
{
    path = q;
    start = q.begin()->start;
    //cout << "Paths starts at " << start;
    end = (*--q.end()).end;
    //cout << "Paths ends at " << end;
    end = path.back().end;
}

/*Path::~Path()
{
    path.clear();
}*/

void Path::add_start_interval(const Interval& i)
{
    assert (i.end <= path.begin()->start || assert_msg(i.end << ">" << path.begin()->start));
    path.push_front(i);
    start = i.start;
}

void Path::add_end_interval(const Interval& i)
{
    assert (i.start >= path.back().end || assert_msg("tried to add interval starting at " << i.start << " to end of path finishing at " << path.back().end));
    path.push_back(i);
    end = i.end;
}

Path Path::subpath(const uint32_t start, const uint32_t len) const
{
    //function now returns the path starting at position start along the path, rather than at position start on 
    //linear PRG, and for length len
    //cout << "find subpath of " << *this << " from " << start << " for length " << len << endl;
    assert(start+len <= length());
    Path p;
    deque<Interval> d;
    uint32_t covered_length = 0;
    uint32_t added_len = 0;
    for (deque<Interval>::const_iterator it=path.begin(); it!=path.end(); ++it)
    {
	if ((covered_length <= start and covered_length + it->length > start and p.path.size() == 0) or (covered_length == start and it->length == 0))
	{
            assert(added_len == 0);
	    d = {Interval(it->start + start - covered_length, min(it->end, it->start + start - covered_length + len - added_len))};
	    p.initialize(d);
	    added_len += min(len - added_len, it->length - start + covered_length);
            ///cout << "added first interval " << p.path.back() << " and added length is now " << added_len << endl;
	} else if (covered_length >= start and added_len <= len)
	{
	    p.add_end_interval(Interval(it->start, min(it->end, it->start + len - added_len)));
	    added_len += min(len - added_len, it->length);
	    //cout << "added interval " << p.path.back() << " and added length is now " << added_len << endl;
	}
	covered_length += it->length;
	//cout << "covered length is now " << covered_length << endl;
	if (added_len >= len)
	{
	    break;
	}
    }
    assert(added_len == len);
    return p;
}

/*bool Path::is_disjoint(const Path& y) const // does the path overlap this path
{
}*/

bool Path::is_branching(const Path& y) const // returns true if the two paths branch together or coalesce apart
{

    // simple case, one ends before the other starts -> return false
    if (end < y.start or y.end < start)
    {
	return false;
    }

    // otherwise, check it out
    bool overlap = false;
    deque<Interval>::const_iterator it, it2;
    for (it=path.begin(); it!=path.end(); ++it)
    {
        if (overlap == true)
        {
            if (it->start != it2->start)
            {
                // had paths which overlapped and now don't
                //cout << "branch" << endl;
                return true; // represent branching paths at this point
            }
	    ++it2;
            if (it2==y.path.end())
            {
		return false;
		break;
	    }
        } else {
            for (it2=y.path.begin(); it2!=y.path.end(); ++it2)
	    {
	        //cout << *it << " " << *it2 << endl;
	        if ((it->end > it2->start and it->start < it2->end) or (*it==*it2))
	        {
	            // then the paths overlap
	            overlap = true;
		    //cout << "overlap" << endl;
	            // first query the previous intervals if they exist, to see if they overlap
	            if (it!=path.begin() && it2!=y.path.begin() && (--it)->end!=(--it2)->end)
		    {
	  	        //cout << "coal" << endl;
		        return true; // represent coalescent paths at this point
		    } else {
			++it2;
		        if (it2==y.path.end())
            		{
			    return false;
			}
		        break; // we will step through intervals from here comparing to path
		    }
	        }
	    }
	}
    }

    return false;
}

bool Path::operator < ( const Path& y) const
{
    std::deque<Interval>::const_iterator it2=y.path.begin();
    std::deque<Interval>::const_iterator it=path.begin();
    while (it!=path.end() and it2!=y.path.end())
    {
        if (!(*it==*it2)) //for the first interval which is not the same in both paths
	{
	    return (*it<*it2); // return interval comparison
	}
	it++;
        it2++;
    }
    if (it==path.end() and it2!=y.path.end())
    {
	// if path is shorter than comparison path, but equal otherwise, return that it is smaller
	return true;
    }

    return false; // shouldn't ever call this
}

bool Path::operator == ( const Path& y) const
{
    if (path.size() != y.path.size()) { return false;}
    std::deque<Interval>::const_iterator it2=y.path.begin();
    for (std::deque<Interval>::const_iterator it=path.begin(); it!=path.end();)
    {
        if (!(*it==*it2)) {return false;}
        it++;
        it2++;
    }
    return true;
}

bool Path::operator != ( const Path& y) const
{
    return (! (path == y.path));
}

std::ostream& operator<< (std::ostream & out, Path const& p) {
    uint32_t num_intervals = p.path.size();
    out << num_intervals << "{";
    for (std::deque<Interval>::const_iterator it=p.path.begin(); it!=p.path.end(); ++it)
    {
        out << *it;
    }
    out << "}";
    return out ;
}

std::istream& operator>> (std::istream & in, Path& p) {
    uint32_t num_intervals;
    in >> num_intervals;
    deque<Interval> d(num_intervals, Interval());
    in.ignore(1,'{');
    for (uint32_t i = 0; i != num_intervals; ++i)
    {
        in >> d[i];
    }
    in.ignore(1,'{');
    p.initialize(d); 
    return in;
}
