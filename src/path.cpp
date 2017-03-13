#include <iostream>
#include <functional>
#include "path.h"
#include "interval.h"
#include "errormessages.h"
#include <cassert>
#include <deque>

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

Path::Path(){
    length = 0; 
}

void Path::initialize(const deque<Interval>& q)
{
    path = q;
    start = q.begin()->start;
    length = 0;
    //cout << "Paths starts at " << start;
    end = (*--q.end()).end;
    //cout << "Paths ends at " << end;
    for (std::deque<Interval>::const_iterator it=path.begin(); it!=path.end(); ++it)
    {
	length += (*it).length;
    }
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
    length += i.length;
    start = i.start;
}

void Path::add_end_interval(const Interval& i)
{
    assert (i.start >= path.back().end || assert_msg("tried to add interval starting at " << i.start << " to end of path finishing at " << path.back().end));
    path.push_back(i);
    length += i.length;
    end = i.end;
}

Path Path::subpath(const uint32_t start, const uint32_t len) const
{
    //cout << "start subpath function with : " << start << " " << len << " for path " << *this << endl;
    Path p;
    deque<Interval> d;
    uint32_t added_len = 0;
    for (deque<Interval>::const_iterator it=path.begin(); it!=path.end(); ++it)
    {
	if ((it->start <= start and it->end > start) or (it->start == start and it->length == 0))
	{
	    //cout << "start of subpath" << endl;
	    // first interval to add
	    d = {Interval(start, min(start+len, it->end))};
	    p.initialize(d);
	    added_len += min(start+len, it->end) - start;
            if (added_len == len)
            { return p;}
	} else if (it->length < len - added_len and p.path.size() > 0){ //check p initialized with a path
	    //cout << "mid subpath interval" << endl;
	    p.add_end_interval(*it);
	    added_len += it->length;
	} else if (it->length >= len - added_len and p.path.size() > 0) { //check p initialised with a path
	    //cout << "end of subpath" << endl;
	    // last interval to add
	    p.add_end_interval(Interval(it->start, it->start + len - added_len));
            //cout << "p is now: " << p << endl;
	    added_len = len;
	    return p;
	} else if (it->start > start and p.path.size() == 0) {
	    // start doesn't lie in an interval of the path
	    return p;
	}
	
    }
    //cout << "after subpath processing" << endl;
    // should only be here if added_len <len
    // return empty path, not half a path
    assert(added_len < len);
    d.clear();
    p.initialize(d);
    return p; // this should never happen
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
