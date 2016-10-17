#include <iostream>
#include <functional>
#include "path.h"
#include "interval.h"
#include "errormessages.h"
#include <cassert>
#include <deque>

using namespace std;

Path::Path(){ }

void Path::initialize(deque<Interval> q)
{
    path = q;
    length = 0;
    start = q.begin()->start;
    //cout << "Paths starts at " << start;
    end = (*--q.end()).end;
    //cout << "Paths ends at " << end;
    for (std::deque<Interval>::iterator it=path.begin(); it!=path.end(); ++it)
    {
	length += (*it).length;
    }
    end = path.back().end;
}

Path::~Path()
{
    path.clear();
}

void Path::add_start_interval(Interval i)
{
    assert (i.end <= path.begin()->start);// && PATH_ERROR);
    path.push_front(i);
    length += i.length;
    start = i.start;
}

void Path::add_end_interval(Interval i)
{
    assert (i.start >= path.back().end);// && PATH_ERROR);
    path.push_back(i);
    length += i.length;
    end = i.end;
}

Path Path::subpath(uint32_t start, uint32_t len) const
{
   Path p;
   uint32_t added_len = 0;
   for (deque<Interval>::const_iterator it=path.begin(); it!=path.end(); ++it)
   {
	if (it->start <= start and it->end > start)
	{
	    // first interval to add
	    deque<Interval> d = {Interval(start, min(start+len, it->end))};
	    p.initialize(d);
	    added_len += min(start+len, it->end) - start;
	} else if (it->length < len - added_len and p.path.size() > 0){ //check p initialized with a path
	    p.add_end_interval(*it);
	    added_len += it->length;
	} else if (it->length >= len - added_len and p.path.size() > 0) { //check p initialised with a path
	    // last interval to add
	    p.add_end_interval(Interval(it->start, it->start + len - added_len));
	    added_len = len;
	    return p;
	}
   }
   assert (added_len == len);// && SUBPATH_ERROR);
   return p; // this should never happen
}

bool Path::operator < ( const Path& y) const
{
    if (start < y.start) { return true; }
    if (start > y.start) { return false; }
    if ( end < y.end ) { return true; }
    if ( end > y.end ) { return false; }
    return false;
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

std::ostream& operator<< (std::ostream & out, Path const& p) {
    out << "{";
    for (std::deque<Interval>::const_iterator it=p.path.begin(); it!=p.path.end(); ++it)
    {
        out << *it;
    }
    out << "}";
    return out ;
}
