#include <iostream>
#include <functional>
#include "path.h"
#include "interval.h"
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
    path.push_front(i);
    length += i.length;
    start = i.start;
}

void Path::add_end_interval(Interval i)
{
    path.push_back(i);
    length += i.length;
    end = i.end;
}

bool Path::operator < ( const Path& y) const
{
    if (start < y.start) { return true; }
    if (start > y.start) { return false; }
    if ( end < y.end ) { return true; }
    if ( end > y.end ) { return false; }
    return false;
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
