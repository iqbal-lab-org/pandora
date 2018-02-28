#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "interval.h"
#include "path.h"
#include <stdint.h>
#include <iostream>

using namespace std;

TEST(PathTest, initialize)
{
    Path p;
    deque<Interval> d = {Interval(0,1), Interval(3,3), Interval(5,10)};
    p.initialize(d);
    EXPECT_EQ(p.path.size(), d.size());
    for(uint i=0; i!=d.size(); ++i)
    {
	EXPECT_EQ(p.path[i], d[i]);
    }
}

TEST(PathTest, length)
{
    Path p;
    deque<Interval> d = {Interval(0,0)};
    p.initialize(d);
    uint j = 0;
    EXPECT_EQ(j, p.length());

    d = {Interval(0,1), Interval(3,3), Interval(5,10)};
    p.initialize(d);
    j = 6;
    EXPECT_EQ(j, p.length());

    d = {Interval(0,1), Interval(3,3)};
    p.initialize(d);
    j = 1;
    EXPECT_EQ(j, p.length());
}

TEST(PathTest,add_start_interval)
{
    deque<Interval> d = {Interval(4,5)};
    Path p;
    p.initialize(d);
    p.add_start_interval(Interval(0,1));
    d.push_front(Interval(0,1));
    EXPECT_ITERABLE_EQ(deque<Interval>, d, p.path);
    EXPECT_DEATH(p.add_start_interval(Interval(3,4)), "");
}

TEST(PathTest,add_end_interval)
{
    deque<Interval> d = {Interval(4,5)};
    Path p;
    p.initialize(d);
    p.add_end_interval(Interval(6,9));
    d.push_back(Interval(6,9));
    EXPECT_ITERABLE_EQ(deque<Interval>, d, p.path);
    EXPECT_DEATH(p.add_end_interval(Interval(0,1)), "");
}

TEST(PathTest, subpath)
{
    deque<Interval> d, d1;
    d = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)};
    Path p, p1;
    p.initialize(d);

    // regular
    p1 = p.subpath(0,3);
    d1 = {Interval(1,3), Interval(4,5)};
    EXPECT_ITERABLE_EQ(deque<Interval>, d1, p1.path);

    // handle zero-length interval
    p1 = p.subpath(1,3);
    d1 = {Interval(2,3), Interval(4,5), Interval(6,6), Interval(9,10)};
    EXPECT_ITERABLE_EQ(deque<Interval>, d1, p1.path);

    // start in another interval
    p1 = p.subpath(2,3);
    d1 = {Interval(4,5), Interval(6,6), Interval(9,11)};
    EXPECT_ITERABLE_EQ(deque<Interval>, d1, p1.path);

    // all in one interval
    p1 = p.subpath(3,3);
    d1 = {Interval(6,6), Interval(9,12)};
    EXPECT_ITERABLE_EQ(deque<Interval>, d1, p1.path);

    // all in one interval
    p1 = p.subpath(4,3);
    d1 = {Interval(10,13)};
    EXPECT_ITERABLE_EQ(deque<Interval>, d1, p1.path);

    // and several null nodes at start of path
    d = {Interval(0,0), Interval(1,1), Interval(3,3), Interval(4,5), Interval(6,6), Interval(9,40)};
    p.initialize(d);
    d1 = {Interval(0,0), Interval(1,1), Interval(3,3), Interval(4,5), Interval(6,6), Interval(9,10)};
    p1 = p.subpath(0,2);
    EXPECT_ITERABLE_EQ(deque<Interval>, d1, p1.path);

    // can't get subpath from a coordinate not in path
    //EXPECT_DEATH(p.subpath(0,3), "");
    // can't get subpath of right length if not enough length left in path from start
    //EXPECT_DEATH(p.subpath(39,3), "");
}

TEST(PathTest, is_branching)
{
    deque<Interval> d, d1;
    d = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)};
    d1 = {Interval(1,3), Interval(4,5), Interval(8,9), Interval(9,40)}; // same number intervals, different intervals
    Path p, p1;
    p.initialize(d);
    p1.initialize(d1);
    EXPECT_EQ(p.is_branching(p), false);
    EXPECT_EQ(p1.is_branching(p1), false);
    EXPECT_EQ(p.is_branching(p1), true);
    EXPECT_EQ(p1.is_branching(p), true);

    d1 = {Interval(4,5), Interval(6,6), Interval(9,47)};
    p1.initialize(d1);
    EXPECT_EQ(p1.is_branching(p1), false);
    EXPECT_EQ(p1.is_branching(p), false);
    EXPECT_EQ(p.is_branching(p1), false);

    d1 = {Interval(0,0), Interval(4,5), Interval(6,6), Interval(9,40)};
    p1.initialize(d1);
    EXPECT_EQ(p.is_branching(p1), true);
    EXPECT_EQ(p1.is_branching(p), true);

    d1 = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(41,50)};
    p1.initialize(d1);
    EXPECT_EQ(p.is_branching(p1), true);
    EXPECT_EQ(p1.is_branching(p), true);
    
    d = {Interval(16810, 16812), Interval(16897, 16904), Interval(16909, 16909), Interval(16914, 16920)};
    d1 = {Interval(16819, 16822), Interval(16897, 16904), Interval(16909, 16909), Interval(16914, 16920)};
    p.initialize(d);
    p1.initialize(d1);
    EXPECT_EQ(p.is_branching(p), false);
    EXPECT_EQ(p1.is_branching(p1), false);
    EXPECT_EQ(p.is_branching(p1), true);
    EXPECT_EQ(p1.is_branching(p), true);

    d = {Interval(37, 52)};
    d1 = {Interval(41, 54), Interval(61, 63)};
    p.initialize(d);
    p1.initialize(d1);
    EXPECT_EQ(p.is_branching(p), false);
    EXPECT_EQ(p1.is_branching(p1), false);
    EXPECT_EQ(p.is_branching(p1), false);
    EXPECT_EQ(p1.is_branching(p), false);
}

TEST(PathTest, less_than)
{
    deque<Interval> d, d1;
    d = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)};
    d1 = {Interval(1,3), Interval(4,5), Interval(8,9), Interval(9,40)}; // same number intervals, different intervals
    Path p, p1;
    p.initialize(d);
    p1.initialize(d1);

    EXPECT_EQ((p<p1), true);
    EXPECT_EQ((p1<p), false);

    d1 = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)}; // identical
    p1.initialize(d1);
    EXPECT_EQ((p<p1), false);
    EXPECT_EQ((p1<p), false);

    d1 = {Interval(1,3), Interval(4,5), Interval(9,40)}; // different number of intervals missing middle one
    p1.initialize(d1);
    EXPECT_EQ((p<p1), true);
    EXPECT_EQ((p1<p), false);

    d1 = {Interval(4,5), Interval(6,6), Interval(9,40)}; // different number of intervals missing first
    p1.initialize(d1);
    EXPECT_EQ((p<p1), true);
    EXPECT_EQ((p1<p), false);

    d1 = {Interval(1,3), Interval(4,6), Interval(6,6), Interval(9,40)}; // different end to one interval
    p1.initialize(d1);
    EXPECT_EQ((p<p1), true);
    EXPECT_EQ((p1<p), false);

    d1 = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(10,40)}; // different start to one interval
    p1.initialize(d1);
    EXPECT_EQ((p<p1), true);
    EXPECT_EQ((p1<p), false);
}

TEST(PathTest, equals)
{
    deque<Interval> d, d1;
    d = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)};
    d1 = {Interval(1,3), Interval(4,5), Interval(8,9), Interval(9,40)}; // same number intervals, different intervals
    Path p, p1;
    p.initialize(d);
    p1.initialize(d1);

    EXPECT_EQ(p,p);
    EXPECT_EQ(p1,p1);
    EXPECT_EQ((p==p1),false);
    EXPECT_EQ((p1==p),false);

    d1 = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)}; // identical
    p1.initialize(d1);
    EXPECT_EQ(p,p1);
    EXPECT_EQ(p1,p);

    d1 = {Interval(1,3), Interval(4,5), Interval(9,40)}; // different number of intervals missing middle one
    p1.initialize(d1);
    EXPECT_EQ((p==p1),false);
    EXPECT_EQ((p1==p),false);

    d1 = {Interval(4,5), Interval(6,6), Interval(9,40)}; // different number of intervals missing first
    p1.initialize(d1);
    EXPECT_EQ((p==p1),false);
    EXPECT_EQ((p1==p),false);

    d1 = {Interval(1,3), Interval(4,6), Interval(6,6), Interval(9,40)}; // different end to one interval
    p1.initialize(d1);
    EXPECT_EQ((p==p1),false);
    EXPECT_EQ((p1==p),false);

    d1 = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(10,40)}; // different start to one interval
    p1.initialize(d1);
    EXPECT_EQ((p==p1),false);
    EXPECT_EQ((p1==p),false);
}

TEST(PathTest, equal_except_null_nodes)
{   
    deque<Interval> d, d1, d2;
    d = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)};
    d1 = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40), Interval(40,40), Interval(59,59)};
    d2 = {Interval(0,0), Interval(1,1), Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)};

    Path p, p1, p2;
    p.initialize(d);
    p1.initialize(d1);
    p2.initialize(d2);

    EXPECT_EQ(equal_except_null_nodes(p,p), true);
    EXPECT_EQ(equal_except_null_nodes(p,p1), true);
    EXPECT_EQ(equal_except_null_nodes(p,p2), true);
    EXPECT_EQ(equal_except_null_nodes(p1,p1), true);
    EXPECT_EQ(equal_except_null_nodes(p1,p2), true);
    EXPECT_EQ(equal_except_null_nodes(p2,p2), true);
}

TEST(PathTest, write)
{
    deque<Interval> d;
    d = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)};
    Path p;
    p.initialize(d);

    stringstream out;
    out << p;
    EXPECT_EQ(out.str(), "4{[1, 3)[4, 5)[6, 6)[9, 40)}"); 
}

TEST(PathTest, read)
{
    deque<Interval> d;
    d = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)};
    Path p, q;
    p.initialize(d);

    stringstream out;
    out << p;

    out >> q;
    EXPECT_EQ(p,q); 
}
