#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "interval.h"
#include "path.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class PathTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    //     // (right before the destructor).
  }
};

TEST_F(PathTest,addStartInterval)
{
    deque<Interval> d = {Interval(4,5)};
    Path p;
    p.initialize(d);
    p.add_start_interval(Interval(0,1));
    d.push_front(Interval(0,1));
    EXPECT_ITERABLE_EQ(deque<Interval>, d, p.path);
    EXPECT_DEATH(p.add_start_interval(Interval(3,4)), "");
}

TEST_F(PathTest,addEndInterval)
{
    deque<Interval> d = {Interval(4,5)};
    Path p;
    p.initialize(d);
    p.add_end_interval(Interval(6,9));
    d.push_back(Interval(6,9));
    EXPECT_ITERABLE_EQ(deque<Interval>, d, p.path);
    EXPECT_DEATH(p.add_end_interval(Interval(0,1)), "");
}

TEST_F(PathTest, subpath)
{
    deque<Interval> d, d1;
    d = {Interval(1,3), Interval(4,5), Interval(6,6), Interval(9,40)};
    Path p, p1;
    p.initialize(d);

    // regular
    p1 = p.subpath(1,3);
    d1 = {Interval(1,3), Interval(4,5)};
    EXPECT_ITERABLE_EQ(deque<Interval>, d1, p1.path);

    // handle zero-length interval
    p1 = p.subpath(2,3);
    d1 = {Interval(2,3), Interval(4,5), Interval(6,6), Interval(9,10)};
    EXPECT_ITERABLE_EQ(deque<Interval>, d1, p1.path);

    // start in another interval
    p1 = p.subpath(4,3);
    d1 = {Interval(4,5), Interval(6,6), Interval(9,11)};
    EXPECT_ITERABLE_EQ(deque<Interval>, d1, p1.path);

    // all in one interval
    p1 = p.subpath(10,3);
    d1 = {Interval(10,13)};
    EXPECT_ITERABLE_EQ(deque<Interval>, d1, p1.path);

    // can't get subpath from a coordinate not in path
    //EXPECT_DEATH(p.subpath(0,3), "");
    // can't get subpath of right length if not enough length left in path from start
    //EXPECT_DEATH(p.subpath(39,3), "");
}

TEST_F(PathTest, lessthan)
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

TEST_F(PathTest, equals)
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


