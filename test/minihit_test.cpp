#include "gtest/gtest.h"
#include "minihit.h"
#include "minimizer.h"
#include "minirecord.h"
#include "interval.h"
#include "path.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class MinimizerHitTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(MinimizerHitTest,create){
    Minimizer* m;
    m = new Minimizer("hello", 0,5);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p = Path();
    p.initialize(d);
    MiniRecord mr = MiniRecord(0,p);
    MinimizerHit mh = MinimizerHit(1, m, mr, 0);
    uint32_t j = 1;
    EXPECT_EQ(j, mh.read_id);
    EXPECT_EQ(Interval(0,5), mh.read_interval);
    j=0;
    EXPECT_EQ(j, mh.prg_id);
    EXPECT_EQ(p, mh.prg_path);
    EXPECT_EQ(j, mh.strand);

    m = new Minimizer("hello",1,6);
    EXPECT_DEATH(MinimizerHit(1, m, mr, 0), "");
}

TEST_F(MinimizerHitTest,equals){
    Minimizer* m;
    m = new Minimizer("hello", 0,5);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p = Path();
    p.initialize(d);
    MiniRecord mr = MiniRecord(0,p);
    MinimizerHit mh1 = MinimizerHit(1, m, mr, 0);

    m = new Minimizer("hello", 0,5);
    d = {Interval(7,9), Interval(11, 14)};
    p.initialize(d);
    mr = MiniRecord(0,p);
    MinimizerHit mh2 = MinimizerHit(1, m, mr, 0);

    EXPECT_EQ(mh1, mh1);
    EXPECT_EQ(mh2, mh2);
    EXPECT_EQ((mh1==mh2), false);
}

TEST_F(MinimizerHitTest,compare){
    set<MinimizerHit> hits;

    Minimizer* m;
    m = new Minimizer("hello", 1,6);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p = Path();
    p.initialize(d);
    MiniRecord mr = MiniRecord(0,p);
    MinimizerHit mh1 = MinimizerHit(1, m, mr, 0);

    MinimizerHit mh2 = MinimizerHit(2, m, mr, 0);

    m = new Minimizer("hello", 0,5);
    MinimizerHit mh3 = MinimizerHit(1, m, mr, 0);

    d = {Interval(6,10), Interval(11, 14)};
    p.initialize(d);
    mr = MiniRecord(0,p);
    MinimizerHit mh4 = MinimizerHit(1, m, mr, 0);

    d = {Interval(6,10), Interval(12, 15)};
    p.initialize(d);
    mr = MiniRecord(0,p);
    MinimizerHit mh5 = MinimizerHit(1, m, mr, 0);

    hits.insert(mh1);
    hits.insert(mh2);
    hits.insert(mh3);
    hits.insert(mh4);
    hits.insert(mh5);
    
    vector<MinimizerHit> expected;
    expected.push_back(mh1);
    expected.push_back(mh2);
    expected.push_back(mh3);
    expected.push_back(mh4);
    expected.push_back(mh5);

    uint32_t j = 0;
    for (set<MinimizerHit>::iterator it=hits.begin(); it!=hits.end(); ++it)
    {
        EXPECT_EQ(expected[j], *it);
        j++;
    }
}

