#include "gtest/gtest.h"
#include "minirecord.h"
#include "path.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class MiniRecordTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(MiniRecordTest,create){
    deque<Interval> v1 = {Interval(0,5)};
    deque<Interval> v2 = {Interval(1,4), Interval(15,17)};
    deque<Interval> v3 = {Interval(1,6)};
    deque<Interval> v4 = {Interval(0,3), Interval(16,18)};

    Path p;
    p.initialize(v1);
    MiniRecord m1(1,p,0);
    uint32_t j=1;
    EXPECT_EQ(j, m1.prg_id);
    EXPECT_EQ(p, m1.path);
    p.initialize(v2);
    MiniRecord m2(2,p,0);
    j=2;
    EXPECT_EQ(j, m2.prg_id);
    EXPECT_EQ(p, m2.path);
    p.initialize(v3);
    MiniRecord m3(3,p,0);
    j=3;
    EXPECT_EQ(j, m3.prg_id);
    EXPECT_EQ(p, m3.path);
    p.initialize(v4);
    MiniRecord m4(4,p,0);
    j=4;
    EXPECT_EQ(j, m4.prg_id);
    EXPECT_EQ(p, m4.path);
}

TEST_F(MiniRecordTest,equals){
    deque<Interval> v1 = {Interval(0,5)};
    deque<Interval> v2 = {Interval(1,4), Interval(15,17)};
    deque<Interval> v3 = {Interval(1,6)};
    deque<Interval> v4 = {Interval(0,3), Interval(16,18)};

    Path p;
    p.initialize(v1);
    MiniRecord m1(1,p,0);
    p.initialize(v2);
    MiniRecord m2(2,p,0);
    p.initialize(v3);
    MiniRecord m3(3,p,0);
    p.initialize(v4);
    MiniRecord m4(4,p,0);

    EXPECT_EQ(m1, m1);
    EXPECT_EQ(m2, m2);
    EXPECT_EQ(m3, m3);
    EXPECT_EQ(m4, m4);
    EXPECT_EQ((m1==m2), false);
    EXPECT_EQ((m3==m2), false);
    EXPECT_EQ((m1==m4), false);
    EXPECT_EQ((m3==m4), false);
}
