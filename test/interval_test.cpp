#include "gtest/gtest.h"
#include <stdint.h>
#include <iostream>
#include "interval.h"

using namespace std;

class IntervalTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    //     // (right before the destructor).
  }
};

TEST_F(IntervalTest,create)
{
    Interval i = Interval(0,0);
    uint32_t j = 0;
    EXPECT_EQ(i.start,j);
    EXPECT_EQ(i.end,j);
    EXPECT_EQ(i.length,j);

    i = Interval(1,9);
    j=1;
    EXPECT_EQ(i.start,j);
    j=9;
    EXPECT_EQ(i.end,j);
    j=8;
    EXPECT_EQ(i.length,j);

    // should fail if end is before start
    EXPECT_DEATH(Interval(9,1), "");
    // input should be non-negative
    EXPECT_DEATH(Interval(-1,10),"");
}

TEST_F(IntervalTest,equals)
{
    Interval i = Interval(1,5);
    Interval j = Interval(1,5);
    EXPECT_EQ(i,j);
    
    Interval k = Interval(0,4);
    EXPECT_EQ((i==k), false);

    i = Interval(0,0);
    j = Interval(1,1);
    EXPECT_EQ((i==j), false);
}
