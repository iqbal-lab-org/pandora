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
    Interval i = Interval(1,9);
    EXPECT_EQ(i.start,1);
    EXPECT_EQ(i.end,9);
    EXPECT_EQ(i.length,8);

    // should fail if end is before start
    EXPECT_DEATH(Interval(9,1), "");
    // input should be non-negative
    EXPECT_DEATH(Interval(-1,10),"");
}

TEST_F(IntervalTest,compare)
{
    Interval i = Interval(1,5);
    Interval j = Interval(1,5);
    EXPECT_EQ(i,j);
    
    Interval k = Interval(0,4);
    EXPECT_EQ(false, (i==k));
}
