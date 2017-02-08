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
    Interval i(0,0);
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

TEST_F(IntervalTest, write)
{
    Interval i(1,5);
    stringstream out;
    out << i;
    EXPECT_EQ(out.str(), "[1, 5)");
}

TEST_F(IntervalTest, read)
{
    Interval i(1,5);
    stringstream out;
    out << i;
    Interval j;
    out >> j;
    EXPECT_EQ(i,j);
}

TEST_F(IntervalTest,equals)
{
    Interval i(1,5);
    Interval j(1,5);
    EXPECT_EQ(i,j);
    EXPECT_EQ(j,i);
    
    Interval k(0,4);
    EXPECT_EQ((i==k), false);
    EXPECT_EQ((k==j), false);

    i = Interval(0,0);
    j = Interval(1,1);
    EXPECT_EQ((i==i), true);
    EXPECT_EQ((j==j), true);
    EXPECT_EQ((i==j), false);
    EXPECT_EQ(i,i);
    EXPECT_EQ(j,j);
}

TEST_F(IntervalTest,notequals)
{
    Interval i(1,5);
    Interval j(1,5);
    EXPECT_EQ((i!=j), false);

    Interval k(0,4);
    EXPECT_EQ((i!=k), true);
    EXPECT_EQ((j!=k), true);
    EXPECT_NE(i,k);
    EXPECT_NE(k,j);

    i = Interval(0,0);
    j = Interval(1,1);
    EXPECT_EQ((i!=i), false);
    EXPECT_EQ((j!=j), false);
    EXPECT_EQ((i!=j), true);
    EXPECT_NE(i,j);
    EXPECT_NE(j,i);
}

TEST_F(IntervalTest,lessthan)
{
    Interval i(1,5);
    Interval j(2,5);
    Interval k(0,4);
    Interval l(0,7);

    EXPECT_EQ((i<i), false);

    EXPECT_EQ((i<j), true);
    EXPECT_EQ((j<i), false);

    EXPECT_EQ((i<k), false);
    EXPECT_EQ((k<i), true);
    EXPECT_EQ((k<j), true);
    EXPECT_EQ((j<k), false);

    EXPECT_EQ((k<l), true);
    EXPECT_EQ((l<i), true);
    EXPECT_EQ((l<j), true);
    EXPECT_EQ((l<k), false);
    EXPECT_EQ((i<l), false);
    EXPECT_EQ((j<l), false);
}
