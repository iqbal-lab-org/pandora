#include "gtest/gtest.h"
#include "localnode.h"
#include "interval.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class LocalNodeTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(LocalNodeTest,create){

    LocalNode ln("hello", Interval(0,5), 0);

    uint32_t j = 1;
    EXPECT_EQ("hello", ln.seq);
    EXPECT_EQ(Interval(0,5), ln.pos);
    j=0;
    EXPECT_EQ(j, ln.id);
}

TEST_F(LocalNodeTest,equals){
    LocalNode ln1("hello", Interval(0,5), 0);
    LocalNode ln2("heppo", Interval(0,5), 0);
    LocalNode ln3("hello", Interval(0,4), 0);
    LocalNode ln4("hello", Interval(0,5), 1);
    // can't compare outNodes bit outside of localGraph
    EXPECT_EQ(ln1, ln1);
    EXPECT_EQ(ln2, ln2);
    EXPECT_EQ(ln3, ln3);
    EXPECT_EQ(ln4, ln4);
    EXPECT_EQ((ln1==ln2), false);
    EXPECT_EQ((ln1==ln3), false);
    EXPECT_EQ((ln1==ln4), false);
    EXPECT_EQ((ln2==ln3), false);
    EXPECT_EQ((ln2==ln4), false);
    EXPECT_EQ((ln3==ln4), false);
}

TEST_F(LocalNodeTest,compare){
}

