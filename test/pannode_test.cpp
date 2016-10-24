#include "gtest/gtest.h"
#include "pannode.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class PanNodeTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(PanNodeTest,create){

    PanNode pn(3);
    uint32_t j=3;
    EXPECT_EQ(j, pn.id);
    pn = PanNode(8);
    j=8;
    EXPECT_EQ(j, pn.id);
}

TEST_F(PanNodeTest,addRead){
    PanNode pn(3);
    pn.add_read(0);
    uint32_t j=1;
    EXPECT_EQ(j, pn.foundReads.size());
    j=0;
    EXPECT_EQ(j, pn.foundReads[0]);
    pn.add_read(0);
    j=2;
    EXPECT_EQ(j, pn.foundReads.size());
    j=0;
    EXPECT_EQ(j, pn.foundReads[0]);
    EXPECT_EQ(j, pn.foundReads[1]);
    pn.add_read(7);
    j=3;
    EXPECT_EQ(j, pn.foundReads.size());
    j=7;
    EXPECT_EQ(j, pn.foundReads[2]);
}

TEST_F(PanNodeTest,equals){
    PanNode pn1(3);
    PanNode pn2(2);
    PanNode pn3(2);

    EXPECT_EQ(pn1, pn1);
    EXPECT_EQ(pn2, pn2);
    EXPECT_EQ(pn3, pn3);
    EXPECT_EQ(pn2, pn3);
    EXPECT_EQ(pn3, pn2);
    EXPECT_EQ((pn1==pn2), false);
    EXPECT_EQ((pn1==pn3), false);
}

