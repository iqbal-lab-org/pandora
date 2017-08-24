#include "gtest/gtest.h"
#include "pannode.h"
#include "minihit.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class PanNodeTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};

TEST_F(PanNodeTest,create){

    PanNode pn(3, "3");
    uint32_t j=3;
    EXPECT_EQ(j, pn.id);
    EXPECT_EQ("3", pn.name);
    EXPECT_EQ((uint)1, pn.covg);
}

TEST_F(PanNodeTest,equals){
    PanNode pn1(3,"3");
    PanNode pn2(2,"2");
    PanNode pn3(2,"2");

    EXPECT_EQ(pn1, pn1);
    EXPECT_EQ(pn2, pn2);
    EXPECT_EQ(pn3, pn3);
    EXPECT_EQ(pn2, pn3);
    EXPECT_EQ(pn3, pn2);
    EXPECT_EQ((pn1==pn2), false);
    EXPECT_EQ((pn1==pn3), false);
}

TEST_F(PanNodeTest,nequals){
    PanNode pn1(3,"3");
    PanNode pn2(2,"2");
    PanNode pn3(2,"2");

    EXPECT_EQ((pn1!=pn2), true);
    EXPECT_EQ((pn2!=pn1), true);
    EXPECT_EQ((pn1!=pn1), false);
    EXPECT_EQ((pn2!=pn2), false);
    EXPECT_EQ((pn3!=pn3), false);
    EXPECT_EQ((pn2!=pn3), false);
}

TEST_F(PanNodeTest,less){
    PanNode pn1(3,"3");
    PanNode pn2(2,"2");
    PanNode pn3(2,"2");

    EXPECT_EQ((pn1<pn1), false);
    EXPECT_EQ((pn2<pn2), false);
    EXPECT_EQ((pn3<pn3), false);
    EXPECT_EQ((pn1<pn3), false);
    EXPECT_EQ((pn1<pn2), false);
    EXPECT_EQ((pn2<pn1), true);
    EXPECT_EQ((pn3<pn1), true);

}
