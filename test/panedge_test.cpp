#include "gtest/gtest.h"
#include "pannode.h"
#include "panedge.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class PanEdgeTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};

TEST_F(PanEdgeTest,create){

    PanNode *pn1, *pn2;
    pn1 = new PanNode(4,3, "three");
    pn2 = new PanNode(7,6, "six");
    PanEdge pe(pn1, pn2, 0);

    EXPECT_EQ((uint)3, pe.from->node_id);
    EXPECT_EQ((uint)6, pe.to->node_id);
    EXPECT_EQ((uint)0, pe.orientation);
    EXPECT_EQ((uint)1, pe.covg);

    EXPECT_DEATH(PanEdge(pn1, pn2, 5), "");
    delete pn1;
    delete pn2;
}

TEST_F(PanEdgeTest,equals){
    PanNode *pn1, *pn2, *pn3;
    pn1 = new PanNode(3,3, "three");
    pn2 = new PanNode(6,6, "six");
    pn3 = new PanNode(9,9, "nine");
    PanEdge pe1(pn1, pn2, 0);
    PanEdge pe2(pn1, pn2, 2);
    PanEdge pe3(pn2, pn1, 0);
    PanEdge pe4(pn1, pn3, 0);
    PanEdge pe5(pn1, pn2, 0);
    PanEdge pe6(pn2, pn1, 3);
    PanEdge pe7(pn2, pn1, 2);
    PanEdge pe8(pn1, pn3, 1);
    PanEdge pe9(pn3, pn1, 1);

    // identity
    EXPECT_EQ(pe1, pe1);
    EXPECT_EQ(pe2, pe2);
    EXPECT_EQ(pe3, pe3);
    EXPECT_EQ(pe4, pe4);
    EXPECT_EQ(pe5, pe5);
    // a separate node defined identically
    EXPECT_EQ(pe1, pe5);
    EXPECT_EQ(pe5, pe1);
    // complement edges
    EXPECT_EQ(pe1, pe6);
    EXPECT_EQ(pe6, pe1);
    EXPECT_EQ(pe2, pe7);
    EXPECT_EQ(pe7, pe2);
    EXPECT_EQ(pe8, pe9);
    EXPECT_EQ(pe9, pe8);
    // not equal
    EXPECT_EQ((pe1==pe2), false);
    EXPECT_EQ((pe1==pe3), false);
    EXPECT_EQ((pe1==pe4), false);
    EXPECT_EQ((pe2==pe1), false);
    EXPECT_EQ((pe3==pe1), false);
    EXPECT_EQ((pe4==pe1), false);

    delete pn1;
    delete pn2;
    delete pn3;
}

TEST_F(PanEdgeTest,nequals){
    PanNode *pn1, *pn2, *pn3;
    pn1 = new PanNode(3,3, "three");
    pn2 = new PanNode(6,6, "six");
    pn3 = new PanNode(9,9, "nine");
    PanEdge pe1(pn1, pn2, 0);
    PanEdge pe2(pn1, pn2, 2);
    PanEdge pe3(pn2, pn1, 0);
    PanEdge pe4(pn1, pn3, 0);
    PanEdge pe5(pn1, pn2, 0);
    PanEdge pe6(pn2, pn1, 3);
    
    EXPECT_EQ((pe1!=pe1), false);
    EXPECT_EQ((pe2!=pe2), false);
    EXPECT_EQ((pe3!=pe3), false);
    EXPECT_EQ((pe4!=pe4), false);
    EXPECT_EQ((pe5!=pe5), false);
    EXPECT_EQ((pe6!=pe6), false);
    EXPECT_EQ((pe1!=pe5), false);
    EXPECT_EQ((pe5!=pe1), false);
    EXPECT_EQ((pe1!=pe6), false);
    EXPECT_EQ((pe6!=pe1), false);
    EXPECT_EQ((pe1!=pe2), true);
    EXPECT_EQ((pe1!=pe3), true);
    EXPECT_EQ((pe1!=pe4), true);
    EXPECT_EQ((pe2!=pe1), true);
    EXPECT_EQ((pe3!=pe1), true);
    EXPECT_EQ((pe4!=pe1), true);
    
    delete pn1;
    delete pn2;
    delete pn3;
}

TEST_F(PanEdgeTest,less){
    PanNode *pn1, *pn2, *pn3;
    pn1 = new PanNode(3,3, "three");
    pn2 = new PanNode(6,6, "six");
    pn3 = new PanNode(9,9, "nine");
    PanEdge pe1(pn1, pn2, 0);
    PanEdge pe2(pn1, pn2, 2);
    PanEdge pe3(pn2, pn1, 0);
    PanEdge pe4(pn1, pn3, 0);
    PanEdge pe5(pn1, pn2, 0);
    
    EXPECT_EQ((pe1<pe1), false);
    EXPECT_EQ((pe2<pe2), false);
    EXPECT_EQ((pe3<pe3), false);
    EXPECT_EQ((pe4<pe4), false);
    EXPECT_EQ((pe5<pe5), false);
    EXPECT_EQ((pe1<pe5), false);
    EXPECT_EQ((pe5<pe1), false);
    EXPECT_EQ((pe1<pe2), true);
    EXPECT_EQ((pe1<pe3), true);
    EXPECT_EQ((pe1<pe4), true);
    EXPECT_EQ((pe2<pe1), false);
    EXPECT_EQ((pe3<pe1), false);
    EXPECT_EQ((pe4<pe1), false);

    delete pn1;
    delete pn2;
    delete pn3;
}
