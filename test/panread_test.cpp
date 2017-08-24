#include "gtest/gtest.h"
#include "pannode.h"
#include "panedge.h"
#include "panread.h"
#include "minihit.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class PanReadTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};

TEST_F(PanReadTest,create){

    PanRead pr(3);
    EXPECT_EQ((uint)3, pr.id);
    EXPECT_EQ((uint)0, pr.edges.size());
    EXPECT_EQ((uint)0, pr.hits.size());
}

TEST_F(PanReadTest,equals){
    PanRead pr1(1);
    PanRead pr2(2);
    EXPECT_EQ(pr1, pr1);
    EXPECT_EQ(pr2, pr2);
    EXPECT_EQ((pr1==pr2), false);
    EXPECT_EQ((pr2==pr1), false);   
}

/*TEST_F(PanReadTest,equals){
    PanNode *pn1, *pn2, *pn3;
    PanEdge *pe1, *pe2, *pe3, *pe4, *pe5;
    pn1 = new PanNode(3, "three");
    pn2 = new PanNode(6, "six");
    pn3 = new PanNode(9, "nine");
    pe1 = new PanEdge(pn1, pn2, 0);
    pe2 = new PanEdge(pn1, pn2, 2);
    pe3 = new PanEdge(pn2, pn3, 0);
    pe4 = new PanEdge(pn3, pn1, 0);
    pe5 = new PanEdge(pn2, pn1, 3);

    PanRead pr1(1);
    pr1.edges = {pe1, pe3};

    // less and edge
    PanRead pr2(2);
    pr2.edges = {pe1};

    // add an edge
    PanRead pr3(3);
    pr3.edges = {pe1, pe3, pe4};

    // with a different edge
    PanRead pr4(4);
    pr4.edges = {pe2, pe3}; // NB not valid combination of orientations...

    // same edges, but one written rev 
    PanRead pr5(5);
    pr5.edges = {pe5, pe3};

    EXPECT_EQ(pr1, pr1);
    EXPECT_EQ(pr2, pr2);
    EXPECT_EQ(pr3, pr3);
    EXPECT_EQ(pr4, pr4);
    EXPECT_EQ(pr5, pr5);
    EXPECT_EQ(pr1, pr5);
    EXPECT_EQ(pr5, pr1);
    EXPECT_EQ((pr1==pr2), false);
    EXPECT_EQ((pr1==pr3), false);
    EXPECT_EQ((pr1==pr4), false);
    EXPECT_EQ((pr2==pr1), false);
    EXPECT_EQ((pr3==pr1), false);
    EXPECT_EQ((pr4==pr1), false);
    EXPECT_EQ((pr2==pr3), false);
    EXPECT_EQ((pr2==pr4), false);
    EXPECT_EQ((pr3==pr4), false);
   
    delete pn1;
    delete pn2;
    delete pn3;
    delete pe1;
    delete pe2;
    delete pe3;
    delete pe4;
    delete pe5;
}*/

TEST_F(PanReadTest,nequals){
}

TEST_F(PanReadTest,less){
}
