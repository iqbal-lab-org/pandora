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
    PanRead pr1(1);
    PanRead pr2(2);
    EXPECT_NE(pr1, pr2);
    EXPECT_NE(pr2, pr1);
    EXPECT_EQ((pr1!=pr1), false);
    EXPECT_EQ((pr2!=pr2), false);
}

TEST_F(PanReadTest,less){
    PanRead pr1(1);
    PanRead pr2(2);
    EXPECT_EQ((pr1<pr1), false);
    EXPECT_EQ((pr2<pr2), false);
    EXPECT_EQ((pr1<pr2), true);
    EXPECT_EQ((pr2<pr1), false);
}

TEST_F(PanReadTest,add_hits)
{
    PanRead pr1(1);
    set<MinimizerHit*, pComp> c;

    pr1.add_hits(4, c);
    EXPECT_EQ((uint)1, pr1.hits.size());
    EXPECT_EQ((uint)0, pr1.hits[4].size());

    // add again and should append pr1.hits
    Interval i(0,5);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MinimizerHit* mh;
    mh = new MinimizerHit(4, i, 0, p, 0, 0);
    c.insert(mh);
    pr1.add_hits(4, c);
    EXPECT_EQ((uint)1, pr1.hits.size());
    EXPECT_EQ((uint)1, pr1.hits[4].size());

    // add to new id
    pr1.add_hits(5, c);
    EXPECT_EQ((uint)2, pr1.hits.size());
    EXPECT_EQ((uint)1, pr1.hits[5].size());
    EXPECT_EQ((uint)1, pr1.hits[4].size());

    delete mh;
}

TEST_F(PanReadTest,get_edge)
{
    PanRead pr1(1);
    PanEdge* e1;
    PanNode *n1, *n2;
    n1 = new PanNode(3,3,"three");
    n2 = new PanNode(6,6,"six");
    e1 = new PanEdge(n1, n2, 0);
    pr1.edges.push_back(e1);

    PanEdge* e2;
    PanNode* n3;
    n3 = new PanNode(7,7,"seven");
    e2 = new PanEdge(n2, n3, 0);
    pr1.edges.push_back(e2);

    PanEdge* e3;
    PanNode* n4;
    n4 = new PanNode(9,9,"nine");
    e3 = new PanEdge(n3, n4, 0);
    pr1.edges.push_back(e3);

    vector<PanEdge*>::iterator it = pr1.get_edge(e3);
    EXPECT_EQ(*it, e3);
    it = pr1.get_edge(e2);
    EXPECT_EQ(*it, e2);
    it = pr1.get_edge(e1);
    EXPECT_EQ(*it, e1);

    delete e1;
    delete e2;
    delete e3;
    delete n1;
    delete n2;
    delete n3;
    delete n4;
}

TEST_F(PanReadTest,get_next_edge)
{
    PanRead pr1(1);
    PanEdge* e1;
    PanNode *n1, *n2;
    n1 = new PanNode(3,3,"three");
    n2 = new PanNode(6,6,"six");
    e1 = new PanEdge(n1, n2, 0);
    pr1.edges.push_back(e1);

    PanEdge* e2;
    PanNode* n3;
    n3 = new PanNode(7,7,"seven");
    e2 = new PanEdge(n2, n3, 0);
    pr1.edges.push_back(e2);

    PanEdge* e3;
    PanNode* n4;
    n4 = new PanNode(9,9,"nine");
    e3 = new PanEdge(n3, n4, 0);
    pr1.edges.push_back(e3);
    
    vector<PanEdge*>::iterator it = pr1.get_next_edge(e3);
    EXPECT_EQ(it, pr1.edges.end());
    it = pr1.get_next_edge(e2);
    EXPECT_EQ(*it, e3);
    it = pr1.get_next_edge(e1);
    EXPECT_EQ(*it, e2);
    
    delete e1;
    delete e2;
    delete e3;
    delete n1;
    delete n2;
    delete n3;
    delete n4;
}

TEST_F(PanReadTest,get_previous_edge)
{
    PanRead pr1(1);
    PanEdge* e1;
    PanNode *n1, *n2;
    n1 = new PanNode(3,3,"three");
    n2 = new PanNode(6,6,"six");
    e1 = new PanEdge(n1, n2, 0);
    pr1.edges.push_back(e1);

    PanEdge* e2;
    PanNode* n3;
    n3 = new PanNode(7,7,"seven");
    e2 = new PanEdge(n2, n3, 0);
    pr1.edges.push_back(e2);

    PanEdge* e3;
    PanNode* n4;
    n4 = new PanNode(9,9,"nine");
    e3 = new PanEdge(n3, n4, 0);
    pr1.edges.push_back(e3);

    vector<PanEdge*>::iterator it = pr1.get_previous_edge(e1);
    EXPECT_EQ(it, pr1.edges.end());
    it = pr1.get_previous_edge(e2);
    EXPECT_EQ(*it, e1);
    it = pr1.get_previous_edge(e3);
    EXPECT_EQ(*it, e2);

    delete e1;
    delete e2;
    delete e3;
    delete n1;
    delete n2;
    delete n3;
    delete n4;
}

TEST_F(PanReadTest,get_other_edge)
{
    PanRead pr1(1);
    PanEdge* e1;
    PanNode *n1, *n2;
    n1 = new PanNode(3,3,"three");
    n2 = new PanNode(6,6,"six");
    e1 = new PanEdge(n1, n2, 0);
    pr1.edges.push_back(e1);

    PanEdge* e2;
    PanNode* n3;
    n3 = new PanNode(7,7,"seven");
    e2 = new PanEdge(n2, n3, 0);
    pr1.edges.push_back(e2);

    PanEdge* e3;
    PanNode* n4;
    n4 = new PanNode(9,9,"nine");
    e3 = new PanEdge(n3, n4, 0);
    pr1.edges.push_back(e3);

    vector<PanEdge*>::iterator it = pr1.get_other_edge(e3, n3);
    EXPECT_EQ(*it, e2);
    it = pr1.get_other_edge(e2, n3);
    EXPECT_EQ(*it, e3);
    it = pr1.get_other_edge(e2, n2);
    EXPECT_EQ(*it, e1);
    it = pr1.get_other_edge(e1, n2);
    EXPECT_EQ(*it, e2);
    it = pr1.get_other_edge(e1, n1);
    EXPECT_EQ(it, pr1.edges.end());
    it = pr1.get_other_edge(e3, n4);
    EXPECT_EQ(it, pr1.edges.end());

    delete e1;
    delete e2;
    delete e3;
    delete n1;
    delete n2;
    delete n3;
    delete n4;
}

TEST_F(PanReadTest,replace_edge1)
{
}

TEST_F(PanReadTest,replace_edge2)
{
}

TEST_F(PanReadTest,remove_edge1)
{
}

TEST_F(PanReadTest,remove_edge2)
{
}

TEST_F(PanReadTest,replace_node)
{
}

TEST_F(PanReadTest,remove_node)
{
}
