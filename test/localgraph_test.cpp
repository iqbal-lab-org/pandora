#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "interval.h"
#include "path.h"
#include "localgraph.h"
#include "localnode.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class LocalGraphTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }

/*  LocalGraph lg0 = LocalGraph();
  lg0.add_node(0,"",Interval(0,0));

  LocalGraph lg1 = LocalGraph();
  lg1.add_node(0,"AGCT", Interval(0,4));

  LocalGraph lg2 = LocalGraph();
  lg2.add_node(0,"A", Interval(0,1));
  lg2.add_node(1,"GC", Interval(4,6));
  lg2.add_node(2,"G", Interval(7,8));
  lg2.add_node(3,"T", Interval(13,14));
  lg2.add_edge(0,1);
  lg2.add_edge(0,2);
  lg2.add_edge(1,3);
  lg2.add_edge(2,3);

  LocalGraph lg3 = LocalGraph();
  lg3.add_node(0,"A", Interval(0,1));
  lg3.add_node(1,"G", Interval(4,5));
  lg3.add_node(2,"C", Interval(8,9));
  lg3.add_node(3,"T", Interval(12,13));
  lg3.add_node(4,"", Interval(16,16));
  lg3.add_node(5,"G", Interval(19,20));
  lg3.add_node(6,"T", Interval(23,24));
  lg3.add_edge(0,1);
  lg3.add_edge(0,5);
  lg3.add_edge(1,2);
  lg3.add_edge(1,3);
  lg3.add_edge(2,4);
  lg3.add_edge(3,4);
  lg3.add_edge(4,6);
  lg3.add_edge(5,6);
*/
};

TEST_F(LocalGraphTest, addNode)
{
    // add node and check it's there
    LocalGraph lg1 = LocalGraph();
    lg1.add_node(0,"AGCT", Interval(0,4));
    LocalNode ln1 = LocalNode("AGCT", Interval(0,4), 0);
    EXPECT_EQ(ln1, *lg1.nodes[0]);

    // add node another time and expect nothing to happen
    lg1.add_node(0,"AGCT", Interval(0,4));
    EXPECT_EQ(ln1, *lg1.nodes[0]);

    // add impossible nodes and expect and error
    EXPECT_DEATH(lg1.add_node(0,"AGGT", Interval(0,4)), "");
    EXPECT_DEATH(lg1.add_node(1,"AGG", Interval(0,4)), "");
}

TEST_F(LocalGraphTest, addEdge)
{
    LocalGraph lg2 = LocalGraph();
    lg2.add_node(0,"A", Interval(0,1));
    lg2.add_node(1,"GC", Interval(4,6));
    lg2.add_node(2,"G", Interval(7,8));
    lg2.add_node(3,"T", Interval(13,14));
    lg2.add_edge(0,1);
    EXPECT_EQ(lg2.nodes[0]->outNodes[0], lg2.nodes[1]);
    lg2.add_edge(0,2);
    lg2.add_edge(1,3);
    lg2.add_edge(2,3);

    // expect failure if a node doesn't exist in the graph
    EXPECT_DEATH(lg2.add_edge(0,4),"");
}

TEST_F(LocalGraphTest, equals)
{
    LocalGraph lg1 = LocalGraph();
    lg1.add_node(0,"AGCT", Interval(0,4));
    EXPECT_EQ(lg1, lg1);

    LocalGraph lg2 = LocalGraph();
    lg2.add_node(0,"A", Interval(0,1));
    lg2.add_node(1,"GC", Interval(4,6));
    lg2.add_node(2,"G", Interval(7,8));
    lg2.add_node(3,"T", Interval(13,14));
    lg2.add_edge(0,1);
    lg2.add_edge(0,2);
    lg2.add_edge(1,3);
    lg2.add_edge(2,3); 
    EXPECT_EQ(lg2, lg2);

    EXPECT_EQ((lg1==lg2), false);

    // order adding shouldn't matter
    LocalGraph lg2p = LocalGraph();
    lg2p.add_node(2,"G", Interval(7,8));
    lg2p.add_node(3,"T", Interval(13,14));
    lg2p.add_node(1,"GC", Interval(4,6));
    lg2p.add_node(0,"A", Interval(0,1));
    lg2p.add_edge(1,3);
    lg2p.add_edge(2,3);
    lg2p.add_edge(0,1);
    lg2p.add_edge(0,2);
    EXPECT_EQ(lg2, lg2p);
}

/*TEST_F(LocalGraphTest, walk)
{
}*/

/*TEST_F(LocalGraphTest,extendPath){
    deque<Interval> d = {Interval(0,1)};
    Path p;
    p.initialize(d);
    LocalGraph lg3;
    lg3.add_node(0,"A", Interval(0,1));
    lg3.add_node(1,"G", Interval(4,5));
    lg3.add_node(2,"C", Interval(8,9));
    lg3.add_node(3,"T", Interval(12,13));
    lg3.add_node(4,"", Interval(16,16));
    lg3.add_node(5,"G", Interval(19,20));
    lg3.add_node(6,"T", Interval(23,24));
    lg3.add_edge(0,1);
    lg3.add_edge(0,5);
    lg3.add_edge(1,2);
    lg3.add_edge(1,3);
    lg3.add_edge(2,4);
    lg3.add_edge(3,4);
    lg3.add_edge(4,6);
    lg3.add_edge(5,6);

    // extend path over branch point
    set<Path> s1 = lg3.extend_path(p);
    set<Path> s_prime;
    deque<Interval> d_prime = {Interval(0,1), Interval(4,5)};
    Path p_prime;
    p_prime.initialize(d_prime);
    s_prime.insert(p_prime);
    d_prime = {Interval(0,1), Interval(19,20)};
    p_prime.initialize(d_prime);
    s_prime.insert(p_prime);
    EXPECT_ITERABLE_EQ(set<Path>, s_prime, s1);

    //extend path again
    set<Path>::iterator it = s1.begin();
    set<Path> s2 = lg3.extend_path(*it);
    s_prime.clear();
    d_prime = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p_prime.initialize(d_prime);
    s_prime.insert(p_prime);
    d_prime = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p_prime.initialize(d_prime);
    s_prime.insert(p_prime);
    EXPECT_ITERABLE_EQ(set<Path>, s_prime, s2);

    // extend path no branch point
    it++;
    s2 = lg3.extend_path(*it);
    s_prime.clear();
    d_prime = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p_prime.initialize(d_prime);
    s_prime.insert(p_prime);
    EXPECT_ITERABLE_EQ(set<Path>, s_prime, s2);

    // extend path at end, should return empty set
    s2 = lg3.extend_path(*s2.begin());
    s_prime.clear();
    EXPECT_ITERABLE_EQ(set<Path>, s_prime, s2);
}*/


