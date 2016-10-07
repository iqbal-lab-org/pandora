#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "interval.h"
#include "path.h"
#include "localgraph.h"
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

TEST_F(LocalGraphTest,extendPath){
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
    set<Path> s = lg3.extend_path(p);
}
