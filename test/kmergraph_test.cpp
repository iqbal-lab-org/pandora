#include "gtest/gtest.h"
#include "test_macro.cpp"

#include "interval.h"
#include "path.h"
#include "kmergraph.h"
#include "kmernode.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class KmerGraphTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }

};

TEST_F(KmerGraphTest, addNode)
{
    // add node and check it's there
    KmerGraph kg;

    deque<Interval> d = {Interval(0,3)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    uint j = 1;
    EXPECT_EQ(j, kg.nodes.size());
    EXPECT_EQ(p, kg.nodes[0]->path);
    j = 0;
    EXPECT_EQ(j, kg.nodes[0]->id);
    j = 2;
    EXPECT_EQ(j, kg.nodes[0]->covg.size());
    j = 0;
    EXPECT_EQ(j, kg.nodes[0]->covg[0]);
    EXPECT_EQ(j, kg.nodes[0]->covg[0]);

    // add node another time and expect nothing to happen
    kg.add_node(p);
    j = 1;
    EXPECT_EQ(j, kg.nodes.size());
    EXPECT_EQ(p, kg.nodes[0]->path);
    j = 0;
    EXPECT_EQ(j, kg.nodes[0]->id);
    j = 2;
    EXPECT_EQ(j, kg.nodes[0]->covg.size());
    j = 0;
    EXPECT_EQ(j, kg.nodes[0]->covg[0]);
    EXPECT_EQ(j, kg.nodes[0]->covg[0]);

    // add a second node and check gets next id
    d = {Interval(1,4)};
    p.initialize(d);
    kg.add_node(p);
    j = 2;
    EXPECT_EQ(j, kg.nodes.size());
    EXPECT_EQ(p, kg.nodes[1]->path);
    j = 0;
    EXPECT_EQ(j, kg.nodes[0]->id);
    EXPECT_EQ(j, kg.nodes[0]->covg[0]);
    EXPECT_EQ(j, kg.nodes[1]->covg[0]);
    j = 1;
    EXPECT_EQ(j, kg.nodes[1]->id);
}

TEST_F(KmerGraphTest, addEdge)
{
    // add edge and check it's there
    KmerGraph kg;

    deque<Interval> d = {Interval(0,3)};
    Path p1,p2, p3;
    p1.initialize(d);
    kg.add_node(p1);
    d = {Interval(1,4)};
    p2.initialize(d);
    kg.add_node(p2);
    uint j = 2;
    EXPECT_EQ(j, kg.nodes.size());

    // first via path constructor
    kg.add_edge(p1, p2);
    j = 1;
    EXPECT_EQ(j, kg.nodes[0]->outNodes.size());
    EXPECT_EQ(j, kg.nodes[1]->inNodes.size());
    j = 0;
    EXPECT_EQ(j, kg.nodes[1]->outNodes.size());
    EXPECT_EQ(j, kg.nodes[0]->inNodes.size());

    // repeat with path constructor, nothing should happen
    kg.add_edge(p1, p2);
    j = 1;
    EXPECT_EQ(j, kg.nodes[0]->outNodes.size());
    EXPECT_EQ(j, kg.nodes[1]->inNodes.size());
    j = 0;
    EXPECT_EQ(j, kg.nodes[1]->outNodes.size());
    EXPECT_EQ(j, kg.nodes[0]->inNodes.size());

    // expect failure if a node doesn't exist in the graph
    d = {Interval(4,7)};
    p3.initialize(d);
    EXPECT_DEATH(kg.add_edge(p1,p3),"");
    EXPECT_DEATH(kg.add_edge(p3,p2),"");

    /*//now test the same with id edge constructor
    KmerGraph kg2;
    kg2.add_node(p1);
    kg2.add_node(p2);
    j = 2;
    EXPECT_EQ(j, kg2.nodes.size());
    kg2.add_edge(0, 1);
    j = 1;
    EXPECT_EQ(j, kg2.nodes[0]->outNodes.size());
    EXPECT_EQ(j, kg2.nodes[1]->inNodes.size());
    j = 0;
    EXPECT_EQ(j, kg2.nodes[1]->outNodes.size());
    EXPECT_EQ(j, kg2.nodes[0]->inNodes.size());

    // repeat with path constructor, nothing should happen
    kg2.add_edge(0, 1);
    j = 1;
    EXPECT_EQ(j, kg2.nodes[0]->outNodes.size());
    EXPECT_EQ(j, kg2.nodes[1]->inNodes.size());
    j = 0;
    EXPECT_EQ(j, kg2.nodes[1]->outNodes.size());
    EXPECT_EQ(j, kg2.nodes[0]->inNodes.size());*/

}

TEST_F(KmerGraphTest, clear)
{
    KmerGraph kg;
    deque<Interval> d = {Interval(0,3)};
    Path p1,p2;
    p1.initialize(d);
    kg.add_node(p1);
    d = {Interval(1,4)};
    p2.initialize(d);
    kg.add_node(p2);
    kg.add_edge(p1,p2);
    uint j = 2;
    EXPECT_EQ(j, kg.nodes.size());

    kg.clear();
    j = 0;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_node(p1);
    kg.add_node(p2);
    kg.add_edge(p1,p2);
    j = 2;
    EXPECT_EQ(j, kg.nodes.size());
}

TEST_F(KmerGraphTest, equals)
{
    KmerGraph kg1, kg2;
    deque<Interval> d = {Interval(0,3)};
    Path p1,p2,p3;
    p1.initialize(d);
    kg1.add_node(p1);
    kg2.add_node(p1);
    d = {Interval(1,4)};
    p2.initialize(d);
    kg1.add_node(p2);
    kg2.add_node(p2);
    kg1.add_edge(p1,p2);
    kg2.add_edge(p1,p2);

    d = {Interval(2,5)};
    p3.initialize(d);
    kg2.add_node(p3);

    // same as themselves, different if different numbers of nodes
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1==kg2), false);
    EXPECT_EQ((kg2==kg1), false);

    kg1.add_node(p3);    
    kg2.add_edge(p1, p3);

    // same as themselves, different if different numbers of edges
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1==kg2), false);
    EXPECT_EQ((kg2==kg1), false);

    kg1.add_edge(p2, p3);

    // same as themselves, different if edges in different places
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1==kg2), false);
    EXPECT_EQ((kg2==kg1), false);
}

TEST_F(KmerGraphTest,findMaxPathSimple)
{
    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);    
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 7;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0],kg.nodes[1]);
    kg.add_edge(kg.nodes[1],kg.nodes[2]);
    kg.add_edge(kg.nodes[0],kg.nodes[3]);
    kg.add_edge(kg.nodes[3],kg.nodes[4]);
    kg.add_edge(kg.nodes[0],kg.nodes[5]);
    kg.add_edge(kg.nodes[2],kg.nodes[6]);
    kg.add_edge(kg.nodes[4],kg.nodes[6]);
    kg.add_edge(kg.nodes[5],kg.nodes[6]);

    kg.nodes[1]->covg[0]+=4;
    kg.nodes[2]->covg[0]+=3;
    kg.nodes[5]->covg[1]+=5;

    kg.num_reads = 5;
    kg.k = 3;

    vector<KmerNode*> mp;
    kg.find_max_path(0, 0.01, mp);
    vector<KmerNode*> exp_order = {kg.nodes[1], kg.nodes[2]};    
    EXPECT_ITERABLE_EQ(vector<KmerNode*>, exp_order, mp);

    mp.clear();
    kg.find_max_path(1, 0.01, mp);
    exp_order = {kg.nodes[5]};
    EXPECT_ITERABLE_EQ(vector<KmerNode*>, exp_order, mp);
}

TEST_F(KmerGraphTest,findMaxPath2Level)
{
    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(8, 9), Interval(16,18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(12,13), Interval(16,18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(16,18), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 10;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0],kg.nodes[1]);
    kg.add_edge(kg.nodes[1],kg.nodes[2]);
    kg.add_edge(kg.nodes[2],kg.nodes[3]);
    kg.add_edge(kg.nodes[0],kg.nodes[4]);
    kg.add_edge(kg.nodes[4],kg.nodes[5]);
    kg.add_edge(kg.nodes[5],kg.nodes[6]);
    kg.add_edge(kg.nodes[3],kg.nodes[7]);
    kg.add_edge(kg.nodes[6],kg.nodes[7]);
    kg.add_edge(kg.nodes[0],kg.nodes[8]);
    kg.add_edge(kg.nodes[7],kg.nodes[9]);
    kg.add_edge(kg.nodes[8],kg.nodes[9]);

    kg.nodes[4]->covg[0]+=4;
    kg.nodes[5]->covg[0]+=3;
    kg.nodes[6]->covg[0]+=5;
    kg.nodes[7]->covg[0]+=4;
    kg.nodes[8]->covg[1]+=5;

    kg.num_reads = 5;
    kg.k = 3;

    vector<KmerNode*> mp;
    kg.find_max_path(0, 0.01, mp);
    vector<KmerNode*> exp_order = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNode*>, exp_order, mp);

    mp.clear();
    kg.find_max_path(1, 0.01, mp);
    exp_order = {kg.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNode*>, exp_order, mp);
}

TEST_F(KmerGraphTest, save){
    KmerGraph kg;
    deque<Interval> d = {Interval(0,3)};
    Path p1,p2;
    p1.initialize(d);
    kg.add_node(p1);
    d = {Interval(1,4)};
    p2.initialize(d);
    kg.add_node(p2);
    kg.add_edge(p1,p2);
    kg.nodes[0]->covg[1] +=5;

    kg.save("../test/test_cases/kmergraph_test.gfa");
}

TEST_F(KmerGraphTest, load){
    KmerGraph kg, read_kg;
    deque<Interval> d = {Interval(0,3)};
    Path p1,p2;
    p1.initialize(d);
    kg.add_node(p1);
    d = {Interval(1,4)};
    p2.initialize(d);
    kg.add_node(p2);
    kg.add_edge(p1,p2);
    kg.nodes[0]->covg[1] +=5;

    read_kg.load("../test/test_cases/kmergraph_test.gfa");
    EXPECT_EQ(kg, read_kg);
}

