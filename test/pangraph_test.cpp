#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "pangraph.h"
#include "pannode.h"
#include "minihit.h"
#include <stdint.h>
#include <cassert>
#include <iostream>

using namespace std;

class PanGraphTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(PanGraphTest, addNode)
{
    set<MinimizerHit*, pComp> mhs;

    // add node and check it's there
    PanGraph pg;
    pg.add_node(0,"0",0, mhs);

    PanNode *pn;
    pn = new PanNode(0,"0");
    pn->add_read(0);
    EXPECT_EQ(*pg.nodes[0], *pn);
    uint32_t j = 1;
    EXPECT_EQ(pg.nodes.size(), j);
    EXPECT_EQ(pg.nodes[0]->foundReads.size(), j);
    j = 0;
    EXPECT_EQ(pg.nodes[0]->foundHits.size(), j);

    // add node again with same read
    pg.add_node(0,"0",0, mhs);
    pn->add_read(0);
    EXPECT_EQ(*pg.nodes[0], *pn);
    j = 1;
    EXPECT_EQ(pg.nodes.size(), j);
    j = 2;
    EXPECT_EQ(pg.nodes[0]->foundReads.size(), j);

    // add node again with different read
    pg.add_node(0,"0",2, mhs);
    pn->add_read(2);
    EXPECT_EQ(*pg.nodes[0], *pn);
    j = 1;
    EXPECT_EQ(pg.nodes.size(), j);
    j = 3;
    EXPECT_EQ(pg.nodes[0]->foundReads.size(), j);
    delete pn;

    // add different node
    pg.add_node(1,"1",2, mhs);
    pn = new PanNode(1,"1");
    pn->add_read(2);
    EXPECT_EQ(*pg.nodes[1], *pn);
    j = 2;
    EXPECT_EQ(pg.nodes.size(), j);
    j = 1;
    EXPECT_EQ(pg.nodes[1]->foundReads.size(), j);
    delete pn;

    // add a node with hits
    Path p;
    deque<Interval> d = {Interval(0,1), Interval(4,7)};
    p.initialize(d);
    MinimizerHit* mh0;
    mh0 = new MinimizerHit(2, Interval(1,5), 2, p, true);
    mhs.insert(mh0);
    d = {Interval(0,1), Interval(5,8)};
    p.initialize(d);
    MinimizerHit* mh1;
    mh1 = new MinimizerHit(2, Interval(1,5), 2, p, true);
    mhs.insert(mh1);
    pg.add_node(2,"2",2, mhs);
    pn = new PanNode(2,"2");
    pn->add_read(2);
    EXPECT_EQ(*pg.nodes[2], *pn);
    j = 3;
    EXPECT_EQ(pg.nodes.size(), j);
    j = 2;
    EXPECT_EQ(pg.nodes[2]->foundHits.size(), j);
    j = 1;
    EXPECT_EQ(pg.nodes[2]->foundReads.size(), j);
    MinimizerHit* mh2;
    mh2 = new MinimizerHit(0, Interval(1,5), 0, p, true);
    mhs.insert(mh2);
    //pg.add_node(0,"0",0, mhs);
    EXPECT_DEATH(pg.add_node(0,"0",0, mhs), "");
    delete pn;
    delete mh0;
    delete mh1;
    delete mh2;
}

TEST_F(PanGraphTest, addEdge)
{
    set<MinimizerHit*, pComp> mhs;
    PanGraph pg;
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",2, mhs);
    pg.add_edge(0,1);

    PanNode *pn1;
    pn1 = new PanNode(0,"0");
    pn1->add_read(0);
    PanNode *pn2;
    pn2 = new PanNode(1,"1");
    pn2->add_read(2);
    pn1->outNodes.push_back(pn2);

    EXPECT_EQ(*pg.nodes[0], *pn1);
    EXPECT_EQ(*pg.nodes[1], *pn2);

    EXPECT_EQ(*(pg.nodes[0]->outNodes[0]), *pn2);

    // expect failure if a node doesn't exist in the graph
    EXPECT_DEATH(pg.add_edge(0,4),"");
    delete pn1;
    delete pn2;
}

TEST_F(PanGraphTest, equals)
{
    set<MinimizerHit*, pComp> mhs;
    PanGraph pg1;
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_node(2,"2",2, mhs);
    pg1.add_edge(0,1);
    pg1.add_edge(1,2);
  
    PanGraph pg2;
    pg2.add_node(1,"1",2, mhs);
    pg2.add_node(0,"0",0, mhs);
    pg2.add_edge(0,1);
    pg2.add_node(2,"2",2, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(1,2);

    // adding nodes and edges in different order should make no difference
    EXPECT_EQ(pg1, pg1);
    EXPECT_EQ(pg2, pg2);
    EXPECT_EQ(pg1, pg2);
    EXPECT_EQ(pg2, pg1);

    // adding an extra edge does make a difference
    pg2.add_edge(0,2);
    EXPECT_EQ((pg1 == pg2), false);

    // having one fewer edge makes a difference
    PanGraph pg3;
    pg3.add_node(1,"1",2, mhs);
    pg3.add_node(0,"0",0, mhs);
    pg3.add_node(2,"2",2, mhs);
    pg3.add_node(1,"1",0, mhs);
    pg3.add_edge(1,2);
    EXPECT_EQ((pg1 == pg3), false);

    // or one extra node
    pg3.add_edge(0,1);
    EXPECT_EQ((pg1 == pg3), true); //adds the missing edge
    pg3.add_node(3,"3",0, mhs);
    EXPECT_EQ((pg1 == pg3), false);

}

TEST_F(PanGraphTest, write_gfa)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg2;
    pg2.add_node(1,"1",2, mhs);
    pg2.add_node(0,"0",0, mhs);
    pg2.add_edge(0,1);
    pg2.add_node(2,"2",2, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(1,2);
    pg2.add_edge(1,2);
    pg2.write_gfa("../test/test_cases/pangraph_test_save.gfa");
}
