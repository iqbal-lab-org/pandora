#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "pangraph.h"
#include "pannode.h"
#include "panedge.h"
#include "panread.h"
#include "pansample.h"
#include "minihit.h"
#include "localPRG.h"
#include <stdint.h>
#include <numeric>
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

TEST_F(PanGraphTest, add_node)
{
    set<MinimizerHit*, pComp> mhs;

    // add node and check it's there
    PanGraph pg;
    pg.add_node(0,"0",1, mhs);

    PanNode *pn;
    pn = new PanNode(0,0,"0");
    uint32_t j = 1;
    EXPECT_EQ(pg.nodes.size(), j);
    EXPECT_EQ(*pg.nodes[0], *pn);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "0");
    EXPECT_EQ(pg.nodes[0]->covg, j);
    EXPECT_EQ(pg.nodes[0]->reads.size(), j);
    EXPECT_EQ(pg.nodes[0]->edges.size(), (uint)0);
    PanRead *pr;
    pr = new PanRead(1);
    EXPECT_EQ(pg.reads.size(), j);
    EXPECT_EQ(*pg.reads[1], *pr);
    EXPECT_EQ(pg.reads[1]->hits.size(), j);
    EXPECT_EQ(pg.reads[1]->hits[0].size(), (uint)0);
    EXPECT_EQ(pg.reads[1]->edges.size(), (uint)0);
    

    // add node again with same read
    pg.add_node(0,"0",1, mhs);

    EXPECT_EQ(pg.nodes.size(), j);
    EXPECT_EQ(*pg.nodes[0], *pn);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "0");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[0]->reads.size(), j);
    EXPECT_EQ(pg.nodes[0]->edges.size(), (uint)0);
    EXPECT_EQ(pg.reads.size(), j);
    EXPECT_EQ(*pg.reads[1], *pr);
    EXPECT_EQ(pg.reads[1]->hits.size(), j);
    EXPECT_EQ(pg.reads[1]->hits[0].size(), (uint)0);
    EXPECT_EQ(pg.reads[1]->edges.size(), (uint)0);

    // add node again with different read
    pg.add_node(0,"0",2, mhs);

    EXPECT_EQ(pg.nodes.size(), j);
    EXPECT_EQ(*pg.nodes[0], *pn);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "0");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)2);
    EXPECT_EQ(pg.nodes[0]->edges.size(), (uint)0);
    EXPECT_EQ(pg.reads.size(), (uint)2);
    EXPECT_EQ(*pg.reads[1], *pr);
    EXPECT_EQ(pg.reads[1]->hits.size(), j);
    EXPECT_EQ(pg.reads[1]->hits[0].size(), (uint)0);
    EXPECT_EQ(pg.reads[1]->edges.size(), (uint)0);
    EXPECT_EQ(pg.reads[2]->hits.size(), j);
    EXPECT_EQ(pg.reads[2]->hits[0].size(), (uint)0);
    EXPECT_EQ(pg.reads[2]->edges.size(), (uint)0);

    delete pn;
    delete pr;

    // add different node
    pg.add_node(1,"1",2, mhs);
    pn = new PanNode(1,1,"1");
    EXPECT_EQ(pg.nodes.size(), (uint)2);
    EXPECT_EQ(*pg.nodes[1], *pn);
    EXPECT_EQ(pg.nodes[1]->node_id, j);
    EXPECT_EQ(pg.nodes[1]->prg_id, j);
    EXPECT_EQ(pg.nodes[1]->name, "1");
    EXPECT_EQ(pg.nodes[1]->covg, j);
    EXPECT_EQ(pg.nodes[1]->reads.size(), j);
    EXPECT_EQ(pg.nodes[1]->edges.size(), (uint)0);
    EXPECT_EQ(pg.reads.size(), (uint)2);
    EXPECT_EQ(pg.reads[2]->hits.size(), (uint)2);
    EXPECT_EQ(pg.reads[2]->hits[1].size(), (uint)0);
    EXPECT_EQ(pg.reads[2]->edges.size(), (uint)0);

    delete pn;

    // add a node with hits
    Path p;
    deque<Interval> d = {Interval(0,1), Interval(4,7)};
    p.initialize(d);
    MinimizerHit* mh0;
    mh0 = new MinimizerHit(2, Interval(1,5), 2, p, 0, true);
    mhs.insert(mh0);
    d = {Interval(0,1), Interval(5,8)};
    p.initialize(d);
    MinimizerHit* mh1;
    mh1 = new MinimizerHit(2, Interval(1,5), 2, p, 0, true);
    mhs.insert(mh1);
    pg.add_node(2,"2",2, mhs);
    pn = new PanNode(2,2,"2");
    EXPECT_EQ(pg.nodes.size(), (uint)3);
    EXPECT_EQ(*pg.nodes[2], *pn);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->name, "2");
    EXPECT_EQ(pg.nodes[2]->covg, j);
    EXPECT_EQ(pg.nodes[2]->reads.size(), j);
    EXPECT_EQ(pg.nodes[2]->edges.size(), (uint)0);
    EXPECT_EQ(pg.reads.size(), (uint)2);
    EXPECT_EQ(pg.reads[2]->hits.size(), (uint)3);
    EXPECT_EQ(pg.reads[2]->hits[2].size(), (uint)2);
    EXPECT_EQ(pg.reads[2]->edges.size(), (uint)0);

    // expect death if some hit doesn't match the prg id expect
    MinimizerHit* mh2;
    mh2 = new MinimizerHit(0, Interval(1,5), 0, p, 0, true);
    mhs.insert(mh2);
    EXPECT_DEATH(pg.add_node(0,"0",0, mhs), "");
    delete pn;
    delete mh0;
    delete mh1;
    delete mh2;
}

TEST_F(PanGraphTest, add_node_sample)
{
    // add node and check it's there
    PanGraph pg;

    LocalPRG* l0;
    l0 = new LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA");
    vector<KmerNode*> kmp;
    
    pg.add_node(0, "zero", "sample", kmp, l0);

    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->edges.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)1);

    EXPECT_EQ(pg.samples.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->edges.size(), (uint)0);

    EXPECT_EQ(pg.edges.size(), (uint)0);
    EXPECT_EQ(pg.reads.size(), (uint)0);

    // add a second time
    pg.add_node(0, "zero", "sample", kmp, l0);
    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->edges.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)1);

    EXPECT_EQ(pg.samples.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->edges.size(), (uint)0);

    EXPECT_EQ(pg.edges.size(), (uint)0);
    EXPECT_EQ(pg.reads.size(), (uint)0);

    // add a node with a different sample
    pg.add_node(0, "zero", "sample1", kmp, l0);
    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->edges.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)2);
    
    EXPECT_EQ(pg.samples.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->edges.size(), (uint)0);
    EXPECT_EQ(pg.samples["sample1"]->name, "sample1");
    EXPECT_EQ(pg.samples["sample1"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample1"]->paths[0].size(), (uint)1);
    EXPECT_EQ(pg.samples["sample1"]->edges.size(), (uint)0);

    EXPECT_EQ(pg.edges.size(), (uint)0);
    EXPECT_EQ(pg.reads.size(), (uint)0);
    
    // add a node with a different prg
    pg.add_node(1, "one", "sample1", kmp, l0);
    EXPECT_EQ(pg.nodes.size(), (uint)2);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->edges.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->name, "one");
    EXPECT_EQ(pg.nodes[1]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[1]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[1]->edges.size(), (uint)0);
    EXPECT_EQ(pg.nodes[1]->samples.size(), (uint)1);
    
    EXPECT_EQ(pg.samples.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->edges.size(), (uint)0);
    EXPECT_EQ(pg.samples["sample1"]->name, "sample1");
    EXPECT_EQ(pg.samples["sample1"]->paths.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample1"]->paths[0].size(), (uint)1);
    EXPECT_EQ(pg.samples["sample1"]->paths[1].size(), (uint)1);
    EXPECT_EQ(pg.samples["sample1"]->edges.size(), (uint)0);

    EXPECT_EQ(pg.edges.size(), (uint)0);
    EXPECT_EQ(pg.reads.size(), (uint)0);

    delete l0;
}



TEST_F(PanGraphTest, add_edge)
{
    set<MinimizerHit*, pComp> mhs;
    PanGraph pg;
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",0, mhs);
    pg.add_edge(0,1,3,0); //++

    PanNode *pn1;
    pn1 = new PanNode(0,0,"0");
    PanNode *pn2;
    pn2 = new PanNode(1,1,"1");
    PanEdge *pe;
    pe = new PanEdge(pn1, pn2, 3);

    EXPECT_EQ(*pg.nodes[0], *pn1);
    EXPECT_EQ(*pg.nodes[1], *pn2);
    EXPECT_EQ(*pg.edges[0], *pe);

    // expect failure if a node doesn't exist in the graph
    EXPECT_DEATH(pg.add_edge(0,4,0,0),"");
    delete pn1;
    delete pn2;
    delete pe;
}

TEST_F(PanGraphTest, equals)
{
    set<MinimizerHit*, pComp> mhs;
    PanGraph pg1;
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_node(2,"2",2, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_edge(1,2,3,0);
  
    PanGraph pg2;
    pg2.add_node(1,"1",2, mhs);
    pg2.add_node(0,"0",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(2,"2",2, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(1,2,3,0);

    // adding nodes and edges in different order should make no difference
    EXPECT_EQ(pg1, pg1);
    EXPECT_EQ(pg2, pg2);
    EXPECT_EQ(pg1, pg2);
    EXPECT_EQ(pg2, pg1);

    // adding an extra edge does make a difference
    pg2.add_edge(0,2,3,0);
    EXPECT_EQ((pg1 == pg2), false);
    EXPECT_EQ((pg2 == pg1), false);

    // having one fewer edge makes a difference
    PanGraph pg3;
    pg3.add_node(1,"1",2, mhs);
    pg3.add_node(0,"0",0, mhs);
    pg3.add_node(2,"2",2, mhs);
    pg3.add_node(1,"1",0, mhs);
    pg3.add_edge(1,2,3,0);
    EXPECT_EQ((pg1 == pg3), false);
    EXPECT_EQ((pg3 == pg1), false);

    // or one extra node
    pg3.add_edge(0,1,3,0);
    EXPECT_EQ((pg1 == pg3), true); //adds the missing edge
    EXPECT_EQ((pg3 == pg1), true); //adds the missing edge
    pg3.add_node(3,"3",0, mhs);
    EXPECT_EQ((pg1 == pg3), false);
    EXPECT_EQ((pg3 == pg1), false);

    // should not break when have a cycle in pangraph
    pg3.add_edge(2,0,3,0);
    EXPECT_EQ(pg3, pg3);      

    // having edges orientated in a complementary way shouldn't make a difference (A->B == B- -> A-)
    pg1.add_edge(2,0,0,0);
    EXPECT_EQ(pg1, pg2);
    EXPECT_EQ(pg2, pg1);
}

TEST_F(PanGraphTest, not_equals)
{
    set<MinimizerHit*, pComp> mhs;
    PanGraph pg1;
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_node(2,"2",2, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_edge(1,2,3,0);

    PanGraph pg2;
    pg2.add_node(1,"1",2, mhs);
    pg2.add_node(0,"0",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(2,"2",2, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(1,2,3,0);

    // adding nodes and edges in different order should make no difference
    EXPECT_EQ((pg1!=pg1), false);
    EXPECT_EQ((pg2!=pg2), false);
    EXPECT_EQ((pg1!=pg2), false);
    EXPECT_EQ((pg2!=pg1), false);

    // adding an extra edge does make a difference
    pg2.add_edge(0,2,3,0);
    EXPECT_EQ((pg1 != pg2), true);
    EXPECT_EQ((pg2 != pg1), true);

    // having one fewer edge makes a difference
    PanGraph pg3;
    pg3.add_node(1,"1",2, mhs);
    pg3.add_node(0,"0",0, mhs);
    pg3.add_node(2,"2",2, mhs);
    pg3.add_node(1,"1",0, mhs);
    pg3.add_edge(1,2,3,0);
    EXPECT_EQ((pg1 != pg3), true);
    EXPECT_EQ((pg3 != pg1), true);

    // or one extra node
    pg3.add_edge(0,1,3,0);
    EXPECT_EQ((pg1 != pg3), false); //adds the missing edge
    EXPECT_EQ((pg3 != pg1), false); //adds the missing edge
    pg3.add_node(3,"3",0, mhs);
    EXPECT_EQ((pg1 != pg3), true);
    EXPECT_EQ((pg3 != pg1), true);

    // should not break when have a cycle in pangraph
    pg3.add_edge(2,0,3,0);
    EXPECT_EQ((pg3!=pg3), false);

    // having edges orientated in a complementary way shouldn't make a difference (A->B == B- -> A-)
    pg1.add_edge(2,0,0,0);
    EXPECT_EQ((pg1!=pg2), false);
    EXPECT_EQ((pg2!=pg1), false);
}

TEST_F(PanGraphTest, remove_edge)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_edge(1,2,3,0);
    pg1.add_node(3,"3",0, mhs);
    pg1.add_edge(2,3,3,0);

    // read 0: 0->1 2->3
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(2,"2",0, mhs);
    pg2.add_node(3,"3",0, mhs);
    pg2.add_edge(2,3,3,0);

    pg1.remove_edge(pg1.edges[1]);
    EXPECT_EQ(pg1, pg2);
}

TEST_F(PanGraphTest, remove_node)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_edge(1,2,3,0);
    pg1.add_node(3,"3",0, mhs);
    pg1.add_edge(2,3,3,0);

    // read 0: 0->1 3
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(3,"3",0, mhs);

    pg1.remove_node(pg1.nodes[2]);
    EXPECT_EQ(pg1, pg2);
}

TEST_F(PanGraphTest, add_shortcut_edge)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_edge(1,2,3,0);
    pg1.add_node(3,"3",0, mhs);
    pg1.add_edge(2,3,3,0);

    // read 0: 0->1->3
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(2,"2",0, mhs);
    pg2.add_edge(1,2,3,0);
    pg2.add_node(3,"3",0, mhs);
    pg2.add_edge(2,3,3,0);
    pg2.add_edge(1,3,3,0);
    pg2.edges[1]->covg = 0;
    pg2.edges[2]->covg = 0;

    vector<PanEdge*>::iterator e1 = pg1.edges.begin()+1, e2 = pg1.edges.begin()+2;
    pg1.add_shortcut_edge(e1, pg1.reads[0]);
    EXPECT_EQ(pg1, pg2);
}

TEST_F(PanGraphTest, split_node_by_edges)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_edge(1,2,3,0);
    pg1.add_node(3,"3",0, mhs);
    pg1.add_edge(2,3,3,0);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(1,"1",1, mhs);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);

    // read 0: 0->1->2->3 
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(2,"2",0, mhs);
    pg2.add_edge(1,2,3,0);
    pg2.add_node(3,"3",0, mhs);
    pg2.add_edge(2,3,3,0);
    // read 1: -4 -> -3 -> -5
    pg2.add_node(4,"4",1, mhs);
    pg2.add_edge(4,3,0,1);
    pg2.add_node(5,"5",1, mhs);
    pg2.add_edge(3,5,0,1);
    pg2.add_node(3,"3",1, mhs);
    // read 2: 0 -> 5 -> 3 -> 4
    pg2.add_edge(0,5,3,2);
    pg2.add_edge(5,3,3,2);
    pg2.add_edge(3,4,3,2);
    pg2.add_node(0,"0",2, mhs);
    pg2.add_node(5,"5",2, mhs);
    pg2.add_node(3,"3",2, mhs);
    pg2.add_node(4,"4",2, mhs); 

    // run split
    pg1.split_node_by_edges(pg1.nodes[1], pg1.edges[4], pg1.edges[0]);
    
    EXPECT_EQ(pg2, pg1);
    EXPECT_EQ(pg1, pg2);
}

TEST_F(PanGraphTest, split_nodes_by_reads)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_edge(1,2,3,0);
    pg1.add_node(3,"3",0, mhs);
    pg1.add_edge(2,3,3,0);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(1,"1",1, mhs);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);

    // read 0: 0->5->2->3 
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(5,"5",0, mhs);
    pg2.add_edge(0,5,3,0);
    pg2.add_node(2,"2",0, mhs);
    pg2.add_edge(5,2,3,0);
    pg2.add_node(3,"3",0, mhs);
    pg2.add_edge(2,3,3,0);
    // read 1: -4 -> -6 -> -1
    pg2.add_node(4,"4",1, mhs);
    pg2.add_node(6,"6",1, mhs);
    pg2.add_edge(4,6,0,1);
    pg2.add_node(1,"1",1, mhs);
    pg2.add_edge(6,1,0,1);
    // read 2: 0 -> 1 -> 6 -> 4
    pg2.add_edge(0,1,3,2);
    pg2.add_edge(1,6,3,2);
    pg2.add_edge(6,4,3,2);
    pg2.add_node(0,"0",2, mhs);
    pg2.add_node(1,"1",2, mhs);
    pg2.add_node(6,"6",2, mhs);
    pg2.add_node(4,"4",2, mhs);

    // run split
    pg1.split_nodes_by_reads(2,0);

    EXPECT_EQ(pg2, pg1);
    EXPECT_EQ(pg1, pg2);

    // second example, with a 3 way split
    PanGraph pg3, pg4;
    // read 0: 0->1->2->3
    pg3.add_node(0,"0",0, mhs);
    pg3.add_node(1,"1",0, mhs);
    pg3.add_edge(0,1,3,0);
    pg3.add_node(2,"2",0, mhs);
    pg3.add_edge(1,2,3,0);
    pg3.add_node(3,"3",0, mhs);
    pg3.add_edge(2,3,3,0);
    // read 1: -4 -> -3 -> -1
    pg3.add_node(4,"4",1, mhs);
    pg3.add_edge(4,3,0,1);
    pg3.add_edge(3,1,0,1);
    pg3.add_node(3,"3",1, mhs);
    pg3.add_node(1,"1",1, mhs);
    // read 2: 0 -> 1 -> 3 -> 4
    pg3.add_edge(0,1,3,2);
    pg3.add_edge(1,3,3,2);
    pg3.add_edge(3,4,3,2);
    pg3.add_node(0,"0",2, mhs);
    pg3.add_node(1,"1",2, mhs);
    pg3.add_node(3,"3",2, mhs);
    pg3.add_node(4,"4",2, mhs); 
    // read 3: 0 -> 1 -> 5 -> 4
    pg3.add_node(0,"0",4, mhs);
    pg3.add_node(1,"1",4, mhs);
    pg3.add_edge(0,1,3,4);
    pg3.add_node(5,"5",4, mhs);
    pg3.add_node(4,"4",4, mhs);
    pg3.add_edge(1,5,3,4);
    pg3.add_edge(5,4,3,4);

    cout << endl << "new test " << endl << endl;
    // read 0: 0->6->2->3
    pg4.add_node(0,"0",0, mhs);
    pg4.add_node(6,"6",0, mhs);
    pg4.add_edge(0,6,3,0);
    pg4.add_node(2,"2",0, mhs);
    pg4.add_edge(6,2,3,0);
    pg4.add_node(3,"3",0, mhs);
    pg4.add_edge(2,3,3,0);
    // read 1: -4 -> -8 -> -7
    pg4.add_node(4,"4",1, mhs);
    pg4.add_node(8,"8",1, mhs);
    pg4.add_edge(4,8,0,1);
    pg4.add_node(7,"7",1, mhs);
    pg4.add_edge(8,7,0,1);
    // read 2: 0 -> 7 -> 8 -> 4
    pg4.add_edge(0,7,3,2);
    pg4.add_edge(7,8,3,2);
    pg4.add_edge(8,4,3,2);
    pg4.add_node(0,"0",2, mhs);
    pg4.add_node(7,"7",2, mhs);
    pg4.add_node(8,"8",2, mhs);
    pg4.add_node(4,"4",2, mhs);
    // read 3: 0 -> 1 -> 5 -> 4
    pg4.add_node(0,"0",4, mhs);
    pg4.add_node(1,"1",4, mhs);
    pg4.add_node(5,"5",4, mhs);
    pg4.add_node(4,"4",4, mhs);
    pg4.add_edge(0,1,3,4);
    pg4.add_edge(1,5,3,4);
    pg4.add_edge(5,4,3,4);
    
    // run split
    pg3.split_nodes_by_reads(2,0);

    EXPECT_EQ(pg4, pg3);
    EXPECT_EQ(pg3, pg4);
}

TEST_F(PanGraphTest, read_clean)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_edge(1,2,3,0);
    pg1.add_node(3,"3",0, mhs);
    pg1.add_edge(2,3,3,0);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    // read 3: 0 -> -5
    pg1.add_node(5,"5",3, mhs);
    pg1.add_edge(0,5,1,3);
    // read 4: 5 -> -2
    pg1.add_edge(5,2,2,4);
    // read 5: 3 -> 5 -> 3
    pg1.add_edge(3,5,3,5);
    pg1.add_edge(5,3,3,5);
    // read 6: 0 -> 1 -> 4 -> 1
    pg1.add_edge(0,1,3,6);
    pg1.add_edge(1,4,3,6);
    pg1.add_edge(4,1,3,6);
    // read 7: 0 -> 2 -> 0 -> 4
    pg1.add_edge(0,2,3,8);
    pg1.add_edge(2,0,3,8);
    pg1.add_edge(0,4,3,8);
    // read 8: 4 -> -3 -> 4 -> 5 -> 1
    pg1.add_edge(4,3,1,7);
    pg1.add_edge(3,4,2,7);
    pg1.add_edge(4,5,3,7);
    pg1.add_edge(5,1,3,7);

    // read 0: 0->1->3
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(2,"2",0, mhs);
    pg2.add_node(3,"3",0, mhs);
    pg2.add_edge(1,3,3,0);
    pg2.add_edge(1,2,3,0); // with covg 0
    pg2.edges.back()->covg -= 1;
    pg2.add_edge(2,3,3,0); // with covg 0
    pg2.edges.back()->covg -= 1;
    // read 1: -4 -> -3 -> -1
    pg2.add_node(4,"4",1, mhs);
    pg2.add_edge(4,3,0,1);
    pg2.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg2.add_edge(0,1,3,2);
    pg2.add_edge(1,3,3,2);
    pg2.add_edge(3,4,3,2);
    // read 3: 0 -> -5
    pg2.add_node(5,"5",3, mhs);
    pg2.add_edge(0,5,1,3);
    // read 4: 5 -> -2
    pg2.add_edge(5,2,2,4);    
    // read 5: 0 -> 5 -> 0
    pg2.add_edge(3,5,3,5);
    pg2.edges.back()->covg -= 1;
    pg2.add_edge(5,3,3,5);
    pg2.edges.back()->covg -= 1;
    // read 6: 0 -> 1 -> 4 -> 1
    pg2.add_edge(0,1,3,6);
    pg2.add_edge(1,4,3,6);
    pg2.edges.back()->covg -= 1;
    pg2.add_edge(4,1,3,6);
    pg2.edges.back()->covg -= 1;
    // read 7: 0 -> 2 -> 0 -> 4 -> 5
    pg2.add_edge(0,2,3,8);
    pg2.edges.back()->covg -= 1;
    pg2.add_edge(2,0,3,8);
    pg2.edges.back()->covg -= 1;
    pg2.add_edge(0,4,3,8);
    pg2.edges.back()->covg -= 1;
    // read 8: 4 -> 3 -> 4 -> 5 -> 1
    pg2.add_edge(4,3,1,7);
    pg2.edges.back()->covg -= 1;
    pg2.add_edge(3,4,2,7);
    pg2.edges.back()->covg -= 1;
    pg2.add_edge(4,5,3,7);
    pg2.edges.back()->covg -= 1;
    pg2.add_edge(5,1,3,7);
    pg2.edges.back()->covg -= 1;

    pg1.read_clean(1);
    EXPECT_EQ(pg1, pg2);
}

TEST_F(PanGraphTest, remove_low_covg_nodes)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_edge(1,2,3,0);
    pg1.add_node(3,"3",0, mhs);
    pg1.add_edge(2,3,3,0);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(1,"1",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    // read 3: 0 -> -5
    pg1.add_node(0,"0",3, mhs);
    pg1.add_node(5,"5",3, mhs);
    pg1.add_edge(0,5,1,3);
    // read 4: 5 -> -1
    pg1.add_node(1,"1",4, mhs);
    pg1.add_node(5,"5",4, mhs);
    pg1.add_edge(5,1,2,4);

    // read 0: 0->1 3
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(3,"3",0, mhs);
    // read 1: -4 -> -3 -> -1
    pg2.add_node(4,"4",1, mhs);
    pg2.add_edge(4,3,0,1);
    pg2.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg2.add_edge(0,1,3,2);
    pg2.add_edge(1,3,3,2);
    pg2.add_edge(3,4,3,2);
    // read 3: 0 -> -5
    pg2.add_node(5,"5",3, mhs);
    pg2.add_edge(0,5,1,3);
    // read 4: 5 -> -1
    pg2.add_edge(5,1,2,4);

    pg1.remove_low_covg_nodes(1);
    EXPECT_EQ(pg1, pg2);
}

TEST_F(PanGraphTest, remove_low_covg_edges)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_edge(1,2,3,0);
    pg1.add_node(3,"3",0, mhs);
    pg1.add_edge(2,3,3,0);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    // read 3: 0 -> -5
    pg1.add_node(5,"5",3, mhs);
    pg1.add_edge(0,5,1,3);
    // read 4: 5 -> -1
    pg1.add_edge(5,1,2,4);

    // read 0: 0->1 2 3
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(2,"2",0, mhs);
    pg2.add_node(3,"3",0, mhs);
    // read 1: -4 -> -3 -> -1
    pg2.add_node(4,"4",1, mhs);
    pg2.add_edge(4,3,0,1);
    pg2.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg2.add_edge(0,1,3,2);
    pg2.add_edge(1,3,3,2);
    pg2.add_edge(3,4,3,2);
    // read 3: 0  -5
    pg2.add_node(5,"5",3, mhs);
    // read 4: 5  -1
    
    pg1.remove_low_covg_edges(1);
    EXPECT_EQ(pg1, pg2);
}

TEST_F(PanGraphTest, clean)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_edge(0,1,3,0);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_edge(1,2,3,0);
    pg1.add_node(3,"3",0, mhs);
    pg1.add_edge(2,3,3,0);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(1,"1",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    // read 3: 0 -> -5
    pg1.add_node(0,"0",3, mhs);
    pg1.add_node(5,"5",3, mhs);
    pg1.add_edge(0,5,1,3);
    // read 4: 5 -> -1
    pg1.add_node(1,"1",4, mhs);
    pg1.add_node(5,"5",4, mhs);
    pg1.add_edge(5,1,2,4);

    // repeats
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(1,"1",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(1,"1",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(1,"1",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(1,"1",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(1,"1",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(4,"4",1, mhs);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(1,"1",1, mhs);
    pg1.add_edge(4,3,0,1);
    pg1.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);
    pg1.add_edge(0,1,3,2);
    pg1.add_edge(1,3,3,2);
    pg1.add_edge(3,4,3,2);

    // read 0: 0->1->3
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(3,"3",0, mhs);
    pg2.add_edge(1,3,3,0);
    // read 1: -4 -> -3 -> -1
    pg2.add_node(4,"4",1, mhs);
    pg2.add_edge(4,3,0,1);
    pg2.add_edge(3,1,0,1);
    // read 2: 0 -> 1 -> 3 -> 4
    pg2.add_edge(0,1,3,2);
    pg2.add_edge(1,3,3,2);
    pg2.add_edge(3,4,3,2);

    pg1.clean(40);
    EXPECT_EQ(pg1, pg2);
}

TEST_F(PanGraphTest, write_gfa)
{
    set<MinimizerHit*, pComp> mhs;

    PanGraph pg2;
    pg2.add_node(1,"1",2, mhs);
    pg2.add_node(0,"0",0, mhs);
    pg2.add_edge(0,1,3,0);
    pg2.add_node(2,"2",2, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_edge(1,2,3,0);
    pg2.add_edge(1,2,3,0);
    pg2.write_gfa("../test/test_cases/pangraph_test_save.gfa");
}

TEST_F(PanGraphTest, save_matrix)
{
    // add node and check it's there
    PanGraph pg;

    LocalPRG* l0;
    l0 = new LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA");
    vector<KmerNode*> kmp;
   
    pg.add_node(0, "zero", "sample1", kmp, l0);
    pg.add_node(0, "zero", "sample1", kmp, l0);
    pg.add_node(0, "zero", "sample2", kmp, l0);
    pg.add_node(1, "one", "sample1", kmp, l0);
    pg.add_node(2, "two", "sample3", kmp, l0);
    
    pg.save_matrix("../test/test_cases/pangraph_test_save.matrix");
}
