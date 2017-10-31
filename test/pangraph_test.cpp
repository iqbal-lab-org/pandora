#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "pangenome/ns.cpp"
#include "pangenome_graph_class.h"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "pangenome/pansample.h"
#include "minihit.h"
#include "localPRG.h"
#include <stdint.h>
#include <numeric>
#include <cassert>
#include <iostream>

using namespace pangenome;

class PangenomeGraphTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(PangenomeGraphTest, add_node)
{
    set<MinimizerHitPtr, pComp> mhs;

    // add node and check it's there
    PGraphTester pg;
    pg.add_node(0,"0",1, mhs);

    NodePtr pn = make_shared<Node>(0,0,"0");
    uint32_t j = 1;
    EXPECT_EQ(pg.nodes.size(), j);
    EXPECT_EQ(*pg.nodes[0], *pn);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "0");
    EXPECT_EQ(pg.nodes[0]->covg, j);
    EXPECT_EQ(pg.nodes[0]->reads.size(), j);
    ReadPtr pr = make_shared<Read>(1);
    EXPECT_EQ(pg.reads.size(), j);
    EXPECT_EQ(*pg.reads[1], *pr);
    EXPECT_EQ(pg.reads[1]->hits.size(), j);
    EXPECT_EQ(pg.reads[1]->hits[0].size(), (uint)0);


    // add node again with same read
    pg.add_node(0,"0",1, mhs);

    EXPECT_EQ(pg.nodes.size(), j);
    EXPECT_EQ(*pg.nodes[0], *pn);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "0");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)2);
    EXPECT_EQ(pg.reads.size(), j);
    EXPECT_EQ(*pg.reads[1], *pr);
    EXPECT_EQ(pg.reads[1]->hits.size(), j);
    EXPECT_EQ(pg.reads[1]->hits[0].size(), (uint)0);

    // add node again with different read
    pg.add_node(0,"0",2, mhs);

    EXPECT_EQ(pg.nodes.size(), j);
    EXPECT_EQ(*pg.nodes[0], *pn);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "0");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)3);
    EXPECT_EQ(pg.reads.size(), (uint)2);
    EXPECT_EQ(*pg.reads[1], *pr);
    EXPECT_EQ(pg.reads[1]->hits.size(), j);
    EXPECT_EQ(pg.reads[1]->hits[0].size(), (uint)0);
    EXPECT_EQ(pg.reads[2]->hits.size(), j);
    EXPECT_EQ(pg.reads[2]->hits[0].size(), (uint)0);

    // add different node
    pg.add_node(1,"1",2, mhs);
    pn = make_shared<Node>(1,1,"1");
    EXPECT_EQ(pg.nodes.size(), (uint)2);
    EXPECT_EQ(*pg.nodes[1], *pn);
    EXPECT_EQ(pg.nodes[1]->node_id, j);
    EXPECT_EQ(pg.nodes[1]->prg_id, j);
    EXPECT_EQ(pg.nodes[1]->name, "1");
    EXPECT_EQ(pg.nodes[1]->covg, j);
    EXPECT_EQ(pg.nodes[1]->reads.size(), j);
    EXPECT_EQ(pg.reads.size(), (uint)2);
    EXPECT_EQ(pg.reads[2]->hits.size(), (uint)2);
    EXPECT_EQ(pg.reads[2]->hits[1].size(), (uint)0);

    // add a node with hits
    Path p;
    deque<Interval> d = {Interval(0,1), Interval(4,7)};
    p.initialize(d);
    MinimizerHitPtr mh0 (make_shared<MinimizerHit>(2, Interval(1,5), 2, p, 0, true));
    mhs.insert(mh0);
    d = {Interval(0,1), Interval(5,8)};
    p.initialize(d);
    MinimizerHitPtr mh1 (make_shared<MinimizerHit>(2, Interval(1,5), 2, p, 0, true));
    mhs.insert(mh1);
    pg.add_node(2,"2",2, mhs);
    pn = make_shared<Node>(2,2,"2");
    EXPECT_EQ(pg.nodes.size(), (uint)3);
    EXPECT_EQ(*pg.nodes[2], *pn);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->name, "2");
    EXPECT_EQ(pg.nodes[2]->covg, j);
    EXPECT_EQ(pg.nodes[2]->reads.size(), j);
    EXPECT_EQ(pg.reads.size(), (uint)2);
    EXPECT_EQ(pg.reads[2]->hits.size(), (uint)3);
    EXPECT_EQ(pg.reads[2]->hits[2].size(), (uint)2);

    // expect death if some hit doesn't match the prg id expect
    MinimizerHitPtr mh2 (make_shared<MinimizerHit>(0, Interval(1,5), 0, p, 0, true));
    mhs.insert(mh2);
    EXPECT_DEATH(pg.add_node(0,"0",0, mhs), "");
}

TEST_F(PangenomeGraphTest, add_node_sample)
{
    // add node and check it's there
    PGraphTester pg;

    LocalPRG* l0;
    l0 = new LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA");
    vector<KmerNodePtr> kmp;
    
    pg.add_node(0, "zero", "sample", kmp, l0);

    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)1);

    EXPECT_EQ(pg.samples.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)1);

    EXPECT_EQ(pg.reads.size(), (uint)0);

    // add a second time
    pg.add_node(0, "zero", "sample", kmp, l0);
    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)1);

    EXPECT_EQ(pg.samples.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)2);

    EXPECT_EQ(pg.reads.size(), (uint)0);

    // add a node with a different sample
    pg.add_node(0, "zero", "sample1", kmp, l0);
    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)2);
    
    EXPECT_EQ(pg.samples.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)2);
    EXPECT_EQ(pg.samples["sample1"]->name, "sample1");
    EXPECT_EQ(pg.samples["sample1"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample1"]->paths[0].size(), (uint)1);

    EXPECT_EQ(pg.reads.size(), (uint)0);
    
    // add a node with a different prg
    pg.add_node(1, "one", "sample1", kmp, l0);
    EXPECT_EQ(pg.nodes.size(), (uint)2);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->name, "one");
    EXPECT_EQ(pg.nodes[1]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[1]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[1]->samples.size(), (uint)1);
    
    EXPECT_EQ(pg.samples.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)2);
    EXPECT_EQ(pg.samples["sample1"]->name, "sample1");
    EXPECT_EQ(pg.samples["sample1"]->paths.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample1"]->paths[0].size(), (uint)1);
    EXPECT_EQ(pg.samples["sample1"]->paths[1].size(), (uint)1);

    EXPECT_EQ(pg.reads.size(), (uint)0);

    delete l0;
}

TEST_F(PangenomeGraphTest, clear)
{
    // read pg
    set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg;
    pg.add_node(0,"0",1, mhs);
    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.reads.size(), (uint)1);
    EXPECT_EQ(pg.samples.size(), (uint)0);
    pg.clear();
    EXPECT_EQ(pg.nodes.size(), (uint)0);
    EXPECT_EQ(pg.reads.size(), (uint)0);
    EXPECT_EQ(pg.samples.size(), (uint)0);

    // sample pg
    LocalPRG* l0;
    l0 = new LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA");
    vector<KmerNodePtr> kmp;
    pg.add_node(0, "zero", "sample", kmp, l0);
    EXPECT_EQ(pg.reads.size(), (uint)0);
    EXPECT_EQ(pg.samples.size(), (uint)1);
    pg.clear();
    EXPECT_EQ(pg.nodes.size(), (uint)0);
    EXPECT_EQ(pg.reads.size(), (uint)0);
    EXPECT_EQ(pg.samples.size(), (uint)0);
    delete l0;
}


TEST_F(PangenomeGraphTest, equals)
{
    set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg1;
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_node(2,"2",2, mhs);

    PGraphTester pg2;
    pg2.add_node(1,"1",2, mhs);
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(2,"2",2, mhs);
    pg2.add_node(1,"1",0, mhs);

    // adding nodes in different order should make no difference
    EXPECT_EQ(pg1, pg1);
    EXPECT_EQ(pg2, pg2);
    EXPECT_EQ(pg1, pg2);
    EXPECT_EQ(pg2, pg1);

    // or one extra node
    pg2.add_node(3,"3",0, mhs);
    EXPECT_EQ((pg1 == pg2), false);
    EXPECT_EQ((pg2 == pg1), false);

    // should not break when have a cycle in pangraph
    pg1.add_node(0,"0",0, mhs);
    EXPECT_EQ(pg1, pg1);
}

TEST_F(PangenomeGraphTest, not_equals)
{
    set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg1;
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_node(2,"2",2, mhs);

    PGraphTester pg2;
    pg2.add_node(1,"1",2, mhs);
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(2,"2",2, mhs);
    pg2.add_node(1,"1",0, mhs);

    // adding nodes in different order should make no difference
    EXPECT_EQ((pg1!=pg1), false);
    EXPECT_EQ((pg2!=pg2), false);
    EXPECT_EQ((pg1!=pg2), false);
    EXPECT_EQ((pg2!=pg1), false);

    // or one extra node
    pg2.add_node(3,"3",0, mhs);
    EXPECT_EQ((pg1 != pg2), true);
    EXPECT_EQ((pg2 != pg1), true);

    // should not break when have a cycle in pangraph
    pg1.add_node(0,"0",0, mhs);
    EXPECT_EQ((pg1!=pg1), false);
}

TEST_F(PangenomeGraphTest, remove_node)
{
    set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_node(3,"3",0, mhs);

    // read 0: 0->1->3
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_node(3,"3",0, mhs);

    pg1.remove_node(pg1.nodes[2]);
    EXPECT_EQ(pg1, pg2);
}


TEST_F(PangenomeGraphTest, remove_low_covg_nodes)
{
    set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg1, pg2, pg3;
    // read 0: 0->1->2->3
    pg1.add_node(0,"0",0, mhs);
    pg1.add_node(1,"1",0, mhs);
    pg1.add_node(2,"2",0, mhs);
    pg1.add_node(3,"3",0, mhs);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(1,"1",1, mhs);
    pg1.add_node(3,"3",1, mhs);
    pg1.add_node(4,"4",1, mhs);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(0,"0",2, mhs);
    pg1.add_node(1,"1",2, mhs);
    pg1.add_node(3,"3",2, mhs);
    pg1.add_node(4,"4",2, mhs);
    // read 3: 0 -> 5
    pg1.add_node(0,"0",3, mhs);
    pg1.add_node(5,"5",3, mhs);
    // read 4: 5 -> 1
    pg1.add_node(5,"5",4, mhs);
    pg1.add_node(1,"1",4, mhs);

    // read 0: 0->1->3
    pg2.add_node(0,"0",0, mhs);
    pg2.add_node(1,"1",0, mhs);
    pg2.add_node(3,"3",0, mhs);
    // read 1: -4 -> -3 -> -1
    pg2.add_node(1,"1",1, mhs);
    pg2.add_node(3,"3",1, mhs);
    pg2.add_node(4,"4",1, mhs);
    // read 2: 0 -> 1 -> 3 -> 4
    pg2.add_node(0,"0",2, mhs);
    pg2.add_node(1,"1",2, mhs);
    pg2.add_node(3,"3",2, mhs);
    pg2.add_node(4,"4",2, mhs);
    // read 3: 0 -> 5
    pg2.add_node(0,"0",3, mhs);
    pg2.add_node(5,"5",3, mhs);
    // read 4: 5 -> 1
    pg2.add_node(5,"5",4, mhs);
    pg2.add_node(1,"1",4, mhs);

    pg1.remove_low_covg_nodes(1);
    EXPECT_EQ(pg1, pg2);

    // read 0: 0->1->3
    pg3.add_node(0,"0",0, mhs);
    pg3.add_node(1,"1",0, mhs);
    pg3.add_node(3,"3",0, mhs);
    // read 1: -4 -> -3 -> -1
    pg3.add_node(1,"1",1, mhs);
    pg3.add_node(3,"3",1, mhs);
    // read 2: 0 -> 1 -> 3 -> 4
    pg3.add_node(0,"0",2, mhs);
    pg3.add_node(1,"1",2, mhs);
    pg3.add_node(3,"3",2, mhs);
    // read 3: 0 -> 5
    pg3.add_node(0,"0",3, mhs);
    // read 4: 5 -> 1
    pg3.add_node(1,"1",4, mhs);

    pg1.remove_low_covg_nodes(2);
    EXPECT_EQ(pg1, pg3);
}

TEST_F(PangenomeGraphTest, add_hits_to_kmergraph)
{
}

TEST_F(PangenomeGraphTest, save_matrix)
{
    // add node and check it's there
    PGraphTester pg;

    LocalPRG* l0;
    l0 = new LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA");
    vector<KmerNodePtr> kmp;
   
    pg.add_node(0, "zero", "sample1", kmp, l0);
    pg.add_node(0, "zero", "sample1", kmp, l0);
    pg.add_node(0, "zero", "sample2", kmp, l0);
    pg.add_node(1, "one", "sample1", kmp, l0);
    pg.add_node(2, "two", "sample3", kmp, l0);
    
    pg.save_matrix("../test/test_cases/pangraph_test_save.matrix");
}
