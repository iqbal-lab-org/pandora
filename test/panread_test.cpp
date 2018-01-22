#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "pangenome/ns.cpp"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "pangenome_graph_class.h"
#include "minihit.h"
#include <stdint.h>
#include <iostream>
#include <set>

using namespace pangenome;

class PangenomeReadTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};

TEST_F(PangenomeReadTest,create){

    Read pr(3);
    EXPECT_EQ((uint)3, pr.id);
    EXPECT_EQ((uint)0, pr.nodes.size());
    EXPECT_EQ((uint)0, pr.node_orientations.size());
    EXPECT_EQ((uint)0, pr.hits.size());
}

TEST_F(PangenomeReadTest,add_hits)
{
    Read pr1(1);
    set<MinimizerHitPtr, pComp> c;

    pr1.add_hits(4, c);
    EXPECT_EQ((uint)1, pr1.hits.size());
    EXPECT_EQ((uint)0, pr1.hits[4].size());

    // add again and should append pr1.hits
    Interval i(0,5);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MinimizerHitPtr mh (make_shared<MinimizerHit>(4, i, 0, p, 0, 0));

    c.insert(mh);
    pr1.add_hits(4, c);
    EXPECT_EQ((uint)1, pr1.hits.size());
    EXPECT_EQ((uint)1, pr1.hits[4].size());

    // add to new id
    pr1.add_hits(5, c);
    EXPECT_EQ((uint)2, pr1.hits.size());
    EXPECT_EQ((uint)1, pr1.hits[5].size());
    EXPECT_EQ((uint)1, pr1.hits[4].size());

}

TEST_F(PangenomeReadTest, find_position)
{
    set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg;

    // read 0: 0->1->2->3->5->0->7->2->3->5->9
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",0, mhs);
    pg.add_node(2,"2",0, mhs);
    pg.add_node(3,"3",0, mhs);
    pg.add_node(5,"5",0, mhs);
    pg.add_node(0,"0",0, mhs);
    pg.add_node(7,"7",0, mhs);
    pg.add_node(2,"2",0, mhs);
    pg.add_node(3,"3",0, mhs);
    pg.add_node(5,"5",0, mhs);
    pg.add_node(9,"9",0, mhs);

    // read 1: 0->1->2
    pg.add_node(0,"0",1, mhs);
    pg.add_node(1,"1",1, mhs);
    pg.add_node(2,"2",1, mhs);

    pg.reads[0]->node_orientations[6] = 1;

    vector<uint16_t> v = {2,3,5};
    vector<bool> b = {0,0,0};
    uint p = pg.reads[0]->find_position(v,b);
    EXPECT_EQ(p,(uint)2);

    // one at the end of the string
    v = {3,5,9};
    p = pg.reads[0]->find_position(v,b);
    EXPECT_EQ(p,(uint)8);

    // one in reverse
    v = {0,5,3};
    b = {1,1,1};
    p = pg.reads[0]->find_position(v,b);
    EXPECT_EQ(p,(uint)3);

    // one overlapping start
    v = {9,0,1};
    b = {0,0,0};
    p = pg.reads[0]->find_position(v,b);
    EXPECT_EQ(p,(uint)0);

    // one in reverse overlapping start
    v = {1,0,9};
    b = {1,1,1};
    p = pg.reads[0]->find_position(v,b);
    EXPECT_EQ(p,(uint)0);

    // one overlapping the end
    b = {0,0,0};
    v = {5,9,9};
    p = pg.reads[0]->find_position(v,b);
    EXPECT_EQ(p,(uint)9);

    // one in reverse overlapping end
    b = {1,1,1};
    v = {0,9,5};
    p = pg.reads[0]->find_position(v,b);
    EXPECT_EQ(p,(uint)9);

    // one not a match
    b = {0,0,0};
    v = {8,8,8};
    p = pg.reads[0]->find_position(v,b);
    EXPECT_EQ(p,std::numeric_limits<uint>::max());

    // one where orientations mean not a match
    v = {3,2,7};
    p = pg.reads[0]->find_position(v,b);
    EXPECT_EQ(p,std::numeric_limits<uint>::max());

    // and when is whole read
    v = {0,1,2};
    p = pg.reads[1]->find_position(v,b);
    EXPECT_EQ(p,(uint)0);
}

TEST_F(PangenomeReadTest,remove_node)
{

   set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg;
    vector<NodePtr> exp_read_nodes;
    vector<bool> exp_read_orientations;

    // read 0: 0->1->2->3
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",0, mhs);
    pg.add_node(2,"2",0, mhs);
    pg.add_node(3,"3",0, mhs);

    // read 1: -4 -> -3 -> -1
    pg.add_node(4,"4",1, mhs);
    pg.add_node(3,"3",1, mhs);
    pg.add_node(1,"1",1, mhs);

    // read 2: 0 -> 1 -> 3 -> 4
    pg.add_node(0,"0",2, mhs);
    pg.add_node(1,"1",2, mhs);
    pg.add_node(3,"3",2, mhs);
    pg.add_node(4,"4",2, mhs);

    // check all as expected
    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[2], pg.nodes[3]};
    exp_read_orientations = {0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[0]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[1]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4]};
    exp_read_orientations = {0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[2]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example with a node replacing an old node which only appears in one read
    pg.reads[0]->remove_node(pg.nodes[2]);

    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[3]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[0]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[1]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4]};
    exp_read_orientations = {0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[2]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example where old node appears in more than one read
    pg.reads[0]->remove_node(pg.nodes[1]);

    exp_read_nodes = {pg.nodes[0], pg.nodes[3]};
    exp_read_orientations = {0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[0]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[1]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4]};
    exp_read_orientations = {0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[2]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example where have actual hit
    Interval i(0,5);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MinimizerHitPtr mh (make_shared<MinimizerHit>(4, i, 0, p, 0, 0));
    set<MinimizerHitPtr, pComp> c;
    c.insert(mh);
    pg.reads[2]->add_hits(4, c);

    pg.reads[2]->remove_node(pg.nodes[4]);

    exp_read_nodes = {pg.nodes[0], pg.nodes[3]};
    exp_read_orientations = {0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[0]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[1]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[3]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[2]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    //example where node appears twice in read
    pg.add_node(1,"1",2, mhs);
    pg.reads[2]->remove_node(pg.nodes[1]);

    exp_read_nodes = {pg.nodes[0], pg.nodes[3]};
    exp_read_orientations = {0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[0]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[1]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[0], pg.nodes[3]};
    exp_read_orientations = {0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[2]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);
}

TEST_F(PangenomeReadTest,remove_node_it)
{
    set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg;
    vector<NodePtr> exp_read_nodes;
    vector<bool> exp_read_orientations;

    // read 0: 0->1->2->3
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",0, mhs);
    pg.add_node(2,"2",0, mhs);
    pg.add_node(3,"3",0, mhs);

    // read 1: -4 -> -3 -> -1
    pg.add_node(4,"4",1, mhs);
    pg.add_node(3,"3",1, mhs);
    pg.add_node(1,"1",1, mhs);

    // read 2: 0 -> 1 -> 3 -> 4
    pg.add_node(0,"0",2, mhs);
    pg.add_node(1,"1",2, mhs);
    pg.add_node(3,"3",2, mhs);
    pg.add_node(4,"4",2, mhs);

    // check all as expected
    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[2], pg.nodes[3]};
    exp_read_orientations = {0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[0]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[1]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4]};
    exp_read_orientations = {0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[2]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example removing a node which only appears in one read
    auto it = find(pg.reads[0]->nodes.begin(), pg.reads[0]->nodes.end(),pg.nodes[2]);
    pg.reads[0]->remove_node(it);

    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[3]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[0]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[1]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4]};
    exp_read_orientations = {0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[2]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example where old node appears in more than one read
    pg.reads[0]->remove_node(pg.nodes[1]);

    exp_read_nodes = {pg.nodes[0], pg.nodes[3]};
    exp_read_orientations = {0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[0]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[1]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4]};
    exp_read_orientations = {0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[2]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example where have actual hit
    Interval i(0,5);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MinimizerHitPtr mh (make_shared<MinimizerHit>(4, i, 0, p, 0, 0));
    set<MinimizerHitPtr, pComp> c;
    c.insert(mh);
    pg.reads[2]->add_hits(4, c);

    pg.reads[2]->remove_node(pg.nodes[4]);

    exp_read_nodes = {pg.nodes[0], pg.nodes[3]};
    exp_read_orientations = {0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[0]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[1]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[0], pg.nodes[1], pg.nodes[3]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[2]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    //example where node appears twice in read
    pg.add_node(1,"1",2, mhs);
    pg.reads[2]->remove_node(pg.nodes[1]);

    exp_read_nodes = {pg.nodes[0], pg.nodes[3]};
    exp_read_orientations = {0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[0]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    exp_read_orientations = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[1]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = {pg.nodes[0], pg.nodes[3]};
    exp_read_orientations = {0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, pg.reads[2]->nodes, exp_read_nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);
}

TEST_F(PangenomeReadTest,replace_node)
{
    set<MinimizerHitPtr, pComp> mhs;

    Graph pg;
    // read 0: 0->1->2->3->1
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",0, mhs);
    pg.add_node(2,"2",0, mhs);
    pg.add_node(3,"3",0, mhs);
    pg.add_node(1,"1",0, mhs);

    // read 1: 4 -> 3 -> 1
    pg.add_node(4,"4",1, mhs);
    pg.add_node(3,"3",1, mhs);
    pg.add_node(1,"1",1, mhs);

    // check what we expect to start with
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[4]->node_id, (uint)4);
    EXPECT_EQ(pg.nodes[4]->covg, (uint)1);

    EXPECT_EQ((uint)2, pg.reads.size());
    vector<NodePtr> read_exp = {pg.nodes[0], pg.nodes[1], pg.nodes[2], pg.nodes[3], pg.nodes[1]};
    vector<bool> read_o_exp = {0,0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, read_exp, pg.reads[0]->nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, read_o_exp, pg.reads[0]->node_orientations);
    read_exp = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    read_o_exp = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, read_exp, pg.reads[1]->nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, read_o_exp, pg.reads[1]->node_orientations);

    // example with a node replacing an old node which only appears in one read
    NodePtr n = make_shared<Node>(2, 5, "2_prime");
    pg.nodes[5] = n;
    auto it = pg.reads[0]->nodes.begin()+2;
    pg.reads[0]->replace_node(it,n);

    EXPECT_EQ((uint)6, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[2]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);

    EXPECT_EQ((uint)2, pg.reads.size());
    read_exp = {pg.nodes[0], pg.nodes[1], pg.nodes[5], pg.nodes[3], pg.nodes[1]};
    read_o_exp = {0,0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, read_exp, pg.reads[0]->nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, read_o_exp, pg.reads[0]->node_orientations);
    read_exp = {pg.nodes[4], pg.nodes[3], pg.nodes[1]};
    read_o_exp = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, read_exp, pg.reads[1]->nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, read_o_exp, pg.reads[1]->node_orientations);

    // example where old node appears in more than one read
    n = make_shared<Node>(3, 6, "3_prime");
    pg.nodes[6] = n;
    it = pg.reads[1]->nodes.begin()+1;
    pg.reads[1]->replace_node(it,n);

    EXPECT_EQ((uint)7, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[2]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[6]->prg_id, (uint)3);

    EXPECT_EQ((uint)2, pg.reads.size());
    read_exp = {pg.nodes[0], pg.nodes[1], pg.nodes[5], pg.nodes[3], pg.nodes[1]};
    read_o_exp = {0,0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, read_exp, pg.reads[0]->nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, read_o_exp, pg.reads[0]->node_orientations);
    read_exp = {pg.nodes[4], pg.nodes[6], pg.nodes[1]};
    read_o_exp = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, read_exp, pg.reads[1]->nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, read_o_exp, pg.reads[1]->node_orientations);

    // example where move actual hits
    Interval i(0,5);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MinimizerHitPtr mh(make_shared<MinimizerHit>(4, i, 0, p, 0, 0));
    set<MinimizerHitPtr, pComp> c;
    c.insert(mh);
    pg.reads[1]->add_hits(4, c);
    EXPECT_EQ(pg.reads[1]->hits[4].size(),(uint)1);

    n = make_shared<Node>(4, 7, "4_prime");
    pg.nodes[7] = n;
    it = pg.reads[1]->nodes.begin();
    pg.reads[1]->replace_node(it,n);

    EXPECT_EQ((uint)8, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[2]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[6]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[7]->prg_id, (uint)4);

    EXPECT_EQ((uint)2, pg.reads.size());
    read_exp = {pg.nodes[0], pg.nodes[1], pg.nodes[5], pg.nodes[3], pg.nodes[1]};
    read_o_exp = {0,0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, read_exp, pg.reads[0]->nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, read_o_exp, pg.reads[0]->node_orientations);
    read_exp = {pg.nodes[7], pg.nodes[6], pg.nodes[1]};
    read_o_exp = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, read_exp, pg.reads[1]->nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, read_o_exp, pg.reads[1]->node_orientations);
    EXPECT_EQ(pg.reads[1]->hits[7].size(),(uint)1);

    //example where node appears twice in read
    n = make_shared<Node>(1, 8, "1_prime");
    pg.nodes[8] = n;
    it = pg.reads[0]->nodes.begin()+4;
    pg.reads[0]->replace_node(it,n);

    EXPECT_EQ((uint)9, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[2]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[6]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[7]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[8]->prg_id, (uint)1);


    EXPECT_EQ((uint)2, pg.reads.size());
    read_exp = {pg.nodes[0], pg.nodes[1], pg.nodes[5], pg.nodes[3], pg.nodes[8]};
    read_o_exp = {0,0,0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, read_exp, pg.reads[0]->nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, read_o_exp, pg.reads[0]->node_orientations);
    read_exp = {pg.nodes[7], pg.nodes[6], pg.nodes[1]};
    read_o_exp = {0,0,0};
    EXPECT_ITERABLE_EQ(vector<NodePtr>, read_exp, pg.reads[1]->nodes);
    EXPECT_ITERABLE_EQ(vector<bool>, read_o_exp, pg.reads[1]->node_orientations);
    EXPECT_EQ(pg.reads[1]->hits[7].size(),(uint)1);

}

TEST_F(PangenomeReadTest,equals){
    Read pr1(1);
    Read pr2(2);
    EXPECT_EQ(pr1, pr1);
    EXPECT_EQ(pr2, pr2);
    EXPECT_EQ((pr1==pr2), false);
    EXPECT_EQ((pr2==pr1), false);
}

TEST_F(PangenomeReadTest,nequals){
    Read pr1(1);
    Read pr2(2);
    EXPECT_NE(pr1, pr2);
    EXPECT_NE(pr2, pr1);
    EXPECT_EQ((pr1!=pr1), false);
    EXPECT_EQ((pr2!=pr2), false);
}

TEST_F(PangenomeReadTest,less){
    Read pr1(1);
    Read pr2(2);
    EXPECT_EQ((pr1<pr1), false);
    EXPECT_EQ((pr2<pr2), false);
    EXPECT_EQ((pr1<pr2), true);
    EXPECT_EQ((pr2<pr1), false);
}
