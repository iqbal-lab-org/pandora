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

/*TEST_F(PangenomeReadTest,replace_node)
{
   set<MinimizerHitPtr, pComp> mhs;

    Graph pg;
    // read 0: 0->1->2->3
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",0, mhs);
    pg.add_edge(0,1,3,0);
    pg.add_node(2,"2",0, mhs);
    pg.add_edge(1,2,3,0);
    pg.add_node(3,"3",0, mhs);
    pg.add_edge(2,3,3,0);
    // read 1: -4 -> -3 -> -1
    pg.add_node(4,"4",1, mhs);
    pg.add_edge(4,3,0,1);
    pg.add_edge(3,1,0,1);
    pg.add_node(3,"3",1, mhs);
    pg.add_node(1,"1",1, mhs);
    // read 2: 0 -> 1 -> 3 -> 4
    pg.add_edge(0,1,3,2);
    pg.add_edge(1,3,3,2);
    pg.add_edge(3,4,3,2);
    pg.add_node(0,"0",2, mhs);
    pg.add_node(1,"1",2, mhs);
    pg.add_node(3,"3",2, mhs);
    pg.add_node(4,"4",2, mhs);

    // check all nodes and edges where expect to start with
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[4]->node_id, (uint)4);
    EXPECT_EQ(pg.nodes[4]->covg, (uint)2);
    EXPECT_EQ((uint)5, pg.edges.size());
    EXPECT_EQ(pg.edges[0]->from->node_id, (uint)0);
    EXPECT_EQ(pg.edges[0]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[0]->covg, (uint)2);
    EXPECT_EQ(pg.edges[1]->from->node_id, (uint)1);
    EXPECT_EQ(pg.edges[1]->to->node_id, (uint)2);
    EXPECT_EQ(pg.edges[1]->covg, (uint)1);
    EXPECT_EQ(pg.edges[2]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[2]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[2]->covg, (uint)1);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)2);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)3, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);
    EXPECT_EQ(pg.edges[2], pg.reads[0]->edges[2]);
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]);
    EXPECT_EQ(pg.edges[3], pg.reads[2]->edges[2]);

    // example with a node replacing an old node which only appears in one read
    Node *n5;
    n5 = new Node(2, 5, "2");
    n5->covg = 0;
    pg.nodes[5] = n5;
 
    pg.reads[0]->replace_node(pg.nodes[2], pg.nodes[5]);
    
    EXPECT_EQ((uint)6, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)0);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[4]->node_id, (uint)4);
    EXPECT_EQ(pg.nodes[4]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[5]->node_id, (uint)5);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[5]->covg, (uint)1);
    EXPECT_EQ((uint)5, pg.edges.size());
    EXPECT_EQ(pg.edges[0]->from->node_id, (uint)0);
    EXPECT_EQ(pg.edges[0]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[0]->covg, (uint)2);
    EXPECT_EQ(pg.edges[1]->from->node_id, (uint)1);
    EXPECT_EQ(pg.edges[1]->to->node_id, (uint)2);
    EXPECT_EQ(pg.edges[1]->covg, (uint)1);
    EXPECT_EQ(pg.edges[2]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[2]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[2]->covg, (uint)1);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)2);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)3, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);
    EXPECT_EQ(pg.edges[2], pg.reads[0]->edges[2]);
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]);
    EXPECT_EQ(pg.edges[3], pg.reads[2]->edges[2]);
    EXPECT_EQ((uint)5, pg.reads[0]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[0].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[2].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[3].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[5].size());
    EXPECT_EQ((uint)3, pg.reads[1]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[3].size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[4].size());
    EXPECT_EQ((uint)4, pg.reads[2]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[0].size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[3].size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[4].size());
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)2);
    EXPECT_EQ(pg.nodes[1]->reads.size(), (uint)3);
    EXPECT_EQ(pg.nodes[2]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[3]->reads.size(), (uint)3);
    EXPECT_EQ(pg.nodes[4]->reads.size(), (uint)2);
    EXPECT_EQ(pg.nodes[5]->reads.size(), (uint)1);

    // example where old node appears in more than one read
    Node *n6;
    n6 = new Node(1, 6, "1");
    n6->covg = 0;
    pg.nodes[6] = n6;
 
    pg.reads[0]->replace_node(pg.nodes[1], pg.nodes[6]);
    
    EXPECT_EQ((uint)7, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)0);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[4]->node_id, (uint)4);
    EXPECT_EQ(pg.nodes[4]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[5]->node_id, (uint)5);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[5]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[6]->node_id, (uint)6);
    EXPECT_EQ(pg.nodes[6]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[6]->covg, (uint)1);
    EXPECT_EQ((uint)5, pg.edges.size());
    EXPECT_EQ(pg.edges[0]->from->node_id, (uint)0);
    EXPECT_EQ(pg.edges[0]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[0]->covg, (uint)2);
    EXPECT_EQ(pg.edges[1]->from->node_id, (uint)1);
    EXPECT_EQ(pg.edges[1]->to->node_id, (uint)2);
    EXPECT_EQ(pg.edges[1]->covg, (uint)1);
    EXPECT_EQ(pg.edges[2]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[2]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[2]->covg, (uint)1);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)2);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)3, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);
    EXPECT_EQ(pg.edges[2], pg.reads[0]->edges[2]);
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]);
    EXPECT_EQ(pg.edges[3], pg.reads[2]->edges[2]);
    EXPECT_EQ((uint)6, pg.reads[0]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[0].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[2].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[3].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[5].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[6].size());
    EXPECT_EQ((uint)3, pg.reads[1]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[3].size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[4].size());
    EXPECT_EQ((uint)4, pg.reads[2]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[0].size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[3].size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[4].size());
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)2);
    EXPECT_EQ(pg.nodes[1]->reads.size(), (uint)2);
    EXPECT_EQ(pg.nodes[2]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[3]->reads.size(), (uint)3);
    EXPECT_EQ(pg.nodes[4]->reads.size(), (uint)2);
    EXPECT_EQ(pg.nodes[5]->reads.size(), (uint)1);
    EXPECT_EQ(pg.nodes[6]->reads.size(), (uint)1);

    // example where move actual hits
    Interval i(0,5);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MinimizerHitPtr mh(make_shared<MinimizerHit>(4, i, 0, p, 0, 0));
    set<MinimizerHitPtr, pComp> c;
    c.insert(mh);
    pg.reads[2]->add_hits(4, c);

    Node *n7;
    n7 = new Node(4, 7, "4");
    n7->covg = 0;
    pg.nodes[7] = n7;
 
    pg.reads[2]->replace_node(pg.nodes[4], pg.nodes[7]);
    
    EXPECT_EQ((uint)8, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)0);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[4]->node_id, (uint)4);
    EXPECT_EQ(pg.nodes[4]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[5]->node_id, (uint)5);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[5]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[6]->node_id, (uint)6);
    EXPECT_EQ(pg.nodes[6]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[6]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[7]->node_id, (uint)7);
    EXPECT_EQ(pg.nodes[7]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[7]->covg, (uint)1);
    EXPECT_EQ((uint)5, pg.edges.size());
    EXPECT_EQ(pg.edges[0]->from->node_id, (uint)0);
    EXPECT_EQ(pg.edges[0]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[0]->covg, (uint)2);
    EXPECT_EQ(pg.edges[1]->from->node_id, (uint)1);
    EXPECT_EQ(pg.edges[1]->to->node_id, (uint)2);
    EXPECT_EQ(pg.edges[1]->covg, (uint)1);
    EXPECT_EQ(pg.edges[2]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[2]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[2]->covg, (uint)1);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)2);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)3, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);
    EXPECT_EQ(pg.edges[2], pg.reads[0]->edges[2]);
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]);
    EXPECT_EQ(pg.edges[3], pg.reads[2]->edges[2]);
    EXPECT_EQ((uint)6, pg.reads[0]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[0].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[2].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[3].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[5].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[6].size());
    EXPECT_EQ((uint)3, pg.reads[1]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[3].size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[4].size());
    EXPECT_EQ((uint)5, pg.reads[2]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[0].size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[3].size());
    EXPECT_EQ((uint)1, pg.reads[2]->hits[4].size());
    EXPECT_EQ((uint)1, pg.reads[2]->hits[7].size());
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)2);
    EXPECT_EQ(pg.nodes[1]->reads.size(), (uint)2);
    EXPECT_EQ(pg.nodes[2]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[3]->reads.size(), (uint)3);
    EXPECT_EQ(pg.nodes[4]->reads.size(), (uint)1);
    EXPECT_EQ(pg.nodes[5]->reads.size(), (uint)1);
    EXPECT_EQ(pg.nodes[6]->reads.size(), (uint)1);
    EXPECT_EQ(pg.nodes[7]->reads.size(), (uint)1);

    //example where node appears twice in read
}*/

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
