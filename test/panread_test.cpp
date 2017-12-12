#include "gtest/gtest.h"
#include "pannode.h"
#include "panedge.h"
#include "panread.h"
#include "pangraph.h"
#include "minihit.h"
#include <stdint.h>
#include <iostream>
#include <set>

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
    set<MinimizerHitPtr, pComp> hits;
    pair<set<MinimizerHitPtr, pComp>::iterator,set<MinimizerHitPtr, pComp>::iterator> mhs(hits.begin(), hits.begin());

    pr1.add_hits(4, mhs);
    EXPECT_EQ((uint)1, pr1.hits.size());
    EXPECT_EQ((uint)0, pr1.hits[4].size());

    // add again and should append pr1.hits
    Interval i(0,5);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MinimizerHitPtr mh (make_shared<MinimizerHit>(4, i, 0, p, 0, 0));

    hits.insert(mh);
    mhs = make_pair(hits.begin(), --hits.end());
    pr1.add_hits(4, mhs);
    EXPECT_EQ((uint)1, pr1.hits.size());
    EXPECT_EQ((uint)1, pr1.hits[4].size());

    // add to new id
    pr1.add_hits(5, mhs);
    EXPECT_EQ((uint)2, pr1.hits.size());
    EXPECT_EQ((uint)1, pr1.hits[5].size());
    EXPECT_EQ((uint)1, pr1.hits[4].size());

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
    set<MinimizerHitPtr, pComp> hits;
    pair<set<MinimizerHitPtr, pComp>::iterator,set<MinimizerHitPtr, pComp>::iterator> mhs(hits.begin(), hits.begin());

    PanGraph pg;
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

    pg.add_edge(2,4,3);
    pg.edges.back()->covg -= 1;
    vector<PanEdge*>::iterator r = pg.reads[0]->replace_edge(pg.edges[2], pg.edges.back());
    // expect to get an iterator to the edge we have just inserted which now lies on read
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3); // note this is really now 2, but not been changed
    EXPECT_EQ(pg.nodes[4]->node_id, (uint)4);
    EXPECT_EQ(pg.nodes[4]->covg, (uint)2); // note this is really now 3
    EXPECT_EQ((uint)6, pg.edges.size());
    EXPECT_EQ(pg.edges[0]->from->node_id, (uint)0);
    EXPECT_EQ(pg.edges[0]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[0]->covg, (uint)2);
    EXPECT_EQ(pg.edges[1]->from->node_id, (uint)1);
    EXPECT_EQ(pg.edges[1]->to->node_id, (uint)2);
    EXPECT_EQ(pg.edges[1]->covg, (uint)1);
    EXPECT_EQ(pg.edges[2]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[2]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[2]->covg, (uint)0);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)2);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ(pg.edges[5]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[5]->to->node_id, (uint)4);
    EXPECT_EQ(pg.edges[5]->covg, (uint)1);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)3, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);    
    EXPECT_EQ(pg.edges[5], pg.reads[0]->edges[2]);
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]); 
    EXPECT_EQ(pg.edges[3], pg.reads[2]->edges[2]);    

    EXPECT_EQ(r, pg.reads[0]->edges.begin()+2);
    EXPECT_EQ(*r, pg.edges.back());

    pg.add_edge(3,2,3);
    pg.edges.back()->covg -= 1;
    r = pg.reads[2]->replace_edge(pg.edges[3], pg.edges.back());
    // expect to get an iterator to the edge we have just inserted which now lies on read
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1); // is really now 2
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3); // note this is really now 2, but not been changed
    EXPECT_EQ(pg.nodes[4]->node_id, (uint)4);
    EXPECT_EQ(pg.nodes[4]->covg, (uint)2); // note this is really now 2
    EXPECT_EQ((uint)7, pg.edges.size());
    EXPECT_EQ(pg.edges[0]->from->node_id, (uint)0);
    EXPECT_EQ(pg.edges[0]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[0]->covg, (uint)2);
    EXPECT_EQ(pg.edges[1]->from->node_id, (uint)1);
    EXPECT_EQ(pg.edges[1]->to->node_id, (uint)2);
    EXPECT_EQ(pg.edges[1]->covg, (uint)1);
    EXPECT_EQ(pg.edges[2]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[2]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[2]->covg, (uint)0);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)1);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ(pg.edges[5]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[5]->to->node_id, (uint)4);
    EXPECT_EQ(pg.edges[5]->covg, (uint)1);
    EXPECT_EQ(pg.edges[6]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[6]->to->node_id, (uint)2);
    EXPECT_EQ(pg.edges[6]->covg, (uint)1);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)3, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);
    EXPECT_EQ(pg.edges[5], pg.reads[0]->edges[2]);
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]);
    EXPECT_EQ(pg.edges[6], pg.reads[2]->edges[2]);  

    EXPECT_EQ(r, pg.reads[2]->edges.begin()+2);
    EXPECT_EQ(*r, pg.edges.back());                                                    
}

TEST_F(PanReadTest,replace_edge2)
{
    set<MinimizerHitPtr, pComp> hits;
    pair<set<MinimizerHitPtr, pComp>::iterator,set<MinimizerHitPtr, pComp>::iterator> mhs(hits.begin(), hits.begin());

    PanGraph pg;
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

    pg.add_edge(2,4,3);
    pg.edges.back()->covg -= 1;
    unordered_multiset<PanRead*>::iterator r = pg.reads[0]->replace_edge(pg.edges[2], pg.edges.back(), pg.edges[2]->reads.begin());
    // expect to get an iterator to the edge we have just inserted which now lies on read
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3); // note this is really now 2, but not been changed
    EXPECT_EQ(pg.nodes[4]->node_id, (uint)4);
    EXPECT_EQ(pg.nodes[4]->covg, (uint)2); // note this is really now 3
    EXPECT_EQ((uint)6, pg.edges.size());
    EXPECT_EQ(pg.edges[0]->from->node_id, (uint)0);
    EXPECT_EQ(pg.edges[0]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[0]->covg, (uint)2);
    EXPECT_EQ(pg.edges[1]->from->node_id, (uint)1);
    EXPECT_EQ(pg.edges[1]->to->node_id, (uint)2);
    EXPECT_EQ(pg.edges[1]->covg, (uint)1);
    EXPECT_EQ(pg.edges[2]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[2]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[2]->covg, (uint)0);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)2);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ(pg.edges[5]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[5]->to->node_id, (uint)4);
    EXPECT_EQ(pg.edges[5]->covg, (uint)1);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)3, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);
    EXPECT_EQ(pg.edges[5], pg.reads[0]->edges[2]);
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]);
    EXPECT_EQ(pg.edges[3], pg.reads[2]->edges[2]);

    EXPECT_EQ(r, pg.edges[2]->reads.end());

    pg.add_edge(3,0,3);
    pg.edges.back()->covg -= 1;
    r = pg.edges[3]->reads.find(pg.reads[2]);
    r = pg.reads[2]->replace_edge(pg.edges[3], pg.edges.back(), pg.edges[3]->reads.begin());
    // expect to get an iterator to the edge we have just inserted which now lies on read
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2); // is really now 3
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1); 
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3); // note this is really now 2, but not been changed
    EXPECT_EQ(pg.nodes[4]->node_id, (uint)4);
    EXPECT_EQ(pg.nodes[4]->covg, (uint)2); // note this is really now 2
    EXPECT_EQ((uint)7, pg.edges.size());
    EXPECT_EQ(pg.edges[0]->from->node_id, (uint)0);
    EXPECT_EQ(pg.edges[0]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[0]->covg, (uint)2);
    EXPECT_EQ(pg.edges[1]->from->node_id, (uint)1);
    EXPECT_EQ(pg.edges[1]->to->node_id, (uint)2);
    EXPECT_EQ(pg.edges[1]->covg, (uint)1);
    EXPECT_EQ(pg.edges[2]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[2]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[2]->covg, (uint)0);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)1);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ(pg.edges[5]->from->node_id, (uint)2);
    EXPECT_EQ(pg.edges[5]->to->node_id, (uint)4);
    EXPECT_EQ(pg.edges[5]->covg, (uint)1);
    EXPECT_EQ(pg.edges[6]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[6]->to->node_id, (uint)0);
    EXPECT_EQ(pg.edges[6]->covg, (uint)1);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)3, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);
    EXPECT_EQ(pg.edges[5], pg.reads[0]->edges[2]);
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]);
    EXPECT_EQ(pg.edges[6], pg.reads[2]->edges[2]);

    EXPECT_EQ(r, pg.edges[3]->reads.begin());
}

TEST_F(PanReadTest,remove_edge1)
{
    set<MinimizerHitPtr, pComp> hits;
    pair<set<MinimizerHitPtr, pComp>::iterator,set<MinimizerHitPtr, pComp>::iterator> mhs(hits.begin(), hits.begin());

    PanGraph pg;
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

    vector<PanEdge*>::iterator r = pg.reads[0]->remove_edge(pg.edges[2]);
    // expect to get an iterator to the edge after the one we deleted
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3); // note this is really now 2, but not been changed
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
    EXPECT_EQ(pg.edges[2]->covg, (uint)0);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)2);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)2, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);    
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]); 
    EXPECT_EQ(pg.edges[3], pg.reads[2]->edges[2]);    

    EXPECT_EQ(r, pg.reads[0]->edges.end());

    r = pg.reads[1]->remove_edge(pg.edges[3]);
    // expect to get an iterator to the edge we have just inserted which now lies on read
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3); // note this is really now 2, but not been changed
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
    EXPECT_EQ(pg.edges[2]->covg, (uint)0);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)1);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)2, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);
    EXPECT_EQ((uint)1, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]);
    EXPECT_EQ(pg.edges[3], pg.reads[2]->edges[2]);

    EXPECT_EQ(r, pg.reads[1]->edges.begin());
    EXPECT_EQ(*r, pg.edges[4]);
}

TEST_F(PanReadTest,remove_edge2)
{
    set<MinimizerHitPtr, pComp> hits;
    pair<set<MinimizerHitPtr, pComp>::iterator,set<MinimizerHitPtr, pComp>::iterator> mhs(hits.begin(), hits.begin());

    PanGraph pg;
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

    unordered_multiset<PanRead*>::iterator r = pg.reads[0]->remove_edge(pg.edges[2], pg.edges[2]->reads.begin());
    // expect to get an iterator to the edge we have just inserted which now lies on read
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3); // note this is really now 2, but not been changed
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
    EXPECT_EQ(pg.edges[2]->covg, (uint)0);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)2);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)2, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)3, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]);
    EXPECT_EQ(pg.edges[3], pg.reads[2]->edges[2]);

    EXPECT_EQ(r, pg.edges[2]->reads.end());

    r = pg.edges[3]->reads.find(pg.reads[2]);
    r = pg.reads[2]->remove_edge(pg.edges[3], r);
    // expect to get an iterator to the read after the one we removed
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1); 
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)3); // note this is really now 2, but not been changed
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
    EXPECT_EQ(pg.edges[2]->covg, (uint)0);
    EXPECT_EQ(pg.edges[3]->from->node_id, (uint)4);
    EXPECT_EQ(pg.edges[3]->to->node_id, (uint)3);
    EXPECT_EQ(pg.edges[3]->covg, (uint)1);
    EXPECT_EQ(pg.edges[4]->from->node_id, (uint)3);
    EXPECT_EQ(pg.edges[4]->to->node_id, (uint)1);
    EXPECT_EQ(pg.edges[4]->covg, (uint)2);
    EXPECT_EQ((uint)3, pg.reads.size());
    EXPECT_EQ((uint)2, pg.reads[0]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[0]->edges[0]);
    EXPECT_EQ(pg.edges[1], pg.reads[0]->edges[1]);
    EXPECT_EQ((uint)2, pg.reads[1]->edges.size());
    EXPECT_EQ(pg.edges[3], pg.reads[1]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[1]->edges[1]);
    EXPECT_EQ((uint)2, pg.reads[2]->edges.size());
    EXPECT_EQ(pg.edges[0], pg.reads[2]->edges[0]);
    EXPECT_EQ(pg.edges[4], pg.reads[2]->edges[1]);

    EXPECT_EQ((r==pg.edges[3]->reads.begin()) or (r==pg.edges[3]->reads.end()), true);
}

TEST_F(PanReadTest,replace_node)
{
    set<MinimizerHitPtr, pComp> hits;
    pair<set<MinimizerHitPtr, pComp>::iterator,set<MinimizerHitPtr, pComp>::iterator> mhs(hits.begin(), hits.begin());

    PanGraph pg;
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
    PanNode *n5;
    n5 = new PanNode(2, 5, "2");
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
    PanNode *n6;
    n6 = new PanNode(1, 6, "1");
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
    auto cl = make_pair(c.begin(), --c.end());
    pg.reads[2]->add_hits(4, cl);

    PanNode *n7;
    n7 = new PanNode(4, 7, "4");
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
}

TEST_F(PanReadTest,remove_node)
{
    set<MinimizerHitPtr, pComp> hits;
    pair<set<MinimizerHitPtr, pComp>::iterator,set<MinimizerHitPtr, pComp>::iterator> mhs(hits.begin(), hits.begin());

    PanGraph pg;
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
    pg.reads[0]->remove_node(pg.nodes[2]);
    
    EXPECT_EQ((uint)5, pg.nodes.size());
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
    EXPECT_EQ((uint)4, pg.reads[0]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[0].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[2].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[3].size());
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

    // example where old node appears in more than one read
    pg.reads[0]->remove_node(pg.nodes[1]);
    
    EXPECT_EQ((uint)5, pg.nodes.size());
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
    EXPECT_EQ((uint)4, pg.reads[0]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[0].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[2].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[3].size());
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

    // example where have actual hit
    Interval i(0,5);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MinimizerHitPtr mh (make_shared<MinimizerHit>(4, i, 0, p, 0, 0));
    set<MinimizerHitPtr, pComp> c;
    c.insert(mh);
    auto cl = make_pair(c.begin(), --c.end());

    pg.reads[2]->add_hits(4, cl);

    pg.reads[2]->remove_node(pg.nodes[4]);
    
    EXPECT_EQ((uint)5, pg.nodes.size());
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
    EXPECT_EQ((uint)4, pg.reads[0]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[0].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[2].size());
    EXPECT_EQ((uint)0, pg.reads[0]->hits[3].size());
    EXPECT_EQ((uint)3, pg.reads[1]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[3].size());
    EXPECT_EQ((uint)0, pg.reads[1]->hits[4].size());
    EXPECT_EQ((uint)4, pg.reads[2]->hits.size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[0].size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[1].size());
    EXPECT_EQ((uint)0, pg.reads[2]->hits[3].size());
    EXPECT_EQ((uint)1, pg.reads[2]->hits[4].size());
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)2);
    EXPECT_EQ(pg.nodes[1]->reads.size(), (uint)2);
    EXPECT_EQ(pg.nodes[2]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[3]->reads.size(), (uint)3);
    EXPECT_EQ(pg.nodes[4]->reads.size(), (uint)1);

    //example where node appears twice in read
}
