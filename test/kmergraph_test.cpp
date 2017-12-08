#include "gtest/gtest.h"
#include "test_macro.cpp"

#include "interval.h"
#include "path.h"
#include "kmergraph.h"
#include "kmernode.h"
#include <stdint.h>
#include <iostream>
#include <utility>

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

TEST_F(KmerGraphTest, add_node)
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
    EXPECT_EQ(j, kg.nodes[0]->num_AT);

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
    EXPECT_EQ(j, kg.nodes[0]->num_AT);

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
    EXPECT_EQ(j, kg.nodes[0]->num_AT);
    j = 1;
    EXPECT_EQ(j, kg.nodes[1]->id);
}

TEST_F(KmerGraphTest, add_node_with_kh)
{
    // add node and check it's there
    KmerGraph kg;

    deque<Interval> d = {Interval(0,3)};
    Path p;
    p.initialize(d);
    uint64_t kh = 469;
    kg.add_node_with_kh(p, kh);
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
    EXPECT_EQ(j, kg.nodes[0]->num_AT);

    EXPECT_EQ(kh, kg.nodes[0]->khash);
}

TEST_F(KmerGraphTest, add_edge)
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

    // now with kmernode constructor
    kg.add_node(p3);
    kg.add_edge(kg.nodes[0], kg.nodes[2]);
    j = 2;
    EXPECT_EQ(j, kg.nodes[0]->outNodes.size());
    j = 1;
    EXPECT_EQ(j, kg.nodes[1]->inNodes.size());
    EXPECT_EQ(j, kg.nodes[2]->inNodes.size());
    j = 0;
    EXPECT_EQ(j, kg.nodes[1]->outNodes.size());
    EXPECT_EQ(j, kg.nodes[0]->inNodes.size());

    // repeat and nothing should happen
    kg.add_edge(kg.nodes[0], kg.nodes[2]);
    j = 2;
    EXPECT_EQ(j, kg.nodes[0]->outNodes.size());
    j = 1;
    EXPECT_EQ(j, kg.nodes[1]->inNodes.size());
    EXPECT_EQ(j, kg.nodes[2]->inNodes.size());
    j = 0;
    EXPECT_EQ(j, kg.nodes[1]->outNodes.size());
    EXPECT_EQ(j, kg.nodes[0]->inNodes.size());
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

TEST_F(KmerGraphTest, copy)
{
    KmerGraph kg1;
    deque<Interval> d = {Interval(0,3)};
    Path p1,p2,p3;
    p1.initialize(d);
    kg1.add_node(p1);
    d = {Interval(1,4)};
    p2.initialize(d);
    kg1.add_node(p2);
    kg1.add_edge(p1,p2);

    KmerGraph kg2(kg1);

    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
}

TEST_F(KmerGraphTest, assign)
{
    KmerGraph kg1;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg1.add_node(p);
    d = {Interval(0,3)};
    Path p1,p2,p3;
    p1.initialize(d);
    kg1.add_node(p1);
    d = {Interval(1,4)};
    p2.initialize(d);
    kg1.add_node(p2);
    kg1.add_edge(p1,p2);
    d = {Interval(11,14)};
    p3.initialize(d);
    kg1.add_node(p3);
    kg1.add_edge(p1,p3);
    d = {Interval(15,18)};
    p1.initialize(d);
    kg1.add_node(p1);
    kg1.add_edge(p2,p1);
    d = {Interval(20,20)};
    p.initialize(d);
    kg1.add_node(p);
    kg1.add_edge(p1,p);
    kg1.add_edge(p3,p);

    KmerGraph kg2 = kg1;
    cout << "original " << kg1 << endl;

    cout << "copy " << kg2 << endl;

    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
}

TEST_F(KmerGraphTest, sort_topologically)
{
    vector<KmerNodePtr> exp_sorted_nodes;
    KmerNodePtr n;

    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = {Interval(24,24)};
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);

    /*for (uint i=0; i!=kg.nodes.size(); ++i)
    {
        cout << *kg.nodes[i] << endl;
    }*/

    kg.add_edge(kg.nodes[0],kg.nodes[1]);
    kg.add_edge(kg.nodes[0],kg.nodes[2]);
    kg.add_edge(kg.nodes[0],kg.nodes[3]);
    kg.add_edge(kg.nodes[1],kg.nodes[4]);
    kg.add_edge(kg.nodes[2],kg.nodes[5]);
    kg.add_edge(kg.nodes[3],kg.nodes[6]);
    kg.add_edge(kg.nodes[4],kg.nodes[6]);
    kg.add_edge(kg.nodes[5],kg.nodes[6]);

    kg.sort_topologically();
    /*for (uint i=0; i!=kg.sorted_nodes.size(); ++i)
    {
	cout << *kg.sorted_nodes[i] << endl;
    }*/

    // for each node, outnodes are further along vector
    vector<KmerNodePtr>::iterator it;
    uint i = 0;
    for (vector<KmerNodePtr>::iterator c=kg.sorted_nodes.begin(); c!=kg.sorted_nodes.end(); ++c)
    {
        for (auto d: (*c)->outNodes)
        {
	    it = c+1;
	    while ((*it)->path != d->path and it != kg.sorted_nodes.end())
	    {	
		it++;
            }
    	    EXPECT_EQ((it != kg.sorted_nodes.end()), true);
        }
	EXPECT_EQ(exp_sorted_nodes[i], *c);
	i++;
    }
}

TEST_F(KmerGraphTest, check)
{
    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);

    kg.add_edge(kg.nodes[0],kg.nodes[1]);
    kg.add_edge(kg.nodes[0],kg.nodes[2]);
    kg.add_edge(kg.nodes[0],kg.nodes[3]);
    kg.add_edge(kg.nodes[1],kg.nodes[4]);
    kg.add_edge(kg.nodes[2],kg.nodes[5]);
    kg.add_edge(kg.nodes[3],kg.nodes[6]);
    kg.add_edge(kg.nodes[4],kg.nodes[6]);
    kg.add_edge(kg.nodes[5],kg.nodes[6]);

    kg.sorted_nodes = {kg.nodes[0], kg.nodes[1], kg.nodes[2], kg.nodes[3], kg.nodes[4], kg.nodes[5], kg.nodes[6]};
    kg.check();
    kg.sorted_nodes = {kg.nodes[0], kg.nodes[1], kg.nodes[4], kg.nodes[3], kg.nodes[2], kg.nodes[5], kg.nodes[6]};
    kg.check();
    kg.sorted_nodes = {kg.nodes[6], kg.nodes[5], kg.nodes[0], kg.nodes[3], kg.nodes[2], kg.nodes[1], kg.nodes[4]};
    EXPECT_DEATH(kg.check(),"");
}

/*TEST_F(KmerGraphTest, get_prev)
{
    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);

    kg.add_edge(kg.nodes[0],kg.nodes[1]);
    kg.add_edge(kg.nodes[0],kg.nodes[2]);
    kg.add_edge(kg.nodes[0],kg.nodes[3]);
    kg.add_edge(kg.nodes[1],kg.nodes[4]);
    kg.add_edge(kg.nodes[2],kg.nodes[5]);
    kg.add_edge(kg.nodes[3],kg.nodes[6]);
    kg.add_edge(kg.nodes[4],kg.nodes[6]);
    kg.add_edge(kg.nodes[5],kg.nodes[6]);

    kg.covgs = {{{0,1,1,0,0,1,0},{0,0,0,0,0,0,0}},{{0,0,0,0,0,0,0},{0,0,1,1,1,0,0}}};

    // read 1 strand 0
    uint16_t prev;
    vector<deque<KmerNodePtr>> prev_paths, prev_paths_exp({{kg.nodes[0], kg.nodes[1]}});
    kg.get_prev(0,0,1,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[0], kg.nodes[2]}};
    kg.get_prev(0,0,2,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[0], kg.nodes[3]}};
    kg.get_prev(0,0,3,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[1], kg.nodes[4]}};
    kg.get_prev(0,0,4,prev, prev_paths);
    EXPECT_EQ((uint)1, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[2], kg.nodes[5]}};
    kg.get_prev(0,0,5,prev, prev_paths);
    EXPECT_EQ((uint)2, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[5], kg.nodes[6]}};
    kg.get_prev(0,0,6,prev, prev_paths);
    EXPECT_EQ((uint)5, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    // read 1 strand 1
    prev_paths.clear();
    prev_paths_exp ={{kg.nodes[0], kg.nodes[1]}};
    kg.get_prev(0,1,1,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[0], kg.nodes[2]}};
    kg.get_prev(0,1,2,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[0], kg.nodes[3]}};
    kg.get_prev(0,1,3,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[0], kg.nodes[1], kg.nodes[4]}};
    kg.get_prev(0,1,4,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[0], kg.nodes[2], kg.nodes[5]}};
    kg.get_prev(0,1,5,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[0], kg.nodes[3], kg.nodes[6]}, {kg.nodes[0], kg.nodes[1], kg.nodes[4], kg.nodes[6]}, {kg.nodes[0], kg.nodes[2], kg.nodes[5], kg.nodes[6]}};
    kg.get_prev(0,1,6,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    // read 2 strand 1
    prev_paths.clear();
    prev_paths_exp ={{kg.nodes[0], kg.nodes[1]}};
    kg.get_prev(1,1,1,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[0], kg.nodes[2]}};
    kg.get_prev(1,1,2,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[0], kg.nodes[3]}};
    kg.get_prev(1,1,3,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[0], kg.nodes[1], kg.nodes[4]}};
    kg.get_prev(1,1,4,prev, prev_paths);
    EXPECT_EQ((uint)0, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[2], kg.nodes[5]}};
    kg.get_prev(1,1,5,prev, prev_paths);
    EXPECT_EQ((uint)2, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);

    prev_paths.clear();
    prev_paths_exp = {{kg.nodes[3], kg.nodes[6]}, {kg.nodes[4], kg.nodes[6]}};
    kg.get_prev(1,1,6,prev, prev_paths);
    EXPECT_EQ((uint)4, prev);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, prev_paths_exp, prev_paths);
}*/

/*TEST_F(KmerGraphTest, get_next)
{
    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);

    kg.add_edge(kg.nodes[0],kg.nodes[1]);
    kg.add_edge(kg.nodes[0],kg.nodes[2]);
    kg.add_edge(kg.nodes[0],kg.nodes[3]);
    kg.add_edge(kg.nodes[1],kg.nodes[4]);
    kg.add_edge(kg.nodes[2],kg.nodes[5]);
    kg.add_edge(kg.nodes[3],kg.nodes[6]);
    kg.add_edge(kg.nodes[4],kg.nodes[6]);
    kg.add_edge(kg.nodes[5],kg.nodes[6]);

    kg.covgs = {{{0,1,1,0,0,1,0},{0,0,0,0,0,0,0}},{{0,0,0,0,0,0,0},{0,0,1,1,1,0,0}}};

    // read 1 strand 0
    uint16_t next;
    vector<deque<KmerNodePtr>> next_paths, next_paths_exp({{kg.nodes[0], kg.nodes[1]}, {kg.nodes[0], kg.nodes[2]}});
    kg.get_next(0,0,0,next, next_paths);
    EXPECT_EQ((uint)2, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[1], kg.nodes[4], kg.nodes[6]}};
    kg.get_next(0,0,1,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[2], kg.nodes[5]}};
    kg.get_next(0,0,2,next, next_paths);
    EXPECT_EQ((uint)5, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[3], kg.nodes[6]}};
    kg.get_next(0,0,3,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[4], kg.nodes[6]}};
    kg.get_next(0,0,4,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[5], kg.nodes[6]}};
    kg.get_next(0,0,5,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    // read 1 strand 1
    next_paths.clear();
    next_paths_exp = {{kg.nodes[0], kg.nodes[3], kg.nodes[6]}, {kg.nodes[0], kg.nodes[1], kg.nodes[4], kg.nodes[6]}, {kg.nodes[0], kg.nodes[2], kg.nodes[5], kg.nodes[6]}};
    kg.get_next(0,1,0,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[1], kg.nodes[4], kg.nodes[6]}};
    kg.get_next(0,1,1,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[2], kg.nodes[5], kg.nodes[6]}};
    kg.get_next(0,1,2,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[3], kg.nodes[6]}};
    kg.get_next(0,1,3,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[4], kg.nodes[6]}};
    kg.get_next(0,1,4,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[5], kg.nodes[6]}};
    kg.get_next(0,1,5,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    // read 2 strand 1
    next_paths.clear();
    next_paths_exp ={{kg.nodes[0], kg.nodes[2]}, {kg.nodes[0], kg.nodes[3]}};
    kg.get_next(1,1,0,next, next_paths);
    EXPECT_EQ((uint)3, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[1], kg.nodes[4]}};
    kg.get_next(1,1,1,next, next_paths);
    EXPECT_EQ((uint)4, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[2], kg.nodes[5], kg.nodes[6]}};
    kg.get_next(1,1,2,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[3], kg.nodes[6]}};
    kg.get_next(1,1,3,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[4], kg.nodes[6]}};
    kg.get_next(1,1,4,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);

    next_paths.clear();
    next_paths_exp = {{kg.nodes[5], kg.nodes[6]}};
    kg.get_next(1,1,5,next, next_paths);
    EXPECT_EQ((uint)6, next);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, next_paths_exp, next_paths);
}*/

/*TEST_F(KmerGraphTest, extend_paths_back)
{
    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);

    // note these are nonsense paths
    vector<deque<KmerNodePtr>> paths_to_extend = {{kg.nodes[0]}, {kg.nodes[0],kg.nodes[1]}};
    vector<deque<KmerNodePtr>> path_extensions = {{kg.nodes[4], kg.nodes[5],kg.nodes[0]}, {kg.nodes[6],kg.nodes[0]}};
    vector<deque<KmerNodePtr>> expected_result = {{kg.nodes[4], kg.nodes[5],kg.nodes[0]},
                                                  {kg.nodes[6], kg.nodes[0]},
                                                  {kg.nodes[4], kg.nodes[5], kg.nodes[0], kg.nodes[1]},
                                                  {kg.nodes[6], kg.nodes[0], kg.nodes[1]}};
    kg.extend_paths_back(paths_to_extend, path_extensions);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, expected_result, paths_to_extend);
}*/

TEST_F(KmerGraphTest, extend_paths_forward)
{
    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);

    // note these are nonsense paths
    vector<deque<KmerNodePtr>> paths_to_extend = {{kg.nodes[0], kg.nodes[1]}, {kg.nodes[1]}};
    vector<deque<KmerNodePtr>> path_extensions = {{kg.nodes[1], kg.nodes[4], kg.nodes[5]}, {kg.nodes[1], kg.nodes[6]}};
    vector<deque<KmerNodePtr>> expected_result = {{kg.nodes[0], kg.nodes[1], kg.nodes[4], kg.nodes[5]},
                                                  {kg.nodes[0],kg.nodes[1], kg.nodes[6]},
                                                  {kg.nodes[1], kg.nodes[4], kg.nodes[5]},
                                                  {kg.nodes[1],kg.nodes[6]}};
    kg.extend_paths_forward(paths_to_extend, path_extensions);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, expected_result, paths_to_extend);
}

TEST_F(KmerGraphTest,find_compatible_paths)
{
    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);

    kg.add_edge(kg.nodes[0],kg.nodes[1]);
    kg.add_edge(kg.nodes[0],kg.nodes[2]);
    kg.add_edge(kg.nodes[0],kg.nodes[3]);
    kg.add_edge(kg.nodes[1],kg.nodes[4]);
    kg.add_edge(kg.nodes[2],kg.nodes[5]);
    kg.add_edge(kg.nodes[3],kg.nodes[6]);
    kg.add_edge(kg.nodes[4],kg.nodes[6]);
    kg.add_edge(kg.nodes[5],kg.nodes[6]);

    kg.covgs = {{{0,1,1,0,0,1,0},{0,0,0,0,0,0,0}},{{0,0,0,0,0,0,0},{0,0,1,1,0,0,0}}};
    kg.nodes[1]->covg[0]+=1;
    kg.nodes[2]->covg[0]+=1;
    kg.nodes[5]->covg[0]+=1;
    kg.nodes[2]->covg[1]+=1;
    kg.nodes[3]->covg[1]+=1;

    vector<deque<KmerNodePtr>> paths;
    vector<deque<KmerNodePtr>> exp_paths = {{kg.nodes[0],kg.nodes[1],kg.nodes[4],kg.nodes[6]},
                                            {kg.nodes[0],kg.nodes[2],kg.nodes[5],kg.nodes[6]},
					    {kg.nodes[0],kg.nodes[3],kg.nodes[6]}};
    kg.find_compatible_paths(1,0,paths);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, exp_paths, paths);

    paths.clear();
    //exp_paths = {{kg.nodes[0],kg.nodes[2],kg.nodes[5],kg.nodes[6]}};
    kg.find_compatible_paths(1,1,paths);
    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, exp_paths, paths);

}

TEST_F(KmerGraphTest,find_all_compatible_paths)
{
    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);

    kg.add_edge(kg.nodes[0],kg.nodes[1]);
    kg.add_edge(kg.nodes[0],kg.nodes[2]);
    kg.add_edge(kg.nodes[0],kg.nodes[3]);
    kg.add_edge(kg.nodes[1],kg.nodes[4]);
    kg.add_edge(kg.nodes[2],kg.nodes[5]);
    kg.add_edge(kg.nodes[3],kg.nodes[6]);
    kg.add_edge(kg.nodes[4],kg.nodes[6]);
    kg.add_edge(kg.nodes[5],kg.nodes[6]);

    kg.covgs = {{{0,1,1,0,0,1,0},{0,0,0,0,0,0,0}},{{0,0,0,0,0,0,0},{0,0,1,0,1,0,0}}};
    kg.nodes[1]->covg[0]+=1;
    kg.nodes[2]->covg[0]+=1;
    kg.nodes[5]->covg[0]+=1;
    kg.nodes[2]->covg[1]+=1;
    kg.nodes[4]->covg[1]+=1;
    vector<deque<KmerNodePtr>> paths;
    vector<deque<KmerNodePtr>> exp_paths = {{kg.nodes[0],kg.nodes[1],kg.nodes[4],kg.nodes[6]},
                                        {kg.nodes[0],kg.nodes[2],kg.nodes[5],kg.nodes[6]}};
    vector<vector<pair<uint16_t,uint16_t>>> hit_pairs;
    vector<vector<pair<uint16_t,uint16_t>>> exp_hit_pairs = {{make_pair(0,1), make_pair(1,2)},
                                                             {make_pair(0,1), make_pair(1,1)}};
    kg.find_all_compatible_paths(paths, hit_pairs, 1,1);

    EXPECT_ITERABLE_EQ(vector<deque<KmerNodePtr>>, exp_paths, paths);
    EXPECT_EQ(exp_hit_pairs.size(), hit_pairs.size());
    EXPECT_EQ(exp_hit_pairs[0].size(), hit_pairs[0].size());
    EXPECT_EQ(exp_hit_pairs[0][0], hit_pairs[0][0]);
    EXPECT_EQ(exp_hit_pairs[0][1], hit_pairs[0][1]);
    EXPECT_EQ(exp_hit_pairs[1].size(), hit_pairs[1].size());
    EXPECT_EQ(exp_hit_pairs[1][0], hit_pairs[1][0]);
    EXPECT_EQ(exp_hit_pairs[1][1], hit_pairs[1][1]);

    //EXPECT_ITERABLE_EQ(vector<vector<pair<uint16_t,uint16_t>>>, exp_hit_pairs, hit_pairs);
}

TEST_F(KmerGraphTest, set_p)
{
    KmerGraph kg;
    EXPECT_DEATH(kg.set_p(0.4),"");
    kg.k = 3;
    EXPECT_DEATH(kg.set_p(0),"");
    EXPECT_DEATH(kg.set_p(1),"");
    kg.set_p(0.5);
    EXPECT_EQ(1/exp(1.5)-0.00001 <= kg.p and 1/exp(1.5)+0.00001 >= kg.p, true);
}

TEST_F(KmerGraphTest,prob)
{
    KmerGraph kg;
    EXPECT_DEATH(kg.prob(0),"");  // out of range, no nodes have been defined
 
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);

    EXPECT_DEATH(kg.prob(0),""); // still no p set   
    kg.k = 3;
    kg.set_p(0.5);
    
    EXPECT_DEATH(kg.prob(0),""); // no num_reads set
    kg.num_reads = 1;

    EXPECT_EQ(kg.nodes.size(), (uint)1);
    EXPECT_EQ(0, kg.prob(0));

    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);

    EXPECT_EQ(kg.nodes.size(), (uint)3);
    cout << kg.prob(1) << endl;
    cout << kg.prob(2) << endl;
    //EXPECT_EQ(0, prob(1));

    EXPECT_EQ(kg.prob(1), kg.prob(1,1));
    EXPECT_EQ(kg.prob(2), kg.prob(2,1));

    cout << kg.prob(1,1) << endl;
    cout << kg.prob(2,1) << endl;
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

    kg.num_reads = 5;
    kg.k = 3;

    vector<KmerNodePtr> mp;
    kg.set_p(0.01);
    kg.find_max_path(mp);
    vector<KmerNodePtr> exp_order = {kg.nodes[1], kg.nodes[2]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kg.nodes[1]->covg[0]-=4;
    kg.nodes[2]->covg[0]-=3;
    kg.nodes[5]->covg[1]+=5;
    kg.set_p(0.01);
    kg.find_max_path(mp);
    exp_order = {kg.nodes[5]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);
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

    kg.num_reads = 5;
    kg.k = 3;

    vector<KmerNodePtr> mp;
    kg.set_p(0.01);
    kg.find_max_path(mp);
    vector<KmerNodePtr> exp_order = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kg.nodes[4]->covg[0]-=4;
    kg.nodes[5]->covg[0]-=3;
    kg.nodes[6]->covg[0]-=5;
    kg.nodes[7]->covg[0]-=4;
    kg.nodes[8]->covg[1]+=5;
    kg.set_p(0.01);
    kg.find_max_path(mp);
    exp_order = {kg.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);
}


TEST_F(KmerGraphTest,find_max_paths_2Level)
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

    kg.num_reads = 10;
    kg.k = 3;

    kg.set_p(0.01);
    vector<vector<KmerNodePtr>> mps = kg.find_max_paths(2);
    EXPECT_EQ((uint)2, mps.size());
    vector<KmerNodePtr> exp_order = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mps[1]);

    exp_order = {kg.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mps[0]);
}

TEST_F(KmerGraphTest,random_paths)
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

    vector<vector<KmerNodePtr>> rps;
    vector<KmerNodePtr> exp_order1 = {kg.nodes[1], kg.nodes[2], kg.nodes[3], kg.nodes[7]};
    vector<KmerNodePtr> exp_order2 = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7]};
    vector<KmerNodePtr> exp_order3 = {kg.nodes[8]};

    rps = kg.get_random_paths(10);
    for (uint i=0; i!=rps.size(); ++i)
    {
	//cout << "new path of length " << rps.size() << ": ";
        for (uint j=0; j!= rps[i].size(); ++j)
	{
	    //cout << rps[i][j]->id << "->";
	    if (rps[i][j]->id == 1)
	    {
		EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order1, rps[i]);
	    } else if (rps[i][j]->id == 4)
	    {
		EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order2, rps[i]);
	    } else if (rps[i][j]->id == 8)
            {
                EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order3, rps[i]);
	    }
	}
	//cout << endl;
    }
}

TEST_F(KmerGraphTest,path_prob)
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

    kg.num_reads = 5;
    kg.k = 3;

    vector<KmerNodePtr> mp;
    kg.set_p(0.01);
    float mp_p = kg.find_max_path(mp);
    vector<KmerNodePtr> exp_order = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7], kg.nodes[9]};
    float exp_p = 0;
    for (uint i=0; i!=exp_order.size(); ++i)
    {
	exp_p += kg.prob(exp_order[i]->id);
    }
    exp_p /= 4;
    EXPECT_EQ(mp_p, exp_p);

    mp.clear();
    kg.nodes[4]->covg[0]-=4;
    kg.nodes[5]->covg[0]-=3;
    kg.nodes[6]->covg[0]-=5;
    kg.nodes[7]->covg[0]-=4;
    kg.nodes[8]->covg[1]+=5;
    kg.set_p(0.01);
    mp_p = kg.find_max_path(mp);
    exp_order = {kg.nodes[8], kg.nodes[9]};
    exp_p = 0;
    for (uint i=0; i!=exp_order.size(); ++i)
    {
        exp_p += kg.prob(exp_order[i]->id);
    }
    EXPECT_EQ(mp_p, exp_p);
}


TEST_F(KmerGraphTest,path_probs)
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

    kg.num_reads = 10;
    kg.k = 3;
    kg.set_p(0.01);
    vector<vector<KmerNodePtr>> mps = kg.find_max_paths(2);
    EXPECT_EQ((uint)2, mps.size());

    // check get right answer
    vector<KmerNodePtr> exp_nodes = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7], kg.nodes[8]};
    float exp_p = 0;
    for (uint i=0; i!=exp_nodes.size(); ++i)
    {
        exp_p += kg.prob(exp_nodes[i]->id, 5);
    }
    cout << exp_p << "/5 = ";
    exp_p /= 5;
    cout << exp_p << endl;
    EXPECT_EQ(kg.prob_paths(mps), exp_p);

}

TEST_F(KmerGraphTest, save_covg_dist){
    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p,p1,p2;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,3)};
    p1.initialize(d);
    kg.add_node(p1);
    d = {Interval(1,4)};
    p2.initialize(d);
    kg.add_node(p2);
    kg.add_edge(p1,p2);
    d = {Interval(4,4)};
    p.initialize(d);
    kg.add_node(p);
    kg.nodes[1]->covg[1] +=5;
    kg.nodes[2]->covg[1] +=4;
    kg.nodes[1]->num_AT = 4;
    kg.nodes[2]->num_AT = 6;

    kg.save_covg_dist("../test/test_cases/kmergraph_test.covg.txt");
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
    EXPECT_EQ((uint)0, kg.nodes[0]->num_AT);

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

