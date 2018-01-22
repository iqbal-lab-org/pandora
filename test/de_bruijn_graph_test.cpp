#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "de_bruijn_graph_class.h"
#include "de_bruijn/node.h"
#include <iostream>
#include <deque>
#include <unordered_set>

using namespace debruijn;

class DeBruijnGraphTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(DeBruijnGraphTest,create)
{
    GraphTester g(5);
    EXPECT_EQ(g.size, (uint)5);
    EXPECT_EQ(g.next_id, (uint)0);
}

TEST_F(DeBruijnGraphTest,add_node)
{
    GraphTester g(3);

    deque<uint16_t> v({4,6,8});
    unordered_multiset<uint32_t> w({0});

    g.add_node(v, 0);

    EXPECT_EQ(g.nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[0], Node(0, v, 0));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w);

    // add same node
    g.add_node(v, 0);
    w.insert(0);

    EXPECT_EQ(g.nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[0], Node(0, v, 0));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w);

    // add same node different read
    g.add_node(v, 7);
    w.insert(7);

    EXPECT_EQ(g.nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[0], Node(0, v, 7));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w);

    // add different node same read
    v = {6,9,3};
    g.add_node(v, 7);
    w.erase(0);

    EXPECT_EQ(g.nodes.size(), (uint)2);
    EXPECT_EQ(*g.nodes[1], Node(1, v, 7));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[1]->hashed_node_ids, v);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[1]->read_ids, w);
}

TEST_F(DeBruijnGraphTest,add_edge)
{
    GraphTester g(3);

    deque<uint16_t> v1({4,6,8});
    deque<uint16_t> v2({5,6,9});
    deque<uint16_t> v3({6,2,2});
    deque<uint16_t> v4({5,8,3});
    NodePtr n1 = g.add_node(v1, 0);
    NodePtr n2 = g.add_node(v2, 0);
    g.add_edge(n1, n2);

    EXPECT_EQ(g.nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[0]->out_nodes.begin(), (uint)1);
    EXPECT_EQ(*g.nodes[1]->out_nodes.begin(), (uint)0);

    // add same edge
    g.add_edge(n1, n2);

    EXPECT_EQ(g.nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[0]->out_nodes.begin(), (uint)1);
    EXPECT_EQ(*g.nodes[1]->out_nodes.begin(), (uint)0);

    // add edge the other way
    g.add_edge(n2, n1);

    EXPECT_EQ(g.nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[0]->out_nodes.begin(), (uint)1);
    EXPECT_EQ(*g.nodes[1]->out_nodes.begin(), (uint)0);

    // add a different edge
    NodePtr n3 = g.add_node(v3, 0);
    g.add_edge(n1, n3);

    EXPECT_EQ(g.nodes.size(), (uint)3);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->out_nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[1]->out_nodes.begin(), (uint)0);
    EXPECT_EQ(*g.nodes[2]->out_nodes.begin(), (uint)0);

    // add a different read edge
    NodePtr n4 = g.add_node(v4, 5);
    g.add_edge(n3, n4);

    EXPECT_EQ(g.nodes.size(), (uint)4);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[3]->out_nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[1]->out_nodes.begin(), (uint)0);
    EXPECT_EQ(*g.nodes[3]->out_nodes.begin(), (uint)2);
}

TEST_F(DeBruijnGraphTest,remove_node)
{
    GraphTester g(3);
    deque<uint16_t> v1({4,6,8});
    deque<uint16_t> v2({6,9,3});

    g.add_node(v1, 0);
    NodePtr n1 = g.add_node(v1, 7);
    NodePtr n2 = g.add_node(v2, 7);
    g.add_edge(n1, n2);
    
    EXPECT_EQ(g.nodes.size(), (uint)2);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[1]->hashed_node_ids, v2);
    unordered_multiset<uint32_t> w1({0,7});
    unordered_multiset<uint32_t> w2({7});
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids,w1);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[1]->read_ids,w2);
    unordered_set<uint16_t> s({1});
    unordered_set<uint16_t> t({0});
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[0]->out_nodes,s);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[1]->out_nodes,t);

    // remove a node
    g.remove_node(1);
    EXPECT_EQ(g.nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids,w1);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(),(uint)0);
}


TEST_F(DeBruijnGraphTest,remove_read_from_node)
{
    GraphTester g(3);
    deque<uint16_t> v1({4,6,8});
    deque<uint16_t> v2({6,9,3});
    deque<uint16_t> v3({1,2,3});

    g.add_node(v1, 0);
    g.add_node(v2, 4);
    g.add_node(v3, 5);
    NodePtr n1 = g.add_node(v1, 7);
    NodePtr n2 = g.add_node(v2, 7);
    g.add_edge(n1, n2);

    EXPECT_EQ(g.nodes.size(), (uint)3);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_EQ(*g.nodes[2], Node(3, v3, 5));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[1]->hashed_node_ids, v2);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[2]->hashed_node_ids, v3);
    unordered_multiset<uint32_t> w1({0,7});
    unordered_multiset<uint32_t> w2({4,7});
    unordered_multiset<uint32_t> w3({5});
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids,w1);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[1]->read_ids,w2);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[2]->read_ids,w3);
    unordered_set<uint16_t> s({1});
    unordered_set<uint16_t> t({0});
    unordered_set<uint16_t> u;
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[0]->out_nodes,s);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[1]->out_nodes,t);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[2]->out_nodes,u);

    // remove a read which doesn't exist - nothing should happen
    g.remove_read_from_node(1,0);
    EXPECT_EQ(g.nodes.size(), (uint)3);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_EQ(*g.nodes[2], Node(3, v3, 5));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[1]->hashed_node_ids, v2);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[2]->hashed_node_ids, v3);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids,w1);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[1]->read_ids,w2);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[2]->read_ids,w3);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[0]->out_nodes,s);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[1]->out_nodes,t);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[2]->out_nodes,u);

    // remove a read from a node which doesn't exist - nothing should happen
    g.remove_read_from_node(0,3);
    EXPECT_EQ(g.nodes.size(), (uint)3);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_EQ(*g.nodes[2], Node(3, v3, 5));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[1]->hashed_node_ids, v2);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[2]->hashed_node_ids, v3);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids,w1);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[1]->read_ids,w2);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[2]->read_ids,w3);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[0]->out_nodes,s);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[1]->out_nodes,t);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[2]->out_nodes,u);

    // remove read from a node where should just change the read id list for node
    g.remove_read_from_node(7,1);
    EXPECT_EQ(g.nodes.size(), (uint)3);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_EQ(*g.nodes[2], Node(3, v3, 5));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[1]->hashed_node_ids, v2);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[2]->hashed_node_ids, v3);
    w2 = {4};
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids,w1);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[1]->read_ids,w2);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[2]->read_ids,w3);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[0]->out_nodes,u);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[1]->out_nodes,u);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[2]->out_nodes,u);

    // remove read from a node where should result in node being removed
    g.remove_read_from_node(5,2);
    EXPECT_EQ(g.nodes.size(), (uint)2);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[1]->hashed_node_ids, v2);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids,w1);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[1]->read_ids,w2);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[0]->out_nodes,u);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[1]->out_nodes,u);

    // continue removing reads until graph empty
    g.remove_read_from_node(0,0);
    EXPECT_EQ(g.nodes.size(), (uint)2);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[1]->hashed_node_ids, v2);
    w1 = {7};
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids,w1);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[1]->read_ids,w2);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[0]->out_nodes,u);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[1]->out_nodes,u);

    g.remove_read_from_node(4,1);
    EXPECT_EQ(g.nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, g.nodes[0]->read_ids,w1);
    EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, g.nodes[0]->out_nodes,u);

    g.remove_read_from_node(7,0);
    EXPECT_EQ(g.nodes.size(), (uint)0);
}


TEST_F(DeBruijnGraphTest,get_leaves)
{
    GraphTester g(3);

    deque<uint16_t> v1({4,6,8});
    deque<uint16_t> v2({5,6,9});
    deque<uint16_t> v3({6,2,2});
    deque<uint16_t> v4({5,8,3});
    deque<uint16_t> v5({5,9,3});

    NodePtr n1 = g.add_node(v1, 0);
    NodePtr n2 = g.add_node(v2, 0);
    g.add_edge(n1, n2);
    NodePtr n3 = g.add_node(v3, 0);
    g.add_edge(n1, n3);
    NodePtr n4 = g.add_node(v4, 5);
    g.add_edge(n3, n4);
    NodePtr n5 = g.add_node(v5, 5);

    EXPECT_EQ(g.nodes.size(), (uint)5);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[3]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[4]->out_nodes.size(), (uint)0);

    unordered_set<uint16_t> l = g.get_leaves();
    unordered_set<uint16_t> l_exp = {1, 3, 4};
    for (auto i : l_exp)
    {
        EXPECT_EQ(l.find(i)!=l.end(), true);
    }
    //EXPECT_ITERABLE_EQ(unordered_set<uint16_t>, l, l_exp);
}

TEST_F(DeBruijnGraphTest, get_leaves2)
{
    Graph dbg_exp(3);
    deque<uint16_t> d = {0,2,4}; //0
    NodePtr n1 = dbg_exp.add_node(d, 0);
    d = {2,4,6};//1
    NodePtr n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1,n2);
    d = {4,6,8};//2
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2,n1);
    d = {6,8,10};//3
    n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1,n2);

    d = {6,8,10};//3
    n2 = dbg_exp.add_node(d, 1);
    d = {8,10,0};//4
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2,n1);
    d = {10,0,2};//5
    n2 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n1,n2);
    d = {0,2,4};//0
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2,n1);

    d = {2,4,6};//1
    n1 = dbg_exp.add_node(d, 2);
    d = {4,6,14};//6
    n2 = dbg_exp.add_node(d, 2);
    dbg_exp.add_edge(n1,n2);

    d = {0,12,6};//7
    n1 = dbg_exp.add_node(d, 3);
    d = {12,6,8};//8
    n2 = dbg_exp.add_node(d, 3);
    dbg_exp.add_edge(n1,n2);

    d = {0,2,4};//0
    n1 = dbg_exp.add_node(d, 4);
    d = {2,4,12};//9
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1,n2);
    d = {4,12,6};//10
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2,n1);
    d = {12,6,8};//8
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1,n2);
    d = {6,8,10};//3
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2,n1);

    d = {12,2,4};//11
    n1 = dbg_exp.add_node(d, 5);
    d = {2,4,12};//9
    n2 = dbg_exp.add_node(d, 5);
    dbg_exp.add_edge(n1,n2);
    d = {4,12,6};//10
    n1 = dbg_exp.add_node(d, 5);
    dbg_exp.add_edge(n2,n1);

    unordered_set<uint16_t> l = dbg_exp.get_leaves();
    unordered_set<uint16_t> l_exp = {6,7,11};
    for (auto i : l_exp)
    {
        EXPECT_EQ(l.find(i)!=l.end(), true);
    }
}

TEST_F(DeBruijnGraphTest,get_unitigs)
{
    // 0 -> 1
    //   \> 2 -> 3
    // 4

    GraphTester g(3);

    deque<uint16_t> v1({4,6,8});
    deque<uint16_t> v2({5,6,9});
    deque<uint16_t> v3({6,2,2});
    deque<uint16_t> v4({5,8,3});
    deque<uint16_t> v5({5,9,3});

    NodePtr n0 = g.add_node(v1, 0);
    NodePtr n1 = g.add_node(v2, 0);
    g.add_edge(n0, n1);
    NodePtr n2 = g.add_node(v3, 0);
    g.add_edge(n0, n2);
    NodePtr n3 = g.add_node(v4, 5);
    g.add_edge(n2, n3);
    NodePtr n4 = g.add_node(v5, 5);

    EXPECT_EQ(g.nodes.size(), (uint)5);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[3]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[4]->out_nodes.size(), (uint)0);

    set<deque<uint16_t>> s = g.get_unitigs();
    deque<uint16_t> d1({1,0,2,3});
    deque<uint16_t> d2({3,2,0,1});
    set<deque<uint16_t>> s_exp1({d1});
    set<deque<uint16_t>> s_exp2({d2});
    d1 = {4};
    s_exp1.insert(d1);
    s_exp2.insert(d1);
    EXPECT_EQ(s.size(), s_exp1.size());
    EXPECT_EQ(s.size(), s_exp2.size());
    auto dit1 = s_exp1.begin();
    auto dit2 = s_exp2.begin();
    for (auto dit : s)
    {
        EXPECT_EQ(( (dit1!=s_exp1.end() and dit==*dit1) or
                    (dit2!=s_exp2.end() and dit==*dit2))      , true);
        if (dit1!=s_exp1.end() and dit==*dit1)
        {
            dit1++;
            dit2 = s_exp2.end();
        } else if (dit2!=s_exp2.end() and dit==*dit2)
        {
            dit2++;
            dit1 = s_exp1.end();
        } else {
            break;
        }
    }
}

TEST_F(DeBruijnGraphTest,extend_unitig)
{
    // 0 -> 1
    //   \> 2 -> 3
    // 4

    GraphTester g(3);

    deque<uint16_t> v1({4,6,8});
    deque<uint16_t> v2({5,6,9});
    deque<uint16_t> v3({6,2,2});
    deque<uint16_t> v4({5,8,3});
    deque<uint16_t> v5({5,9,3});

    NodePtr n1 = g.add_node(v1, 0);
    NodePtr n2 = g.add_node(v2, 0);
    g.add_edge(n1, n2);
    NodePtr n3 = g.add_node(v3, 0);
    g.add_edge(n1, n3);
    NodePtr n4 = g.add_node(v4, 5);
    g.add_edge(n3, n4);
    NodePtr n5 = g.add_node(v5, 5);

    EXPECT_EQ(g.nodes.size(), (uint)5);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[3]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[4]->out_nodes.size(), (uint)0);

    deque<uint16_t> d = {0};
    g.extend_unitig(d);
    deque<uint16_t> d_exp = {1,0,2,3};
    deque<uint16_t> d_exp1 = {3,2,0,1};
    EXPECT_EQ((d == d_exp or d == d_exp1), true);
    //EXPECT_ITERABLE_EQ(deque<uint16_t>, d, d_exp);

    d = {1};
    g.extend_unitig(d);
    EXPECT_EQ((d == d_exp or d == d_exp1), true);
    //EXPECT_ITERABLE_EQ(deque<uint16_t>, d, d_exp);

    d = {2};
    g.extend_unitig(d);
    EXPECT_EQ((d == d_exp or d == d_exp1), true);
    //EXPECT_ITERABLE_EQ(deque<uint16_t>, d, d_exp);

    d = {3};
    g.extend_unitig(d);
    EXPECT_EQ((d == d_exp or d == d_exp1), true);
    //EXPECT_ITERABLE_EQ(deque<uint16_t>, d, d_exp);

    d = {4};
    g.extend_unitig(d);
    d_exp = {4};
    EXPECT_ITERABLE_EQ(deque<uint16_t>, d, d_exp);

    // and check doesn't break if there is a cycle
    g.nodes.clear();
    g.next_id = 0;
    v1 = {0,1,2};
    v2 = {1,2,3};
    v3 = {2,3,4};
    v4 = {3,4,5};
    v5 = {4,5,0};
    deque<uint16_t> v6({5,0,1});

    n1 = g.add_node(v1, 0);
    n2 = g.add_node(v2, 0);
    g.add_edge(n1, n2);
    n3 = g.add_node(v3, 0);
    g.add_edge(n2, n3);
    n4 = g.add_node(v4, 0);
    g.add_edge(n3, n4);
    n5 = g.add_node(v5, 0);
    g.add_edge(n4, n5);
    NodePtr n6 = g.add_node(v6, 0);
    g.add_edge(n5, n6);
    g.add_edge(n6, n1);

    EXPECT_EQ(g.nodes.size(), (uint)6);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[2]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[3]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[4]->out_nodes.size(), (uint)2);

    d = {1};
    g.extend_unitig(d);
    d_exp = {0,1,2,3,4,5,0};
    d_exp1 = {2,1,0,5,4,3,2};
    EXPECT_EQ((d == d_exp or d == d_exp1),true);
}

TEST_F(DeBruijnGraphTest,equals)
{
    GraphTester g1(3);

    deque<uint16_t> v1({4,6,8});
    deque<uint16_t> v2({5,6,9});
    deque<uint16_t> v3({6,2,2});
    deque<uint16_t> v4({5,8,3});
    deque<uint16_t> v5({5,9,3});

    NodePtr n1 = g1.add_node(v1, 0);
    NodePtr n2 = g1.add_node(v2, 0);
    g1.add_edge(n1, n2);
    NodePtr n3 = g1.add_node(v3, 0);
    g1.add_edge(n1, n3);
    NodePtr n4 = g1.add_node(v4, 5);
    g1.add_edge(n3, n4);
    NodePtr n5 = g1.add_node(v5, 5);

    GraphTester g2(3);

    NodePtr m2 = g2.add_node(v2, 0);
    EXPECT_NE(g1, g2);
    NodePtr m3 = g2.add_node(v3, 0);
    EXPECT_NE(g1, g2);
    NodePtr m5 = g2.add_node(v5, 5);
    EXPECT_NE(g1, g2);
    NodePtr m4 = g2.add_node(v4, 5);
    EXPECT_NE(g1, g2);
    g2.add_edge(m3, m4);
    EXPECT_NE(g1, g2);
    NodePtr m1 = g2.add_node(v1, 0);
    EXPECT_NE(g1, g2);
    g2.add_edge(m1, m2);
    EXPECT_NE(g1, g2);
    g2.add_edge(m1, m3);

    // shouldn't matter that nodes and edges added in different order
    EXPECT_EQ(g1, g2);
    EXPECT_EQ(g2, g1);


    // an extra node does matter
    deque<uint16_t> v6({0,0,3});
    NodePtr m6 = g2.add_node(v6, 0);
    EXPECT_NE(g1, g2);
    EXPECT_NE(g2, g1);
    g2.remove_node(5);

    // an extra edge does matter
    g2.add_edge(m5, m3);
    EXPECT_NE(g1, g2);
    EXPECT_NE(g2, g1);

}

