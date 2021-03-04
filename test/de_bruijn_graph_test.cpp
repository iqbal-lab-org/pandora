#include <iostream>
#include <deque>
#include <unordered_set>
#include <unordered_map>

#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "de_bruijn_graph_class.h"
#include "test_helpers.h"

using namespace debruijn;

TEST(DeBruijnGraphCreate, Initialize_SetsSizeAndNextId)
{
    GraphTester g(5);
    EXPECT_EQ(g.size, (uint)5);
    EXPECT_EQ(g.next_id, (uint)0);
}

TEST(DeBruijnGraphAddNode, AddNode_NodeHashInIndex)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v = { 4, 6, 8 };
    uint32_t read_id = 0;
    g.add_node(v, read_id);

    bool found = g.node_hash.find(v) != g.node_hash.end();
    EXPECT_TRUE(found);
}

TEST(DeBruijnGraphAddNode, AddNode_NodeIdInIndex)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v = { 4, 6, 8 };
    uint32_t read_id = 0;
    g.add_node(v, read_id);

    bool found = g.nodes.find(0) != g.nodes.end();
    EXPECT_TRUE(found);
}

TEST(DeBruijnGraphAddNode, AddNode_NodePropertiesCorrect)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v = { 4, 6, 8 };
    uint32_t read_id = 0;
    std::unordered_multiset<uint32_t> w = { read_id };
    g.add_node(v, read_id);

    EXPECT_EQ(*g.nodes[0], Node(0, v, 0));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w);
}

TEST(DeBruijnGraphAddNode, AddNodeTwiceForSameRead_NodeReadsMultisetContainsReadTwice)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v = { 4, 6, 8 };
    uint32_t read_id = 0;
    std::unordered_multiset<uint32_t> w = { read_id, read_id };
    g.add_node(v, read_id);
    g.add_node(v, read_id);

    EXPECT_EQ(*g.nodes[0], Node(0, v, 0));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w);
}

TEST(DeBruijnGraphAddNode, AddNodeTwiceForDifferentRead_NodeReadsMultisetContainsReads)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v = { 4, 6, 8 };
    uint32_t read_id = 0;
    g.add_node(v, read_id);
    std::unordered_multiset<uint32_t> w = { read_id };
    read_id = 7;
    g.add_node(v, read_id);
    w.insert(read_id);

    EXPECT_EQ(*g.nodes[0], Node(0, v, 0));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w);
}

TEST(DeBruijnGraphAddNode, AddTwoNodes_SecondNodeInIndex)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v = { 4, 6, 8 };
    uint32_t read_id = 0;
    g.add_node(v, read_id);

    v = { 6, 9, 3 };
    read_id = 7;
    g.add_node(v, read_id);

    bool found = g.nodes.find(1) != g.nodes.end();
    EXPECT_TRUE(found);
}

TEST(DeBruijnGraphAddNode, AddTwoNodes_SecondNodePropertiesCorrect)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v = { 4, 6, 8 };
    uint32_t read_id = 0;
    g.add_node(v, read_id);

    v = { 6, 9, 3 };
    read_id = 7;
    g.add_node(v, read_id);
    std::unordered_multiset<uint32_t> w = { read_id };

    EXPECT_EQ(*g.nodes[1], Node(1, v, 7));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[1]->hashed_node_ids, v);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[1]->read_ids, w);
}

TEST(DeBruijnGraphAddEdge, AddEdgeNodesOverlapForwards_EdgeAdded)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v1 = { 4, 6, 8 };
    std::deque<uint_least32_t> v2 = { 6, 8, 9 };
    OrientedNodePtr n1 = g.add_node(v1, 0);
    OrientedNodePtr n2 = g.add_node(v2, 0);
    g.add_edge(n1, n2);

    bool found_outnode_n1
        = n1.first->out_nodes.find(n2.first->id) != n1.first->out_nodes.end();
    bool found_innode_n2
        = n2.first->in_nodes.find(n1.first->id) != n2.first->in_nodes.end();
    bool found_innode_n1
        = n1.first->in_nodes.find(n2.first->id) != n1.first->in_nodes.end();
    bool found_outnode_n2
        = n2.first->out_nodes.find(n1.first->id) != n2.first->out_nodes.end();
    EXPECT_TRUE(found_outnode_n1);
    EXPECT_TRUE(found_innode_n2);
    EXPECT_FALSE(found_innode_n1);
    EXPECT_FALSE(found_outnode_n2);
}

TEST(DeBruijnGraphAddEdge, AddEdgeFirstForwardSecondRC_EdgeAdded)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v1 = { 4, 6, 8 };
    std::deque<uint_least32_t> v3 = { 6, 8, 9 };
    std::deque<uint_least32_t> v2 = { 8, 9, 7 };
    OrientedNodePtr n1 = g.add_node(v1, 0);
    OrientedNodePtr n2 = g.add_node(v2, 0);
    OrientedNodePtr n3 = g.add_node(v3, 0);
    g.add_edge(n1, n3);

    bool found_outnode_n1
        = n1.first->out_nodes.find(n3.first->id) != n1.first->out_nodes.end();
    bool found_innode_n3
        = n3.first->in_nodes.find(n1.first->id) != n3.first->in_nodes.end();
    bool found_innode_n1
        = n1.first->in_nodes.find(n3.first->id) != n1.first->in_nodes.end();
    bool found_outnode_n3
        = n3.first->out_nodes.find(n1.first->id) != n3.first->out_nodes.end();
    EXPECT_TRUE(found_outnode_n1);
    EXPECT_FALSE(found_innode_n3);
    EXPECT_FALSE(found_innode_n1);
    EXPECT_TRUE(found_outnode_n3);
}

TEST(DeBruijnGraphAddEdge, AddEdgeFirstRCSecondForward_EdgeAdded)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v1 = { 9, 7, 5 };
    std::deque<uint_least32_t> v2 = { 4, 6, 8 };
    std::deque<uint_least32_t> v3 = { 6, 8, 9 };
    OrientedNodePtr n1 = g.add_node(v1, 0);
    OrientedNodePtr n2 = g.add_node(v2, 0);
    OrientedNodePtr n3 = g.add_node(v3, 0);
    g.add_edge(n2, n3);

    bool found_outnode_n2
        = n2.first->out_nodes.find(n3.first->id) != n2.first->out_nodes.end();
    bool found_innode_n3
        = n3.first->in_nodes.find(n2.first->id) != n3.first->in_nodes.end();
    bool found_innode_n2
        = n2.first->in_nodes.find(n3.first->id) != n2.first->in_nodes.end();
    bool found_outnode_n3
        = n3.first->out_nodes.find(n2.first->id) != n3.first->out_nodes.end();
    EXPECT_FALSE(found_outnode_n2);
    EXPECT_TRUE(found_innode_n3);
    EXPECT_TRUE(found_innode_n2);
    EXPECT_FALSE(found_outnode_n3);
}

TEST(DeBruijnGraphAddEdge, AddEdgeNodesBothRC_EdgeAdded)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v1 = { 4, 6, 8 };
    std::deque<uint_least32_t> v2 = { 6, 8, 9 };
    std::deque<uint_least32_t> v1_ = { 9, 7, 5 };
    std::deque<uint_least32_t> v2_ = { 8, 9, 7 };
    OrientedNodePtr n1 = g.add_node(v1, 0);
    OrientedNodePtr n2 = g.add_node(v2, 0);
    n1 = g.add_node(v1_, 0);
    n2 = g.add_node(v2_, 0);
    g.add_edge(n2, n1);

    bool found_outnode_n1
        = n1.first->out_nodes.find(n2.first->id) != n1.first->out_nodes.end();
    bool found_innode_n2
        = n2.first->in_nodes.find(n1.first->id) != n2.first->in_nodes.end();
    bool found_innode_n1
        = n1.first->in_nodes.find(n2.first->id) != n1.first->in_nodes.end();
    bool found_outnode_n2
        = n2.first->out_nodes.find(n1.first->id) != n2.first->out_nodes.end();
    EXPECT_TRUE(found_outnode_n1);
    EXPECT_TRUE(found_innode_n2);
    EXPECT_FALSE(found_innode_n1);
    EXPECT_FALSE(found_outnode_n2);
}

TEST(DeBruijnGraphAddEdge, AddEdgeNoOverlap_FatalRuntimeError)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v1 = { 4, 6, 8 };
    std::deque<uint_least32_t> v2 = { 6, 0, 9 };
    OrientedNodePtr n1 = g.add_node(v1, 0);
    OrientedNodePtr n2 = g.add_node(v2, 0);

    ASSERT_EXCEPTION(g.add_edge(n1, n2), FatalRuntimeError, "Error adding edge to de Bruijn Graph");
}

TEST(DeBruijnGraphTest, remove_node)
{
    GraphTester g(3);
    std::deque<uint_least32_t> v1 = { 4, 6, 8 };
    std::deque<uint_least32_t> v2 = { 6, 8, 3 };

    g.add_node(v1, 0);
    OrientedNodePtr n1 = g.add_node(v1, 7);
    OrientedNodePtr n2 = g.add_node(v2, 7);
    g.add_edge(n1, n2);

    EXPECT_EQ(g.nodes.size(), (uint)2);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[1]->hashed_node_ids, v2);
    std::unordered_multiset<uint32_t> w1 = { 0, 7 };
    std::unordered_multiset<uint32_t> w2 = { 7 };
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w1);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[1]->read_ids, w2);
    std::unordered_set<uint32_t> s = { 1 };
    std::unordered_set<uint32_t> t = { 0 };
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->out_nodes, s);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->in_nodes, t);

    // remove a node
    g.remove_node(1);
    EXPECT_EQ(g.nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w1);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)0);
}

TEST(DeBruijnGraphAddEdge, AddEdgeTwice_EdgeAddedOnce)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v1 = { 4, 6, 8 };
    std::deque<uint_least32_t> v2 = { 6, 8, 9 };
    OrientedNodePtr n1 = g.add_node(v1, 0);
    OrientedNodePtr n2 = g.add_node(v2, 0);
    g.add_edge(n1, n2);

    EXPECT_EQ(n1.first->out_nodes.size(), (uint)1);
    EXPECT_EQ(n2.first->out_nodes.size(), (uint)0);
    EXPECT_EQ(n1.first->in_nodes.size(), (uint)0);
    EXPECT_EQ(n2.first->in_nodes.size(), (uint)1);
}

TEST(DeBruijnGraphTest, remove_read_from_node)
{
    GraphTester g(3);
    std::deque<uint_least32_t> v1 = { 4, 6, 8 };
    std::deque<uint_least32_t> v2 = { 6, 8, 3 };
    std::deque<uint_least32_t> v3 = { 1, 2, 3 };

    g.add_node(v1, 0);
    g.add_node(v2, 4);
    g.add_node(v3, 5);
    OrientedNodePtr n1 = g.add_node(v1, 7);
    OrientedNodePtr n2 = g.add_node(v2, 7);
    g.add_edge(n1, n2);

    EXPECT_EQ(g.nodes.size(), (uint)3);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_EQ(*g.nodes[2], Node(3, v3, 5));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[1]->hashed_node_ids, v2);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[2]->hashed_node_ids, v3);
    std::unordered_multiset<uint32_t> w1 = { 0, 7 };
    std::unordered_multiset<uint32_t> w2 = { 4, 7 };
    std::unordered_multiset<uint32_t> w3 = { 5 };
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w1);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[1]->read_ids, w2);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[2]->read_ids, w3);
    std::unordered_set<uint32_t> s = { 1 };
    std::unordered_set<uint32_t> t = { 0 };
    std::unordered_set<uint32_t> u;
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->out_nodes, s);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->in_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->in_nodes, t);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[2]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[2]->in_nodes, u);

    // remove a read which doesn't exist - nothing should happen
    g.remove_read_from_node(1, 0);
    EXPECT_EQ(g.nodes.size(), (uint)3);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_EQ(*g.nodes[2], Node(3, v3, 5));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[1]->hashed_node_ids, v2);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[2]->hashed_node_ids, v3);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w1);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[1]->read_ids, w2);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[2]->read_ids, w3);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->out_nodes, s);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->in_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->in_nodes, t);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[2]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[2]->in_nodes, u);

    // remove a read from a node which doesn't exist - nothing should happen
    g.remove_read_from_node(0, 3);
    EXPECT_EQ(g.nodes.size(), (uint)3);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_EQ(*g.nodes[2], Node(3, v3, 5));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[1]->hashed_node_ids, v2);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[2]->hashed_node_ids, v3);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w1);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[1]->read_ids, w2);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[2]->read_ids, w3);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->out_nodes, s);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->in_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->in_nodes, t);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[2]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[2]->in_nodes, u);

    // remove read from a node where should just change the read id list for node
    g.remove_read_from_node(7, 1);
    EXPECT_EQ(g.nodes.size(), (uint)3);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_EQ(*g.nodes[2], Node(3, v3, 5));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[1]->hashed_node_ids, v2);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[2]->hashed_node_ids, v3);
    w2 = { 4 };
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w1);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[1]->read_ids, w2);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[2]->read_ids, w3);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->in_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->in_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[2]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[2]->in_nodes, u);

    // remove read from a node where should result in node being removed
    g.remove_read_from_node(5, 2);
    EXPECT_EQ(g.nodes.size(), (uint)2);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[1]->hashed_node_ids, v2);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w1);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[1]->read_ids, w2);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->in_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->in_nodes, u);

    // continue removing reads until graph empty
    g.remove_read_from_node(0, 0);
    EXPECT_EQ(g.nodes.size(), (uint)2);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_EQ(*g.nodes[1], Node(1, v2, 7));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[1]->hashed_node_ids, v2);
    w1 = { 7 };
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w1);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[1]->read_ids, w2);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->in_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[1]->in_nodes, u);

    g.remove_read_from_node(4, 1);
    EXPECT_EQ(g.nodes.size(), (uint)1);
    EXPECT_EQ(*g.nodes[0], Node(0, v1, 7));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, g.nodes[0]->hashed_node_ids, v1);
    EXPECT_ITERABLE_EQ(std::unordered_multiset<uint32_t>, g.nodes[0]->read_ids, w1);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->out_nodes, u);
    EXPECT_ITERABLE_EQ(std::unordered_set<uint32_t>, g.nodes[0]->in_nodes, u);

    g.remove_read_from_node(7, 0);
    EXPECT_EQ(g.nodes.size(), (uint)0);
}

TEST(DeBruijnGraphTest, get_leaves)
{
    GraphTester g(3);

    std::deque<uint_least32_t> v1 = { 4, 1, 8 };
    std::deque<uint_least32_t> v2 = { 1, 8, 9 };
    std::deque<uint_least32_t> v3 = { 1, 8, 2 };
    std::deque<uint_least32_t> v4 = { 8, 2, 4 };
    std::deque<uint_least32_t> v5 = { 2, 4, 3 };

    OrientedNodePtr n1 = g.add_node(v1, 0);
    OrientedNodePtr n2 = g.add_node(v2, 0);
    g.add_edge(n1, n2);
    OrientedNodePtr n3 = g.add_node(v3, 0);
    g.add_edge(n1, n3);
    OrientedNodePtr n4 = g.add_node(v4, 5);
    g.add_edge(n3, n4);
    OrientedNodePtr n5 = g.add_node(v5, 5);

    std::unordered_set<uint32_t> l = g.get_leaves();
    std::unordered_set<uint32_t> l_exp = { 1, 3, 4 };
    for (const auto& i : l_exp) {
        EXPECT_EQ(l.find(i) != l.end(), true);
    }
    // EXPECT_ITERABLE_EQ(std::unordered_set<uint_least32_t>, l, l_exp);
}

TEST(DeBruijnGraphTest, get_leaves2)
{
    Graph dbg_exp(3);
    std::deque<uint_least32_t> d = { 0, 2, 4 }; // 0
    OrientedNodePtr n1 = dbg_exp.add_node(d, 0);
    d = { 2, 4, 6 }; // 1
    OrientedNodePtr n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1, n2);
    d = { 4, 6, 8 }; // 2
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2, n1);
    d = { 6, 8, 10 }; // 3
    n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1, n2);

    d = { 6, 8, 10 }; // 3
    n2 = dbg_exp.add_node(d, 1);
    d = { 8, 10, 0 }; // 4
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2, n1);
    d = { 10, 0, 2 }; // 5
    n2 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n1, n2);
    d = { 0, 2, 4 }; // 0
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2, n1);

    d = { 2, 4, 6 }; // 1
    n1 = dbg_exp.add_node(d, 2);
    d = { 4, 6, 14 }; // 6
    n2 = dbg_exp.add_node(d, 2);
    dbg_exp.add_edge(n1, n2);

    d = { 0, 12, 6 }; // 7
    n1 = dbg_exp.add_node(d, 3);
    d = { 12, 6, 8 }; // 8
    n2 = dbg_exp.add_node(d, 3);
    dbg_exp.add_edge(n1, n2);

    d = { 0, 2, 4 }; // 0
    n1 = dbg_exp.add_node(d, 4);
    d = { 2, 4, 12 }; // 9
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1, n2);
    d = { 4, 12, 6 }; // 10
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);
    d = { 12, 6, 8 }; // 8
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1, n2);
    d = { 6, 8, 10 }; // 3
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);

    d = { 12, 2, 4 }; // 11
    n1 = dbg_exp.add_node(d, 5);
    d = { 2, 4, 12 }; // 9
    n2 = dbg_exp.add_node(d, 5);
    dbg_exp.add_edge(n1, n2);
    d = { 4, 12, 6 }; // 10
    n1 = dbg_exp.add_node(d, 5);
    dbg_exp.add_edge(n2, n1);

    std::unordered_set<uint32_t> l = dbg_exp.get_leaves();
    std::unordered_set<uint32_t> l_exp = { 6, 7, 11 };
    for (const auto& i : l_exp) {
        EXPECT_EQ(l.find(i) != l.end(), true);
    }
}

TEST(DeBruijnGraphGetUnitigs, OneBubble_ThreeTigs)
{
    // 0 -> 1 -> 2 ------> 3 -> 4 -> 5 -> 0
    //             \> 6 -/

    // 012->123->234------> 345->450
    //    \>126->263->634/>

    GraphTester g(3);

    std::deque<uint_least32_t> v1 = { 0, 2, 4 };
    std::deque<uint_least32_t> v2 = { 2, 4, 6 };
    std::deque<uint_least32_t> v3 = { 4, 6, 8 };
    std::deque<uint_least32_t> v4 = { 6, 8, 10 };
    std::deque<uint_least32_t> v5 = { 8, 10, 0 };

    OrientedNodePtr n0 = g.add_node(v1, 0);
    OrientedNodePtr n1 = g.add_node(v2, 0);
    g.add_edge(n0, n1);
    OrientedNodePtr n2 = g.add_node(v3, 0);
    g.add_edge(n1, n2);
    OrientedNodePtr n3 = g.add_node(v4, 0);
    g.add_edge(n2, n3);
    OrientedNodePtr n4 = g.add_node(v5, 0);
    g.add_edge(n3, n4);

    std::deque<uint_least32_t> v6 = { 2, 4, 12 };
    std::deque<uint_least32_t> v7 = { 4, 12, 6 };
    std::deque<uint_least32_t> v8 = { 12, 6, 8 };

    n0 = g.add_node(v1, 1);
    n1 = g.add_node(v6, 1);
    g.add_edge(n0, n1);
    n2 = g.add_node(v7, 1);
    g.add_edge(n1, n2);
    n3 = g.add_node(v8, 1);
    g.add_edge(n2, n3);
    n4 = g.add_node(v4, 1);
    g.add_edge(n3, n4);

    std::set<std::deque<uint32_t>> s = g.get_unitigs();
    EXPECT_EQ(s.size(), (uint)3);

    std::set<std::deque<uint32_t>> s_exp;
    std::deque<uint32_t> d = { 0, 1, 2, 3 };
    s_exp.insert(d);
    d = { 0, 5, 6, 7, 3 };
    s_exp.insert(d);
    d = { 3, 4 };
    s_exp.insert(d);

    EXPECT_ITERABLE_EQ(std::set<std::deque<uint32_t>>, s, s_exp);
}

TEST(DeBruijnGraphTest, get_unitigs)
{
    // 0 -> 1
    //   \> 2 -> 3
    // 4

    GraphTester g(3);

    std::deque<uint_least32_t> v1 = { 4, 6, 8 };
    std::deque<uint_least32_t> v2 = { 6, 8, 9 };
    std::deque<uint_least32_t> v3 = { 6, 8, 2 };
    std::deque<uint_least32_t> v4 = { 8, 2, 3 };
    std::deque<uint_least32_t> v5 = { 5, 9, 3 };

    OrientedNodePtr n0 = g.add_node(v1, 0);
    OrientedNodePtr n1 = g.add_node(v2, 0);
    g.add_edge(n0, n1);
    OrientedNodePtr n2 = g.add_node(v3, 0);
    g.add_edge(n0, n2);
    OrientedNodePtr n3 = g.add_node(v4, 5);
    g.add_edge(n2, n3);
    OrientedNodePtr n4 = g.add_node(v5, 5);

    EXPECT_EQ(g.nodes.size(), (uint)5);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[0]->in_nodes.size(), (uint)0);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)0);
    EXPECT_EQ(g.nodes[1]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[3]->out_nodes.size(), (uint)0);
    EXPECT_EQ(g.nodes[3]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[4]->out_nodes.size(), (uint)0);
    EXPECT_EQ(g.nodes[4]->in_nodes.size(), (uint)0);

    std::set<std::deque<uint32_t>> s = g.get_unitigs();
    std::deque<uint32_t> d1 = { 0, 2, 3 };
    std::deque<uint32_t> d2 = { 0, 1 };
    std::set<std::deque<uint32_t>> s_exp = { d1, d2 };
    d1 = { 4 };
    s_exp.insert(d1);
    EXPECT_EQ(s.size(), s_exp.size());
    EXPECT_ITERABLE_EQ(std::set<std::deque<uint32_t>>, s, s_exp);
}

TEST(DeBruijnGraphTest, extend_unitig)
{
    // 0 -> 1
    //   \> 2 -> 3
    // 4

    GraphTester g(3);

    std::deque<uint_least32_t> v1 = { 4, 6, 8 };
    std::deque<uint_least32_t> v2 = { 6, 8, 9 };
    std::deque<uint_least32_t> v3 = { 6, 8, 2 };
    std::deque<uint_least32_t> v4 = { 8, 2, 3 };
    std::deque<uint_least32_t> v5 = { 5, 9, 3 };

    OrientedNodePtr n1 = g.add_node(v1, 0);
    OrientedNodePtr n2 = g.add_node(v2, 0);
    g.add_edge(n1, n2);
    OrientedNodePtr n3 = g.add_node(v3, 0);
    g.add_edge(n1, n3);
    OrientedNodePtr n4 = g.add_node(v4, 5);
    g.add_edge(n3, n4);
    OrientedNodePtr n5 = g.add_node(v5, 5);

    EXPECT_EQ(g.nodes.size(), (uint)5);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)2);
    EXPECT_EQ(g.nodes[0]->in_nodes.size(), (uint)0);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)0);
    EXPECT_EQ(g.nodes[1]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[3]->out_nodes.size(), (uint)0);
    EXPECT_EQ(g.nodes[3]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[4]->out_nodes.size(), (uint)0);
    EXPECT_EQ(g.nodes[4]->in_nodes.size(), (uint)0);

    std::deque<uint32_t> d = { 0 };
    g.extend_unitig(d);
    std::deque<uint32_t> d_exp = { 0 };
    EXPECT_ITERABLE_EQ(std::deque<uint32_t>, d, d_exp);

    d = { 1 };
    g.extend_unitig(d);
    d_exp = { 0, 1 };
    EXPECT_ITERABLE_EQ(std::deque<uint32_t>, d, d_exp);

    d = { 2 };
    d_exp = { 0, 2, 3 };
    g.extend_unitig(d);
    EXPECT_ITERABLE_EQ(std::deque<uint32_t>, d, d_exp);

    d = { 3 };
    g.extend_unitig(d);
    EXPECT_ITERABLE_EQ(std::deque<uint32_t>, d, d_exp);

    d = { 4 };
    g.extend_unitig(d);
    d_exp = { 4 };
    EXPECT_ITERABLE_EQ(std::deque<uint32_t>, d, d_exp);

    // and check doesn't break if there is a cycle
    g.nodes.clear();
    g.next_id = 0;
    v1 = { 0, 1, 2 };
    v2 = { 1, 2, 3 };
    v3 = { 2, 3, 4 };
    v4 = { 3, 4, 5 };
    v5 = { 4, 5, 0 };
    std::deque<uint_least32_t> v6 = { 5, 0, 1 };

    n1 = g.add_node(v1, 0);
    n2 = g.add_node(v2, 0);
    g.add_edge(n1, n2);
    n3 = g.add_node(v3, 0);
    g.add_edge(n2, n3);
    n4 = g.add_node(v4, 0);
    g.add_edge(n3, n4);
    n5 = g.add_node(v5, 0);
    g.add_edge(n4, n5);
    OrientedNodePtr n6 = g.add_node(v6, 0);
    g.add_edge(n5, n6);
    g.add_edge(n6, n1);

    EXPECT_EQ(g.nodes.size(), (uint)6);
    EXPECT_EQ(g.nodes[0]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[0]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[1]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[1]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[2]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[3]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[3]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[4]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[4]->in_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[5]->out_nodes.size(), (uint)1);
    EXPECT_EQ(g.nodes[5]->in_nodes.size(), (uint)1);

    d = { 1 };
    g.extend_unitig(d);
    d_exp = { 1, 2, 3, 4, 5, 0 };
    EXPECT_ITERABLE_EQ(std::deque<uint32_t>, d, d_exp);
}

TEST(DeBruijnGraphTest, equals)
{
    GraphTester g1(3);

    std::deque<uint_least32_t> v1 = { 4, 6, 8 };
    std::deque<uint_least32_t> v2 = { 6, 8, 9 };
    std::deque<uint_least32_t> v3 = { 6, 8, 2 };
    std::deque<uint_least32_t> v4 = { 8, 2, 3 };
    std::deque<uint_least32_t> v5 = { 5, 6, 8 };

    OrientedNodePtr n1 = g1.add_node(v1, 0);
    OrientedNodePtr n2 = g1.add_node(v2, 0);
    g1.add_edge(n1, n2);
    OrientedNodePtr n3 = g1.add_node(v3, 0);
    g1.add_edge(n1, n3);
    OrientedNodePtr n4 = g1.add_node(v4, 5);
    g1.add_edge(n3, n4);
    OrientedNodePtr n5 = g1.add_node(v5, 5);

    GraphTester g2(3);

    OrientedNodePtr m2 = g2.add_node(v2, 0);
    EXPECT_NE(g1, g2);
    OrientedNodePtr m3 = g2.add_node(v3, 0);
    EXPECT_NE(g1, g2);
    OrientedNodePtr m5 = g2.add_node(v5, 5);
    EXPECT_NE(g1, g2);
    OrientedNodePtr m4 = g2.add_node(v4, 5);
    EXPECT_NE(g1, g2);
    g2.add_edge(m3, m4);
    EXPECT_NE(g1, g2);
    OrientedNodePtr m1 = g2.add_node(v1, 0);
    EXPECT_NE(g1, g2);
    g2.add_edge(m1, m2);
    EXPECT_NE(g1, g2);
    g2.add_edge(m1, m3);

    // shouldn't matter that nodes and edges added in different order
    EXPECT_EQ(g1, g2);
    EXPECT_EQ(g2, g1);

    // an extra node does matter
    std::deque<uint_least32_t> v6 = { 0, 0, 3 };
    OrientedNodePtr m6 = g2.add_node(v6, 0);

    EXPECT_NE(g1, g2);
    EXPECT_NE(g2, g1);

    g2.remove_node(5);
    EXPECT_EQ(g1, g2);
    EXPECT_EQ(g2, g1);

    // an extra edge does matter
    g2.add_edge(m5, m3);
    EXPECT_NE(g1, g2);
    EXPECT_NE(g2, g1);
}
