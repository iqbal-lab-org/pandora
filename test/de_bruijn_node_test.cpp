#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "de_bruijn/node.h"
#include <iostream>
#include <deque>
#include <unordered_set>

using namespace debruijn;

class DeBruijnNodeTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(DeBruijnNodeTest,create)
{
    deque<uint16_t> v({4,6,8});
    unordered_multiset<uint32_t> w({0});
    Node n(2, v, 0);
    EXPECT_EQ(n.id, (uint)2);
    EXPECT_ITERABLE_EQ(deque<uint16_t>, n.hashed_node_ids, v);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, n.read_ids, w);
}

TEST_F(DeBruijnNodeTest,equals)
{
    deque<uint16_t> v({4,6,8});
    deque<uint16_t> w({4,7,8});
    Node n1(2, v, 0);
    Node n2(2, v, 5);
    Node n3(3, v, 0);
    Node n4(2, w, 0);

    EXPECT_EQ(n1, n1);
    EXPECT_EQ(n2, n2); 
    EXPECT_EQ(n3, n3);
    EXPECT_EQ(n4, n4);

    EXPECT_EQ(n1, n2);
    EXPECT_EQ(n1, n3);
    EXPECT_EQ(n2, n1);
    EXPECT_EQ(n3, n1);
    EXPECT_EQ(n2, n3);
    EXPECT_EQ(n3, n2);

    EXPECT_NE(n1, n4);
    EXPECT_NE(n4, n1);
    EXPECT_NE(n2, n4);
    EXPECT_NE(n4, n2);
    EXPECT_NE(n3, n4);
    EXPECT_NE(n4, n3);
}

