#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "de_bruijn/node.h"
#include <iostream>
#include <deque>
#include <unordered_set>


using namespace debruijn;

TEST(DeBruijnNodeTest, create) {
    deque<uint_least32_t> v({4, 6, 8});
    unordered_multiset<uint32_t> w({0});
    Node n(2, v, 0);
    EXPECT_EQ(n.id, (uint) 2);
    EXPECT_ITERABLE_EQ(deque<uint_least32_t>, n.hashed_node_ids, v);
    EXPECT_ITERABLE_EQ(unordered_multiset<uint32_t>, n.read_ids, w);
}

TEST(DeBruijnNodeTest, equals) {
    deque<uint_least32_t> v({4, 7, 8});
    deque<uint_least32_t> w({4, 6, 8});
    deque<uint_least32_t> y({9, 6, 5});

    Node n1(2, v, 0);
    Node n2(2, v, 5);
    Node n3(3, v, 0);
    Node n4(2, w, 0);
    Node n5(2, y, 0);


    EXPECT_EQ(n1, n1);
    EXPECT_EQ(n2, n2);
    EXPECT_EQ(n3, n3);
    EXPECT_EQ(n4, n4);
    EXPECT_EQ(n5, n5);

    EXPECT_EQ(n1, n2);
    EXPECT_EQ(n1, n3);
    EXPECT_EQ(n2, n1);
    EXPECT_EQ(n3, n1);
    EXPECT_EQ(n2, n3);
    EXPECT_EQ(n3, n2);
    EXPECT_EQ(n1, n5);
    EXPECT_EQ(n2, n5);
    EXPECT_EQ(n3, n5);
    EXPECT_EQ(n5, n1);
    EXPECT_EQ(n5, n2);
    EXPECT_EQ(n5, n3);


    EXPECT_NE(n1, n4);
    EXPECT_NE(n4, n1);
    EXPECT_NE(n2, n4);
    EXPECT_NE(n4, n2);
    EXPECT_NE(n3, n4);
    EXPECT_NE(n4, n3);
    EXPECT_NE(n5, n4);
    EXPECT_NE(n4, n5);
}

