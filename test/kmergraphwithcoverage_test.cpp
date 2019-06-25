#include "gtest/gtest.h"
#include "test_macro.cpp"

#include "interval.h"
#include "prg/path.h"
#include "kmergraphwithcoverage.h"
#include "kmernode.h"
#include "localPRG.h"
#include <stdint.h>
#include <iostream>
#include <cmath>


//using namespace std;

/*TEST(KmerGraphWithCoverageTest, equals) {
    KmerGraphWithCoverage kg1, kg2;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kg1.add_node(p1);
    auto m1 = kg2.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg1.add_node(p2);
    auto m2 = kg2.add_node(p2);
    kg1.add_edge(n1, n2);
    kg2.add_edge(m1, m2);

    d = {Interval(2, 5)};
    p3.initialize(d);
    auto m3 = kg2.add_node(p3);

    // same as themselves, different if different numbers of nodes
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1 == kg2), false);
    EXPECT_EQ((kg2 == kg1), false);

    auto n3 = kg1.add_node(p3);
    kg2.add_edge(m1, m3);

    // same as themselves, different if different numbers of edges
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1 == kg2), false);
    EXPECT_EQ((kg2 == kg1), false);

    kg1.add_edge(n2, n3);

    // same as themselves, different if edges in different places
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1 == kg2), false);
    EXPECT_EQ((kg2 == kg1), false);
}

TEST(KmerGraphWithCoverageTest, copy) {
    KmerGraphWithCoverage kg1;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kg1.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg1.add_node(p2);
    kg1.add_edge(n1, n2);

    KmerGraphWithCoverage kg2(kg1);

    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
}

TEST(KmerGraphWithCoverageTest, assign) {
    KmerGraphWithCoverage kg1;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    auto n = kg1.add_node(p);
    d = {Interval(0, 3)};
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kg1.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg1.add_node(p2);
    kg1.add_edge(n1, n2);
    d = {Interval(11, 14)};
    p3.initialize(d);
    auto n3 = kg1.add_node(p3);
    kg1.add_edge(n1, n3);
    d = {Interval(15, 18)};
    p1.initialize(d);
    n1 = kg1.add_node(p1);
    kg1.add_edge(n2, n1);
    d = {Interval(20, 20)};
    p.initialize(d);
    n = kg1.add_node(p);
    kg1.add_edge(n1, n);
    kg1.add_edge(n3, n);

    KmerGraphWithCoverage kg2 = kg1;

    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
}

TEST(KmerGraphWithCoverageTest, save) {
    auto l = std::make_shared<LocalPRG>(LocalPRG(1, "test localPRG", "ACGT"));

    KmerGraphWithCoverage kg;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);
    EXPECT_EQ((uint) 0, kg.nodes[0]->num_AT);

    //kg.save("../test/test_cases/kmergraph_test.gfa");
    kg.save("kmergraph_test.gfa", l);
}

TEST(KmerGraphWithCoverageTest, save_no_prg) {
    KmerGraphWithCoverage kg;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);
    EXPECT_EQ((uint) 0, kg.nodes[0]->num_AT);

    //kg.save("../test/test_cases/kmergraph_test.gfa");
    kg.save("kmergraph_test2.gfa");
}

TEST(KmerGraphWithCoverageTest, load) {
    KmerGraphWithCoverage kg, read_kg;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);

    //read_kg.load("../test/test_cases/kmergraph_test.gfa");
    read_kg.load("kmergraph_test2.gfa");
    EXPECT_EQ(kg, read_kg);
}

TEST(KmerGraphWithCoverageTest, load_prg) {
    KmerGraphWithCoverage read_kg;
    EXPECT_DEATH(read_kg.load("kmergraph_test.gfa"), "");
}*/
