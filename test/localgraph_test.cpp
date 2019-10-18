#include "interval.h"
#include "localPRG.h"
#include "localgraph.h"
#include "localnode.h"
#include "prg/path.h"
#include "test_macro.cpp"
#include "gtest/gtest.h"
#include <iostream>
#include <stdint.h>

using namespace std;

TEST(LocalGraphTest, add_node)
{
    // add node and check it's there
    LocalGraph lg1;
    lg1.add_node(0, "AGCT", Interval(0, 4));
    LocalNode ln1("AGCT", Interval(0, 4), 0);
    EXPECT_EQ(ln1, *lg1.nodes[0]);

    // add node another time and expect nothing to happen
    lg1.add_node(0, "AGCT", Interval(0, 4));
    EXPECT_EQ(ln1, *lg1.nodes[0]);

    // add impossible nodes and expect and error
    EXPECT_DEATH(lg1.add_node(0, "AGGT", Interval(0, 4)), "");
    EXPECT_DEATH(lg1.add_node(1, "AGG", Interval(0, 4)), "");
}

TEST(LocalGraphTest, add_edge)
{
    LocalGraph lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(7, 8));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    EXPECT_EQ(lg2.nodes[0]->outNodes[0], lg2.nodes[1]);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);

    // expect failure if a node doesn't exist in the graph
    EXPECT_DEATH(lg2.add_edge(0, 4), "");
}

TEST(LocalGraphTest, equals)
{
    LocalGraph lg1;
    lg1.add_node(0, "AGCT", Interval(0, 4));
    EXPECT_EQ(lg1, lg1);

    LocalGraph lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(7, 8));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);
    EXPECT_EQ(lg2, lg2);

    EXPECT_EQ((lg1 == lg2), false);

    // order adding shouldn't matter
    LocalGraph lg2p;
    lg2p.add_node(2, "G", Interval(7, 8));
    lg2p.add_node(3, "T", Interval(13, 14));
    lg2p.add_node(1, "GC", Interval(4, 6));
    lg2p.add_node(0, "A", Interval(0, 1));
    lg2p.add_edge(1, 3);
    lg2p.add_edge(2, 3);
    lg2p.add_edge(0, 1);
    lg2p.add_edge(0, 2);
    EXPECT_EQ(lg2, lg2p);

    // missing an edge does
    LocalGraph lg2q;
    lg2q.add_node(2, "G", Interval(7, 8));
    lg2q.add_node(3, "T", Interval(13, 14));
    lg2q.add_node(1, "GC", Interval(4, 6));
    lg2q.add_node(0, "A", Interval(0, 1));
    lg2q.add_edge(1, 3);
    lg2q.add_edge(2, 3);
    lg2q.add_edge(0, 1);
    EXPECT_EQ((lg2 == lg2q), false);

    // adding an extra edge does
    lg2p.add_edge(0, 2);
    lg2p.add_edge(0, 3);
    EXPECT_EQ((lg2 == lg2q), false);

    // adding an extra node
    LocalGraph lg2r;
    lg2r.add_node(2, "G", Interval(7, 8));
    lg2r.add_node(3, "T", Interval(13, 14));
    lg2r.add_node(1, "GC", Interval(4, 6));
    lg2r.add_node(0, "A", Interval(0, 1));
    lg2r.add_edge(1, 3);
    lg2r.add_edge(2, 3);
    lg2r.add_edge(0, 1);
    lg2r.add_edge(0, 2);
    lg2r.add_node(4, "T", Interval(15, 16));
    EXPECT_EQ((lg2 == lg2r), false);
}

TEST(LocalGraphTest, not_equals)
{
    LocalGraph lg1;
    lg1.add_node(0, "AGCT", Interval(0, 4));
    EXPECT_EQ((lg1 != lg1), false);

    LocalGraph lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(7, 8));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);
    EXPECT_EQ((lg2 != lg2), false);

    EXPECT_EQ((lg1 != lg2), true);

    // order adding shouldn't matter
    LocalGraph lg2p;
    lg2p.add_node(2, "G", Interval(7, 8));
    lg2p.add_node(3, "T", Interval(13, 14));
    lg2p.add_node(1, "GC", Interval(4, 6));
    lg2p.add_node(0, "A", Interval(0, 1));
    lg2p.add_edge(1, 3);
    lg2p.add_edge(2, 3);
    lg2p.add_edge(0, 1);
    lg2p.add_edge(0, 2);
    EXPECT_EQ((lg2 != lg2p), false);

    // missing an edge does
    LocalGraph lg2q;
    lg2q.add_node(2, "G", Interval(7, 8));
    lg2q.add_node(3, "T", Interval(13, 14));
    lg2q.add_node(1, "GC", Interval(4, 6));
    lg2q.add_node(0, "A", Interval(0, 1));
    lg2q.add_edge(1, 3);
    lg2q.add_edge(2, 3);
    lg2q.add_edge(0, 1);
    EXPECT_EQ((lg2 != lg2q), true);

    // adding an extra edge does
    lg2p.add_edge(0, 2);
    lg2p.add_edge(0, 3);
    EXPECT_EQ((lg2 != lg2q), true);

    // adding an extra node
    LocalGraph lg2r;
    lg2r.add_node(2, "G", Interval(7, 8));
    lg2r.add_node(3, "T", Interval(13, 14));
    lg2r.add_node(1, "GC", Interval(4, 6));
    lg2r.add_node(0, "A", Interval(0, 1));
    lg2r.add_edge(1, 3);
    lg2r.add_edge(2, 3);
    lg2r.add_edge(0, 1);
    lg2r.add_edge(0, 2);
    lg2r.add_node(4, "T", Interval(15, 16));
    EXPECT_EQ((lg2 != lg2r), true);
}

TEST(LocalGraphTest, write_gfa)
{
    LocalGraph lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(7, 8));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);

    lg2.write_gfa("localgraph_test.gfa");
}

TEST(LocalGraphTest, read_gfa)
{
    LocalGraph lg2, read_lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(7, 8));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);

    read_lg2.read_gfa("localgraph_test.gfa");
    EXPECT_EQ(lg2, read_lg2);
}

TEST(LocalGraphTest, walk)
{
    LocalGraph lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(7, 8));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);

    // simple case, there are 2 paths of length 3
    vector<PathPtr> q1;
    prg::Path p;
    deque<Interval> d = { Interval(0, 1), Interval(4, 6) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    d = { Interval(0, 1), Interval(7, 8), Interval(13, 14) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    vector<PathPtr> p1 = lg2.walk(0, 0, 3);
    bool equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // but only one can be extended to a path of length 4
    q1.clear();
    d = { Interval(0, 1), Interval(4, 6), Interval(13, 14) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    p1 = lg2.walk(0, 0, 4);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // for even simpler path of length 1
    q1.clear();
    d = { Interval(0, 1) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    p1 = lg2.walk(0, 0, 1);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // no paths of length 5
    q1.clear();
    p1 = lg2.walk(0, 0, 5);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // 1 path starting from middle var site
    q1.clear();
    d = { Interval(4, 6), Interval(13, 14) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    p1 = lg2.walk(1, 4, 3);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // test on a slightly more complex graph
    LocalGraph lg3;
    lg3.add_node(0, "A", Interval(0, 1));
    lg3.add_node(1, "G", Interval(4, 5));
    lg3.add_node(2, "C", Interval(8, 9));
    lg3.add_node(3, "T", Interval(12, 13));
    lg3.add_node(4, "", Interval(16, 16));
    lg3.add_node(5, "G", Interval(19, 20));
    lg3.add_node(6, "T", Interval(23, 24));
    lg3.add_edge(0, 1);
    lg3.add_edge(0, 5);
    lg3.add_edge(1, 2);
    lg3.add_edge(1, 3);
    lg3.add_edge(2, 4);
    lg3.add_edge(3, 4);
    lg3.add_edge(4, 6);
    lg3.add_edge(5, 6);

    q1.clear();
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9), Interval(16, 16),
        Interval(23, 24) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13), Interval(16, 16),
        Interval(23, 24) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    p1 = lg3.walk(0, 0, 4);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // also want to allow walks starting from an empty node, including the empty node
    q1.clear();
    d = { Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    p1 = lg3.walk(4, 16, 1);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);
}

TEST(LocalGraphTest, walk_back)
{
    LocalGraph lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(7, 8));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);

    // simple case, there are 2 paths of length 3
    vector<PathPtr> q1;
    prg::Path p;
    deque<Interval> d = { Interval(4, 6), Interval(13, 14) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    d = { Interval(0, 1), Interval(7, 8), Interval(13, 14) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    vector<PathPtr> p1 = lg2.walk_back(3, 14, 3);
    bool equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // but only one can be extended to a path of length 4
    q1.clear();
    d = { Interval(0, 1), Interval(4, 6), Interval(13, 14) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    p1 = lg2.walk_back(3, 14, 4);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // for even simpler path of length 1
    q1.clear();
    d = { Interval(0, 1) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    p1 = lg2.walk_back(0, 1, 1);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // no paths of length 5
    q1.clear();
    p1 = lg2.walk_back(3, 14, 5);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // 1 path starting from middle var site
    q1.clear();
    d = { Interval(0, 1), Interval(4, 6) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    p1 = lg2.walk_back(1, 6, 3);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // test on a slightly more complex graph
    LocalGraph lg3;
    lg3.add_node(0, "A", Interval(0, 1));
    lg3.add_node(1, "G", Interval(4, 5));
    lg3.add_node(2, "C", Interval(8, 9));
    lg3.add_node(3, "T", Interval(12, 13));
    lg3.add_node(4, "", Interval(16, 16));
    lg3.add_node(5, "G", Interval(19, 20));
    lg3.add_node(6, "T", Interval(23, 24));
    lg3.add_edge(0, 1);
    lg3.add_edge(0, 5);
    lg3.add_edge(1, 2);
    lg3.add_edge(1, 3);
    lg3.add_edge(2, 4);
    lg3.add_edge(3, 4);
    lg3.add_edge(4, 6);
    lg3.add_edge(5, 6);

    q1.clear();
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9), Interval(16, 16),
        Interval(23, 24) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13), Interval(16, 16),
        Interval(23, 24) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    p1 = lg3.walk_back(6, 24, 4);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    // also want to allow walks starting from an empty node, including the empty node
    q1.clear();
    d = { Interval(8, 9), Interval(16, 16) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    d = { Interval(12, 13), Interval(16, 16) };
    p.initialize(d);
    q1.push_back(std::make_shared<prg::Path>(p));
    p1 = lg3.walk_back(4, 16, 1);
    equal = std::equal(q1.begin(), q1.end(), p1.begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);
}

TEST(LocalGraphTest, nodes_along_string)
{
    LocalGraph lg2, read_lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(7, 8));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);

    vector<LocalNodePtr> v_exp = { lg2.nodes[0], lg2.nodes[1], lg2.nodes[3] };
    vector<LocalNodePtr> v = lg2.nodes_along_string("AGCT");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);

    v_exp = { lg2.nodes[0], lg2.nodes[2], lg2.nodes[3] };
    v = lg2.nodes_along_string("AGT");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);

    v_exp = { lg2.nodes[0], lg2.nodes[1] };
    v = lg2.nodes_along_string("AGC");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);

    v_exp = { lg2.nodes[0], lg2.nodes[1], lg2.nodes[3] };
    v = lg2.nodes_along_string("AGC", true);
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);

    v_exp = { lg2.nodes[0], lg2.nodes[1] };
    v = lg2.nodes_along_string("AgC");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);

    // check for simple prgs
    LocalGraph lg1, read_lg1;
    lg1.add_node(0, "AGTTCGTAGACCAACGCGCT", Interval(0, 20));
    v_exp = { lg1.nodes[0] };
    v = lg1.nodes_along_string("AGTTCGTagACCAACGCGCT");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);
    v_exp = {};
    v = lg1.nodes_along_string("AGTTCGTAGACCAACGCGGT");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);

    // check where have a substring equal to a whole string
    LocalGraph lg3;
    lg3.add_node(0, "A", Interval(0, 1));
    lg3.add_node(1, "GC", Interval(4, 6));
    lg3.add_node(2, "G", Interval(7, 8));
    lg3.add_node(3, "C", Interval(13, 14));
    lg3.add_edge(0, 1);
    lg3.add_edge(0, 2);
    lg3.add_edge(1, 3);
    lg3.add_edge(2, 3);

    v_exp = { lg3.nodes[0], lg3.nodes[1] };
    v = lg3.nodes_along_string("AGC");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);

    v_exp = { lg3.nodes[0], lg3.nodes[2], lg3.nodes[3] };
    v = lg3.nodes_along_string("AGC", true);
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);
}

TEST(LocalGraphTest, top_path)
{
    LocalGraph lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(7, 8));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);

    vector<LocalNodePtr> v_exp = { lg2.nodes[0], lg2.nodes[1], lg2.nodes[3] };
    vector<LocalNodePtr> v = lg2.top_path();
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);

    LocalPRG lp3 = LocalPRG(3, "3", "T 5 G 7 C 8 T 7  6 G 5 TATG");
    v_exp = { lp3.prg.nodes[0], lp3.prg.nodes[1], lp3.prg.nodes[2], lp3.prg.nodes[4],
        lp3.prg.nodes[6] };
    v = lp3.prg.top_path();
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);
}

TEST(LocalGraphTest, bottom_path)
{
    LocalGraph lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(7, 8));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);

    vector<LocalNodePtr> v_exp = { lg2.nodes[0], lg2.nodes[2], lg2.nodes[3] };
    vector<LocalNodePtr> v = lg2.bottom_path();
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);

    LocalPRG lp3 = LocalPRG(3, "3", "T 5 G 7 C 8 T 7  6 G 5 TATG");
    v_exp = { lp3.prg.nodes[0], lp3.prg.nodes[5], lp3.prg.nodes[6] };
    v = lp3.prg.bottom_path();
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, v_exp, v);
}
