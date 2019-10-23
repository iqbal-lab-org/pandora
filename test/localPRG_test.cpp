#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "localPRG.h"
#include "minimizer.h"
#include "minirecord.h"
#include "minihit.h"
#include "interval.h"
#include "prg/path.h"
#include "localgraph.h"
#include "localnode.h"
#include "index.h"
#include "inthash.h"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "utils.h"
#include "seq.h"
#include "kmergraph.h"
#include "kmernode.h"
#include <stdint.h>
#include <iostream>


using namespace std;

TEST(LocalPRGTest, create) {
    LocalPRG l0(0, "empty", "");
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    uint32_t j = 0;
    EXPECT_EQ(j, l0.id);
    EXPECT_EQ("empty", l0.name);
    EXPECT_EQ("", l0.seq);
    j = 1;
    EXPECT_EQ(j, l1.id);
    EXPECT_EQ("simple", l1.name);
    EXPECT_EQ("AGCT", l1.seq);
    j = 2;
    EXPECT_EQ(j, l2.id);
    EXPECT_EQ("varsite", l2.name);
    EXPECT_EQ("A 5 GC 6 G 5 T", l2.seq);
    j = 3;
    EXPECT_EQ(j, l3.id);
    EXPECT_EQ("nested varsite", l3.name);
    EXPECT_EQ("A 5 G 7 C 8 T 7  6 G 5 T", l3.seq);
}

TEST(LocalPRGTest, isalphaEmptyString) {
    LocalPRG l0(0, "empty", "");
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    bool a0 = l0.isalpha_string("");
    bool a1 = l1.isalpha_string("");
    bool a2 = l2.isalpha_string("");
    bool a3 = l3.isalpha_string("");
    EXPECT_EQ(a0, 1) << "isalpha_string thinks the empty string is not alphabetic for l0";
    EXPECT_EQ(a1, 1) << "isalpha_string thinks the empty string is not alphabetic for l1";
    EXPECT_EQ(a2, 1) << "isalpha_string thinks the empty string is not alphabetic for l2";
    EXPECT_EQ(a3, 1) << "isalpha_string thinks the empty string is not alphabetic for l3";
}

TEST(LocalPRGTest, isalphaSpaceString) {
    LocalPRG l0(0, "empty", "");
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    bool a0 = l0.isalpha_string("AGCT T");
    bool a1 = l1.isalpha_string("AGCT T");
    bool a2 = l2.isalpha_string("AGCT T");
    bool a3 = l3.isalpha_string("AGCT T");
    EXPECT_EQ(a0, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l0";
    EXPECT_EQ(a1, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l1";
    EXPECT_EQ(a2, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l2";
    EXPECT_EQ(a3, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l3";
}

TEST(LocalPRGTest, isalphaNumberString) {
    LocalPRG l0(0, "empty", "");
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    bool a0 = l0.isalpha_string("AGCT 8 T");
    bool a1 = l1.isalpha_string("AGCT 8 T");
    bool a2 = l2.isalpha_string("AGCT 8 T");
    bool a3 = l3.isalpha_string("AGCT 8 T");
    EXPECT_EQ(a0, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l0";
    EXPECT_EQ(a1, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l1";
    EXPECT_EQ(a2, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l2";
    EXPECT_EQ(a3, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l3";
}

TEST(LocalPRGTest, string_along_path) {
    LocalPRG l0(0, "empty", "");
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    // empty interval
    deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    EXPECT_EQ("", l0.string_along_path(p));
    EXPECT_EQ("", l1.string_along_path(p));
    EXPECT_EQ("", l2.string_along_path(p));
    EXPECT_EQ("", l3.string_along_path(p));

    // positive length interval
    d = {Interval(1, 3)};
    p.initialize(d);
    EXPECT_EQ("GC", l1.string_along_path(p));
    EXPECT_EQ(" 5", l2.string_along_path(p));
    EXPECT_EQ(" 5", l3.string_along_path(p));

    // multiple intervals
    d = {Interval(0, 1), Interval(2, 3)};
    p.initialize(d);
    EXPECT_EQ("AC", l1.string_along_path(p));
    EXPECT_EQ("A5", l2.string_along_path(p));
    EXPECT_EQ("A5", l3.string_along_path(p));

    // including empty interval
    d = {Interval(0, 1), Interval(2, 2)};
    p.initialize(d);
    EXPECT_EQ("A", l1.string_along_path(p));
    EXPECT_EQ("A", l2.string_along_path(p));
    EXPECT_EQ("A", l3.string_along_path(p));

    // forbidden paths
    d = {Interval(2, 3), Interval(13, 25)};
    p.initialize(d);
    EXPECT_DEATH(l1.string_along_path(p), "");
    EXPECT_DEATH(l1.string_along_path(p), "");
    EXPECT_DEATH(l2.string_along_path(p), "");
    EXPECT_DEATH(l3.string_along_path(p), "");
}

TEST(LocalPRGTest, string_along_localpath) {
    LocalPRG l0(0, "empty", "");
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");

    // empty interval
    vector<LocalNodePtr> p = {l0.prg.nodes[0]};
    EXPECT_EQ("", l0.string_along_path(p));
    p = {l1.prg.nodes[0]};
    EXPECT_EQ("AGCT", l1.string_along_path(p));

    // extract from intervals
    p = {l2.prg.nodes[0], l2.prg.nodes[1]};
    EXPECT_EQ("AGC", l2.string_along_path(p));
    p = {l2.prg.nodes[0], l2.prg.nodes[2], l2.prg.nodes[3]};
    EXPECT_EQ("AGT", l2.string_along_path(p));
}

//core function called to test the LocalPRG::nodes_along_path() method in both tests that follow:
//TEST(LocalPRGTest, nodes_along_path_without_memoization) and
//TEST(LocalPRGTest, nodes_along_path_with_memoization)
void nodes_along_path_core_test() {
    LocalPRG l0(0, "empty", "");
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    // empty interval expects no nodes along
    deque<Interval> d = {Interval(0, 0)};
    prg::Path p1;
    p1.initialize(d);
    prg::Path p2;
    p2.initialize(d);
    prg::Path p3;
    p3.initialize(d);
    vector<LocalNodePtr> v;

    //EXPECT_EQ(v, l0.nodes_along_path(p));
    EXPECT_EQ(v, l1.nodes_along_path(p1));
    EXPECT_EQ(v, l2.nodes_along_path(p2));
    EXPECT_EQ(v, l3.nodes_along_path(p3));

    // positive length interval
    d = {Interval(1, 3)};
    p1.initialize(d);
    p2.initialize(d);
    p3.initialize(d);
    uint32_t j = 1;
    EXPECT_EQ(j, l1.nodes_along_path(p1).size());
    j = 0;
    EXPECT_EQ(j, l1.nodes_along_path(p1)[0]->id);
    EXPECT_EQ(j, l2.nodes_along_path(p2).size()); // no nodes in this interval
    EXPECT_EQ(j, l3.nodes_along_path(p3).size());
    // different interval
    d = {Interval(4, 5)};
    p1.initialize(d);
    p2.initialize(d);
    p3.initialize(d);
    j = 1;
    EXPECT_EQ(j, l2.nodes_along_path(p2).size());
    EXPECT_EQ(j, l3.nodes_along_path(p3).size());
    EXPECT_EQ(j, l2.nodes_along_path(p2)[0]->id);
    EXPECT_EQ(j, l3.nodes_along_path(p3)[0]->id);

    // multiple intervals
    d = {Interval(0, 1), Interval(4, 5)};
    p1.initialize(d);
    p2.initialize(d);
    p3.initialize(d);
    j = 1;
    EXPECT_EQ(j, l1.nodes_along_path(p1).size());
    j = 2;
    EXPECT_EQ(j, l2.nodes_along_path(p2).size());
    EXPECT_EQ(j, l3.nodes_along_path(p3).size());
    j = 0;
    EXPECT_EQ(j, l1.nodes_along_path(p1)[0]->id);
    EXPECT_EQ(j, l2.nodes_along_path(p2)[0]->id);
    EXPECT_EQ(j, l3.nodes_along_path(p3)[0]->id);
    j = 1;
    EXPECT_EQ(j, l2.nodes_along_path(p2)[1]->id);
    EXPECT_EQ(j, l3.nodes_along_path(p3)[1]->id);

    // including empty interval
    d = {Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p1.initialize(d);
    p2.initialize(d);
    p3.initialize(d);
    j = 3;
    vector<LocalNodePtr> w = l3.nodes_along_path(p3);
    EXPECT_EQ(j, w.size());
    EXPECT_EQ(j, l3.nodes_along_path(p3)[0]->id);
    j = 4;
    EXPECT_EQ(j, l3.nodes_along_path(p3)[1]->id);
    j = 6;
    EXPECT_EQ(j, l3.nodes_along_path(p3)[2]->id);

    // a path with an empty node at end
    d = {Interval(12, 13), Interval(16, 16), Interval(23, 23)};
    p1.initialize(d);
    p2.initialize(d);
    p3.initialize(d);
    j = 3;
    w = l3.nodes_along_path(p3);
    EXPECT_EQ(j, w.size());
    EXPECT_EQ(j, l3.nodes_along_path(p3)[0]->id);
    j = 4;
    EXPECT_EQ(j, l3.nodes_along_path(p3)[1]->id);
    j = 6;
    EXPECT_EQ(j, l3.nodes_along_path(p3)[2]->id);

    // and a path which ends on a null node
    d = {Interval(12, 13), Interval(16, 16)};
    p1.initialize(d);
    p2.initialize(d);
    p3.initialize(d);
    j = 2;
    w = l3.nodes_along_path(p3);
    EXPECT_EQ(j, w.size());
    j = 3;
    EXPECT_EQ(j, l3.nodes_along_path(p3)[0]->id);
    j = 4;
    EXPECT_EQ(j, l3.nodes_along_path(p3)[1]->id);

    // and a path that can't really exist still works
    d = {Interval(12, 13), Interval(19, 20)};
    p1.initialize(d);
    p2.initialize(d);
    p3.initialize(d);
    j = 2;
    EXPECT_EQ(j, l3.nodes_along_path(p3).size());
    j = 3;
    EXPECT_EQ(j, l3.nodes_along_path(p3)[0]->id);
    j = 5;
    EXPECT_EQ(j, l3.nodes_along_path(p3)[1]->id);
}

TEST(LocalPRGTest, nodes_along_path_without_memoization) {
    LocalPRG::do_path_memoization_in_nodes_along_path_method = false;
    nodes_along_path_core_test();
}

TEST(LocalPRGTest, nodes_along_path_with_memoization) {
    LocalPRG::do_path_memoization_in_nodes_along_path_method = true;
    nodes_along_path_core_test();
    LocalPRG::do_path_memoization_in_nodes_along_path_method = false;
}

TEST(LocalPRGTest, split_by_siteNoSites) {
    LocalPRG l0(0, "empty", "");
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    vector<Interval> v0, v1;
    v0.push_back(Interval(0, 0));
    EXPECT_ITERABLE_EQ(vector<Interval>, v0,
                       l0.split_by_site(Interval(0, 0)));// << "Failed to split empty string with input Interval";
    v1.push_back(Interval(0, 4));
    EXPECT_ITERABLE_EQ(vector<Interval>, v1,
                       l1.split_by_site(Interval(0, 4)));// << "Failed to split string with input Interval";
    v1.clear();
    v1.push_back(Interval(0, 2));
    EXPECT_ITERABLE_EQ(vector<Interval>, v1,
                       l1.split_by_site(Interval(0, 2)));// << "Failed to split string with short input Interval";
    v1.clear();
    v1.push_back(Interval(1, 3));
    EXPECT_ITERABLE_EQ(vector<Interval>, v1,
                       l1.split_by_site(Interval(1, 3)));// << "Failed to split string with middle input Interval";
}

TEST(LocalPRGTest, split_by_siteSite) {
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");

    vector<Interval> v2;
    v2.push_back(Interval(0, 1));
    l2.next_site = 5;
    EXPECT_ITERABLE_EQ(vector<Interval>, v2,
                       l2.split_by_site(Interval(0, 1)));// << "Failed to split string in Interval" << Interval(0,1);
    EXPECT_ITERABLE_EQ(vector<Interval>, v2,
                       l2.split_by_site(Interval(0, 2)));// << "Failed to split string in Interval" << Interval(0,2);
    EXPECT_ITERABLE_EQ(vector<Interval>, v2,
                       l2.split_by_site(Interval(0, 3)));// << "Failed to split string in Interval" << Interval(0,3);
    //EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,4)));// << "Failed to split string in Interval" << Interval(0,4);
    v2.push_back(Interval(4, 6));
    EXPECT_ITERABLE_EQ(vector<Interval>, v2,
                       l2.split_by_site(Interval(0, 6)));// << "Failed to split string in Interval" << Interval(0,6);
    EXPECT_ITERABLE_EQ(vector<Interval>, v2,
                       l2.split_by_site(Interval(0, 7)));// << "Failed to split string in Interval" << Interval(0,7);
    EXPECT_ITERABLE_EQ(vector<Interval>, v2,
                       l2.split_by_site(Interval(0, 8)));// << "Failed to split string in Interval" << Interval(0,8);
    v2.push_back(Interval(9, 10));
    EXPECT_ITERABLE_EQ(vector<Interval>, v2,
                       l2.split_by_site(Interval(0, 10)));// << "Failed to split string in Interval" << Interval(0,10);
    EXPECT_ITERABLE_EQ(vector<Interval>, v2,
                       l2.split_by_site(Interval(0, 11)));// << "Failed to split string in Interval" << Interval(0,11);
    EXPECT_ITERABLE_EQ(vector<Interval>, v2,
                       l2.split_by_site(Interval(0, 12)));// << "Failed to split string in Interval" << Interval(0,12);
    v2.push_back(Interval(13, 14));
    EXPECT_ITERABLE_EQ(vector<Interval>, v2,
                       l2.split_by_site(Interval(0, 14)));// << "Failed to split string in Interval" << Interval(0,14);
    v2.clear();
    v2.push_back(Interval(5, 6));
    EXPECT_ITERABLE_EQ(vector<Interval>, v2, l2.split_by_site(
            Interval(5, 8)));// << "Failed to split string in mid Interval" << Interval(5,8);
}

TEST(LocalPRGTest, split_by_siteNestedSite) {
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");
    LocalPRG l4(4, "nested varsite start immediately", " 5 G 7 C 8 T 7  6 G 5 ");

    vector<Interval> v3;
    v3.push_back(Interval(0, 1));
    l3.next_site = 5;
    EXPECT_ITERABLE_EQ(vector<Interval>, v3,
                       l3.split_by_site(Interval(0, 1)));// << "Failed to split string in Interval" << Interval(0,1);
    EXPECT_ITERABLE_EQ(vector<Interval>, v3,
                       l3.split_by_site(Interval(0, 2)));// << "Failed to split string in Interval" << Interval(0,2);
    EXPECT_ITERABLE_EQ(vector<Interval>, v3,
                       l3.split_by_site(Interval(0, 3)));// << "Failed to split string in Interval" << Interval(0,3);
    //EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,4)));// << "Failed to split string in Interval" << Interval(0,4);
    v3.push_back(Interval(4, 16));
    EXPECT_ITERABLE_EQ(vector<Interval>, v3,
                       l3.split_by_site(Interval(0, 16)));// << "Failed to split string in Interval" << Interval(0,6);
    EXPECT_ITERABLE_EQ(vector<Interval>, v3,
                       l3.split_by_site(Interval(0, 17)));// << "Failed to split string in Interval" << Interval(0,7);
    EXPECT_ITERABLE_EQ(vector<Interval>, v3,
                       l3.split_by_site(Interval(0, 18)));// << "Failed to split string in Interval" << Interval(0,8);
    v3.push_back(Interval(19, 20));
    EXPECT_ITERABLE_EQ(vector<Interval>, v3,
                       l3.split_by_site(Interval(0, 20)));// << "Failed to split string in Interval" << Interval(0,10);
    EXPECT_ITERABLE_EQ(vector<Interval>, v3,
                       l3.split_by_site(Interval(0, 21)));// << "Failed to split string in Interval" << Interval(0,11);
    EXPECT_ITERABLE_EQ(vector<Interval>, v3,
                       l3.split_by_site(Interval(0, 22)));// << "Failed to split string in Interval" << Interval(0,12);
    v3.push_back(Interval(23, 24));
    EXPECT_ITERABLE_EQ(vector<Interval>, v3,
                       l3.split_by_site(Interval(0, 24)));// << "Failed to split string in Interval" << Interval(0,14);
    l3.next_site = 7;
    v3.clear();
    v3.push_back(Interval(4, 5));
    v3.push_back(Interval(8, 9));
    v3.push_back(Interval(12, 13));
    v3.push_back(Interval(16, 16));
    EXPECT_ITERABLE_EQ(vector<Interval>, v3, l3.split_by_site(
            Interval(4, 16)));// << "Failed to split string in mid Interval" << Interval(5,8);

    vector<Interval> v4;
    v4.push_back(Interval(0, 0));
    l4.next_site = 5;
    EXPECT_ITERABLE_EQ(vector<Interval>, v4,
                       l4.split_by_site(Interval(0, 1)));// << "Failed to split string in Interval" << Interval(0,1);
    EXPECT_ITERABLE_EQ(vector<Interval>, v4,
                       l4.split_by_site(Interval(0, 2)));// << "Failed to split string in Interval" << Interval(0,2);
    v4.push_back(Interval(3, 15));
    EXPECT_ITERABLE_EQ(vector<Interval>, v4,
                       l4.split_by_site(Interval(0, 15)));// << "Failed to split string in Interval" << Interval(0,6);
    EXPECT_ITERABLE_EQ(vector<Interval>, v4,
                       l4.split_by_site(Interval(0, 16)));// << "Failed to split string in Interval" << Interval(0,7);
    EXPECT_ITERABLE_EQ(vector<Interval>, v4,
                       l4.split_by_site(Interval(0, 17)));// << "Failed to split string in Interval" << Interval(0,8);
    v4.push_back(Interval(18, 19));
    EXPECT_ITERABLE_EQ(vector<Interval>, v4,
                       l4.split_by_site(Interval(0, 19)));// << "Failed to split string in Interval" << Interval(0,10);
    EXPECT_ITERABLE_EQ(vector<Interval>, v4,
                       l4.split_by_site(Interval(0, 20)));// << "Failed to split string in Interval" << Interval(0,11);
    l4.next_site = 7;
    v4.clear();
    v4.push_back(Interval(0, 4));
    v4.push_back(Interval(7, 8));
    v4.push_back(Interval(11, 12));
    v4.push_back(Interval(15, 22));
    EXPECT_ITERABLE_EQ(vector<Interval>, v4, l4.split_by_site(
            Interval(0, 22)));// << "Failed to split string in mid Interval" << Interval(5,8);

}

TEST(LocalPRGTest, build_graph) {
    LocalPRG l0(0, "empty", "");
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    LocalGraph lg0;
    lg0.add_node(0, "", Interval(0, 0));
    EXPECT_EQ(lg0, l0.prg);

    LocalGraph lg1;
    lg1.add_node(0, "AGCT", Interval(0, 4));
    EXPECT_EQ(lg1, l1.prg);

    LocalGraph lg2;
    lg2.add_node(0, "A", Interval(0, 1));
    lg2.add_node(1, "GC", Interval(4, 6));
    lg2.add_node(2, "G", Interval(9, 10));
    lg2.add_node(3, "T", Interval(13, 14));
    lg2.add_edge(0, 1);
    lg2.add_edge(0, 2);
    lg2.add_edge(1, 3);
    lg2.add_edge(2, 3);
    EXPECT_EQ(lg2, l2.prg);

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
    EXPECT_EQ(lg3, l3.prg);
}

TEST(LocalPRGTest, shift) {
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "AT 5 G 7 C 8 T 7  6 G 5 T");
    //LocalPRG l4(4, "much more complex", "TCATTC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AGCTG");
    LocalPRG l5(5, "one with lots of null at start and end, and a long stretch in between",
                " 5  7  9  11 AGTTCTGAAACATTGCGCGTGAGATCTCTG 12 T 11  10 A 9  8 C 7  6 G 5 ");
    LocalPRG l6(6, "one representing a possible deletion at end", "GATCTCTAG 5 TTATG 6  5 ");

    deque<Interval> d = {Interval(0, 3)};
    prg::Path p, q;
    p.initialize(d);
    d = {Interval(1, 4)};
    q.initialize(d);
    vector<PathPtr> v_exp = { std::make_shared<prg::Path>(q) };
    bool equal = std::equal(v_exp.begin(), v_exp.end(), l1.shift(p).begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);
    v_exp.clear();
    equal = std::equal(v_exp.begin(), v_exp.end(), l1.shift(q).begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    d = {Interval(0, 1), Interval(4, 6)};
    p.initialize(d);
    d = {Interval(4, 6), Interval(13, 14)};
    q.initialize(d);
    v_exp = { std::make_shared<prg::Path>(q) };
    equal = std::equal(v_exp.begin(), v_exp.end(), l2.shift(p).begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);
    v_exp.clear();
    equal = std::equal(v_exp.begin(), v_exp.end(), l2.shift(q).begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    v_exp.clear();
    d = {Interval(0, 2)};
    p.initialize(d);
    d = {Interval(1, 2), Interval(5, 6)};
    q.initialize(d);
    v_exp.push_back(std::make_shared<prg::Path>(q));
    d = {Interval(1, 2), Interval(20, 21)};
    q.initialize(d);
    v_exp.push_back(std::make_shared<prg::Path>(q));
    equal = std::equal(v_exp.begin(), v_exp.end(), l3.shift(p).begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    v_exp.clear();
    d = {Interval(1, 2), Interval(5, 6)};
    p.initialize(d);
    d = {Interval(5, 6), Interval(9, 10)};
    q.initialize(d);
    v_exp.push_back(std::make_shared<prg::Path>(q));
    d = {Interval(5, 6), Interval(13, 14)};
    q.initialize(d);
    v_exp.push_back(std::make_shared<prg::Path>(q));
    equal = std::equal(v_exp.begin(), v_exp.end(), l3.shift(p).begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    v_exp.clear();
    d = {Interval(0, 0), Interval(3, 3), Interval(6, 6), Interval(9, 9), Interval(13, 18)};
    p.initialize(d);
    d = {Interval(14, 19)};
    q.initialize(d);
    v_exp.push_back(std::make_shared<prg::Path>(q));
    equal = std::equal(v_exp.begin(), v_exp.end(), l5.shift(p).begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    v_exp.clear();
    d = {Interval(3, 8)};
    p.initialize(d);
    d = {Interval(4, 9), Interval(20, 20), Interval(23, 23)};
    q.initialize(d);
    v_exp.push_back(std::make_shared<prg::Path>(q));
    d = {Interval(4, 9)};
    q.initialize(d);
    v_exp.push_back(std::make_shared<prg::Path>(q));
    equal = std::equal(v_exp.begin(), v_exp.end(), l6.shift(p).begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);

    v_exp.clear();
    d = {Interval(4, 9)};
    p.initialize(d);
    d = {Interval(5, 9), Interval(12, 13)};
    q.initialize(d);
    v_exp.push_back(std::make_shared<prg::Path>(q));
    equal = std::equal(v_exp.begin(), v_exp.end(), l6.shift(p).begin(), ComparePathPtr());
    EXPECT_EQ(equal, true);
}

TEST(LocalPRGTest, minimizer_sketch) {
    // note this is a bad test
    LocalPRG l0(0, "empty", "");
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");
    LocalPRG l4(4, "much more complex", "TCATTC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AGCTG");
    LocalPRG l5(5, "one with lots of null at start and end, and a long stretch in between",
                " 5  7  9  11 AGTTCTGAAACATTGCGCGTGAGATCTCTG 12 T 11  10 A 9  8 C 7  6 G 5 ");
    LocalPRG l6(2, "too short for w and k", "A 5 GC 6 G 5 T");

    auto index = std::make_shared<Index>();
    KmerHash hash;

    l0.minimizer_sketch(index, 1, 3);
    uint32_t j = 0;
    EXPECT_EQ(j, index->minhash.size());

    l1.minimizer_sketch(index, 2, 3);
    j = 1;
    EXPECT_EQ(j, index->minhash.size());
    l1.minimizer_sketch(index, 1, 3);
    EXPECT_EQ(j, index->minhash.size());
    j = 2;
    pair<uint64_t, uint64_t> kh = hash.kmerhash("AGC", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());

    index->clear();
    l2.minimizer_sketch(index, 2, 3);
    j = 1;
    EXPECT_EQ(j, index->minhash.size());
    l2.minimizer_sketch(index, 1, 3);
    j = 2;
    EXPECT_EQ(j, index->minhash.size());
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    j = 1;
    kh = hash.kmerhash("AGT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());

    index->clear();
    l3.minimizer_sketch(index, 2, 3);
    j = 2;
    EXPECT_EQ(j, index->minhash.size());
    l3.minimizer_sketch(index, 1, 3);
    j = 3;
    EXPECT_EQ(j, index->minhash.size());
    j = 2;
    kh = hash.kmerhash("AGC", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size()); //AGC
    kh = hash.kmerhash("AGT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size()); //AGTx2
    j = 1;
    kh = hash.kmerhash("GTT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());

    index->clear();
    l4.minimizer_sketch(index, 1, 3);
    j = 16;
    EXPECT_EQ(j, index->minhash.size());
    j = 5;
    kh = hash.kmerhash("TCA", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    j = 4;
    kh = hash.kmerhash("CTA", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    j = 3;
    kh = hash.kmerhash("ACT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("CAA", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("AAG", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("TCT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("AGC", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    j = 2;
    kh = hash.kmerhash("TTC", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("CAC", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("CTC", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    j = 1;
    kh = hash.kmerhash("CAT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("ATT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("GTC", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("GTT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("TGT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("CTG", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());

    index->clear();
    l4.minimizer_sketch(index, 3, 3);
    j = 10;
    EXPECT_EQ(j, index->minhash.size());
    j = 4;
    kh = hash.kmerhash("CTA", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    j = 3;
    kh = hash.kmerhash("CTT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    j = 2;
    kh = hash.kmerhash("CAC", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    j = 1;
    kh = hash.kmerhash("ATT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("ACT", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("TCA", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("AAC", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("GTC", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("GAG", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
    kh = hash.kmerhash("CTG", 3);
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());

    index->clear();
    l5.minimizer_sketch(index, 4, 5);
    EXPECT_EQ((index->minhash.size() > 2), true);

    index->clear();
    l6.minimizer_sketch(index, 2, 4);
    EXPECT_EQ((uint)0, index->minhash.size());

    index->clear();
}

struct MiniPos {
    bool operator()(Minimizer lhs, Minimizer rhs) {
        return (lhs.pos_of_kmer_in_read.start) < (rhs.pos_of_kmer_in_read.start);
    }
};

TEST(LocalPRGTest, minimizer_sketch_SameAsSeqw1) {
    std::string st = "ATGGCAATCCGAATCTTCGCGATACTTTTCTCCATTTTTTCTCTTGCCACTTTCGCGCATGCGCAAGAAGGCACGCTAGAACGTTCTGACTGGAGGAAGTTTTTCAGCGAATTTCAAGCCAAAGGCACGATAGTTGTGGCAGACGAACGCCAAGCGGATCGTGCCATGTTGGTTTTTGATCCTGTGCGATCGAAGAAACGCTACTCGCCTGCATCGACATTCAAGATACCTCATACACTTTTTGCACTTGATGCAGGCGCTGTTCGTGATGAGTTCCAGATTTTTCGATGGGACGGCGTTAACAGGGGCTTTGCAGGCCACAATCAAGACCAAGATTTGCGATCAGCAATGCGGAATTCTACTGTTTGGGTGTATGAGCTATTTGCAAAGGAAATTGGTGATGACAAAGCTCGGCGCTATTTGAAGAAAATCGACTATGGCAACGCCGATCCTTCGACAAGTAATGGCGATTACTGTATAGAAGGCAGCCTTGCAATCTCGGCGCAGGAGCAAATTGCATTTCTCAGGAAGCTCTATCGTAACGAGCTGCCCTTTCGGGTAGAACATCAGCGCTTGGTCAAGGATCTCATGATTGTGGAAGCCGGTCGCAACTGGATACTGCGTGCAAAGACGGGCTGGGAAGGCCGTATGGGTTGGTGGGTAGGATGGGTTGAGTGGCCGACTGGCTCCGTATTCTTCGCACTGAATATTGATACGCCAAACAGAATGGATGATCTTTTCAAGAGGGAGGCAATCGTGCGGGCAATCCTT";

    auto index = std::make_shared<Index>();
    LocalPRG l(0, "prg", st);
    l.minimizer_sketch(index, 1, 15);

    Seq s = Seq(0, "read", st, 1, 15);

    EXPECT_EQ(l.kmer_prg.nodes.size(), s.sketch.size() + 2);

    std::set<Minimizer, MiniPos> sketch(s.sketch.begin(), s.sketch.end());
    auto lit = l.kmer_prg.sorted_nodes.begin();
    lit++;

    for (auto sit = sketch.begin(); sit != sketch.end(); ++sit) {
        EXPECT_EQ((*sit).pos_of_kmer_in_read, (*lit)->path[0]);
        ++lit;
    }
}

TEST(LocalPRGTest, minimizer_sketch_SameAsSeqw5) {
    string st = "ATGGCAATCCGAATCTTCGCGATACTTTTCTCCATTTTTTCTCTTGCCACTTTCGCGCATGCGCAAGAAGGCACGCTAGAACGTTCTGACTGGAGGAAGTTTTTCAGCGAATTTCAAGCCAAAGGCACGATAGTTGTGGCAGACGAACGCCAAGCGGATCGTGCCATGTTGGTTTTTGATCCTGTGCGATCGAAGAAACGCTACTCGCCTGCATCGACATTCAAGATACCTCATACACTTTTTGCACTTGATGCAGGCGCTGTTCGTGATGAGTTCCAGATTTTTCGATGGGACGGCGTTAACAGGGGCTTTGCAGGCCACAATCAAGACCAAGATTTGCGATCAGCAATGCGGAATTCTACTGTTTGGGTGTATGAGCTATTTGCAAAGGAAATTGGTGATGACAAAGCTCGGCGCTATTTGAAGAAAATCGACTATGGCAACGCCGATCCTTCGACAAGTAATGGCGATTACTGTATAGAAGGCAGCCTTGCAATCTCGGCGCAGGAGCAAATTGCATTTCTCAGGAAGCTCTATCGTAACGAGCTGCCCTTTCGGGTAGAACATCAGCGCTTGGTCAAGGATCTCATGATTGTGGAAGCCGGTCGCAACTGGATACTGCGTGCAAAGACGGGCTGGGAAGGCCGTATGGGTTGGTGGGTAGGATGGGTTGAGTGGCCGACTGGCTCCGTATTCTTCGCACTGAATATTGATACGCCAAACAGAATGGATGATCTTTTCAAGAGGGAGGCAATCGTGCGGGCAATCCTT";

    auto index = std::make_shared<Index>();
    LocalPRG l(0, "prg", st);
    l.minimizer_sketch(index, 5, 15);

    Seq s = Seq(0, "read", st, 5, 15);

    EXPECT_EQ(l.kmer_prg.nodes.size(), s.sketch.size() + 2);

    set<Minimizer, MiniPos> sketch(s.sketch.begin(), s.sketch.end());
    auto lit = l.kmer_prg.sorted_nodes.begin();
    lit++;

    for (auto sit = sketch.begin(); sit != sketch.end(); ++sit) {
        EXPECT_EQ((*sit).pos_of_kmer_in_read, (*lit)->path[0]);
        ++lit;
    }
}

TEST(LocalPRGTest, minimizer_sketch_SameAsSeqw10) {
    string st = "ATGGCAATCCGAATCTTCGCGATACTTTTCTCCATTTTTTCTCTTGCCACTTTCGCGCATGCGCAAGAAGGCACGCTAGAACGTTCTGACTGGAGGAAGTTTTTCAGCGAATTTCAAGCCAAAGGCACGATAGTTGTGGCAGACGAACGCCAAGCGGATCGTGCCATGTTGGTTTTTGATCCTGTGCGATCGAAGAAACGCTACTCGCCTGCATCGACATTCAAGATACCTCATACACTTTTTGCACTTGATGCAGGCGCTGTTCGTGATGAGTTCCAGATTTTTCGATGGGACGGCGTTAACAGGGGCTTTGCAGGCCACAATCAAGACCAAGATTTGCGATCAGCAATGCGGAATTCTACTGTTTGGGTGTATGAGCTATTTGCAAAGGAAATTGGTGATGACAAAGCTCGGCGCTATTTGAAGAAAATCGACTATGGCAACGCCGATCCTTCGACAAGTAATGGCGATTACTGTATAGAAGGCAGCCTTGCAATCTCGGCGCAGGAGCAAATTGCATTTCTCAGGAAGCTCTATCGTAACGAGCTGCCCTTTCGGGTAGAACATCAGCGCTTGGTCAAGGATCTCATGATTGTGGAAGCCGGTCGCAACTGGATACTGCGTGCAAAGACGGGCTGGGAAGGCCGTATGGGTTGGTGGGTAGGATGGGTTGAGTGGCCGACTGGCTCCGTATTCTTCGCACTGAATATTGATACGCCAAACAGAATGGATGATCTTTTCAAGAGGGAGGCAATCGTGCGGGCAATCCTT";

    auto index = std::make_shared<Index>();
    LocalPRG l(0, "prg", st);
    l.minimizer_sketch(index, 10, 15);

    Seq s = Seq(0, "read", st, 10, 15);

    EXPECT_EQ(l.kmer_prg.nodes.size(), s.sketch.size() + 2);

    set<Minimizer, MiniPos> sketch(s.sketch.begin(), s.sketch.end());
    auto lit = l.kmer_prg.sorted_nodes.begin();
    lit++;

    for (auto sit = sketch.begin(); sit != sketch.end(); ++sit) {
        EXPECT_EQ((*sit).pos_of_kmer_in_read, (*lit)->path[0]);
        ++lit;
    }
}

TEST(LocalPRGTest, minimizer_sketch_SameAsSeqw15) {
    string st = "ATGGCAATCCGAATCTTCGCGATACTTTTCTCCATTTTTTCTCTTGCCACTTTCGCGCATGCGCAAGAAGGCACGCTAGAACGTTCTGACTGGAGGAAGTTTTTCAGCGAATTTCAAGCCAAAGGCACGATAGTTGTGGCAGACGAACGCCAAGCGGATCGTGCCATGTTGGTTTTTGATCCTGTGCGATCGAAGAAACGCTACTCGCCTGCATCGACATTCAAGATACCTCATACACTTTTTGCACTTGATGCAGGCGCTGTTCGTGATGAGTTCCAGATTTTTCGATGGGACGGCGTTAACAGGGGCTTTGCAGGCCACAATCAAGACCAAGATTTGCGATCAGCAATGCGGAATTCTACTGTTTGGGTGTATGAGCTATTTGCAAAGGAAATTGGTGATGACAAAGCTCGGCGCTATTTGAAGAAAATCGACTATGGCAACGCCGATCCTTCGACAAGTAATGGCGATTACTGTATAGAAGGCAGCCTTGCAATCTCGGCGCAGGAGCAAATTGCATTTCTCAGGAAGCTCTATCGTAACGAGCTGCCCTTTCGGGTAGAACATCAGCGCTTGGTCAAGGATCTCATGATTGTGGAAGCCGGTCGCAACTGGATACTGCGTGCAAAGACGGGCTGGGAAGGCCGTATGGGTTGGTGGGTAGGATGGGTTGAGTGGCCGACTGGCTCCGTATTCTTCGCACTGAATATTGATACGCCAAACAGAATGGATGATCTTTTCAAGAGGGAGGCAATCGTGCGGGCAATCCTT";

    auto index = std::make_shared<Index>();
    LocalPRG l(0, "prg", st);
    l.minimizer_sketch(index, 15, 15);

    Seq s = Seq(0, "read", st, 15, 15);

    EXPECT_EQ(l.kmer_prg.nodes.size(), s.sketch.size() + 2);

    set<Minimizer, MiniPos> sketch(s.sketch.begin(), s.sketch.end());
    auto lit = l.kmer_prg.sorted_nodes.begin();
    lit++;

    for (auto sit = sketch.begin(); sit != sketch.end(); ++sit) {
        EXPECT_EQ((*sit).pos_of_kmer_in_read, (*lit)->path[0]);
        ++lit;
    }
}

TEST(LocalPRGTest, localnode_path_from_kmernode_path) {
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");
    LocalPRG l4(4, "much more complex", "TC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AG");

    auto index = std::make_shared<Index>();

    KmerHash hash;

    l3.minimizer_sketch(index, 2, 3);
    //vector<KmerNodePtr> kmp = {l3.kmer_prg.nodes[0], l3.kmer_prg.nodes[1], l3.kmer_prg.nodes[2], l3.kmer_prg.nodes[4]};
    vector<KmerNodePtr> kmp = {l3.kmer_prg.nodes[2], l3.kmer_prg.nodes[4]};
    vector<LocalNodePtr> lmp = l3.localnode_path_from_kmernode_path(kmp);
    vector<LocalNodePtr> lmp_exp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                    l3.prg.nodes[6]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, lmp_exp, lmp);
    lmp = l3.localnode_path_from_kmernode_path(kmp, 2);
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, lmp_exp, lmp);

    index->clear();
    l4.minimizer_sketch(index, 3, 3);
    //kmp = {l4.kmer_prg.nodes[0], l4.kmer_prg.nodes[1], l4.kmer_prg.nodes[3], l4.kmer_prg.nodes[7], l4.kmer_prg.nodes[9], l4.kmer_prg.nodes[11], l4.kmer_prg.nodes[13]};
    kmp = {l4.kmer_prg.nodes[3], l4.kmer_prg.nodes[7]};
    lmp = l4.localnode_path_from_kmernode_path(kmp, 2);
    lmp_exp = {l4.prg.nodes[0], l4.prg.nodes[1], l4.prg.nodes[3], l4.prg.nodes[4], l4.prg.nodes[6]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, lmp_exp, lmp);
    lmp = l4.localnode_path_from_kmernode_path(kmp, 3);
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, lmp_exp, lmp);
}

TEST(LocalPRGTest, kmernode_path_from_localnode_path) {
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");
    LocalPRG l4(4, "much more complex", "TC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AG");
    LocalPRG l5(5, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");

    auto index = std::make_shared<Index>();
    KmerHash hash;

    l3.minimizer_sketch(index, 2, 3);
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6]};

    vector<KmerNodePtr> kmp = l3.kmernode_path_from_localnode_path(lmp);
    sort(kmp.begin(), kmp.end());

    vector<KmerNodePtr> kmp_exp = {l3.kmer_prg.nodes[0], l3.kmer_prg.nodes[1], l3.kmer_prg.nodes[2],
                                   l3.kmer_prg.nodes[4]};
    sort(kmp_exp.begin(), kmp_exp.end());
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, kmp_exp, kmp);

    index->clear();
    l4.minimizer_sketch(index, 3, 3);
    lmp = {l4.prg.nodes[0], l4.prg.nodes[1], l4.prg.nodes[3], l4.prg.nodes[4], l4.prg.nodes[6]};

    kmp = l4.kmernode_path_from_localnode_path(lmp);
    sort(kmp.begin(), kmp.end());

    kmp_exp = {l4.kmer_prg.nodes[0], l4.kmer_prg.nodes[1], l4.kmer_prg.nodes[3], l4.kmer_prg.nodes[7],
               l4.kmer_prg.nodes[9], l4.kmer_prg.nodes[11], l4.kmer_prg.nodes[13]};
    sort(kmp_exp.begin(), kmp_exp.end());
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, kmp_exp, kmp);

    // case where we don't have start and end point in localpath, so need to consider whether kmer overlaps
    index->clear();
    l5.minimizer_sketch(index, 2, 3);
    lmp = {l5.prg.nodes[1], l5.prg.nodes[2], l5.prg.nodes[4], l5.prg.nodes[6], l5.prg.nodes[7]};

    kmp = l5.kmernode_path_from_localnode_path(lmp);
    sort(kmp.begin(), kmp.end());

    kmp_exp = {l5.kmer_prg.nodes[1], l5.kmer_prg.nodes[2], l5.kmer_prg.nodes[6], l5.kmer_prg.nodes[8],
               l5.kmer_prg.nodes[10], l5.kmer_prg.nodes[12], l5.kmer_prg.nodes[13]};
    sort(kmp_exp.begin(), kmp_exp.end());

    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, kmp_exp, kmp);
}


TEST(GetCovgsAlongLocalnodePathTest, emptyPanNodeReturnsEmpty) {
    const auto prg_id { 3 };
    auto local_prg_ptr { std::make_shared<LocalPRG>(prg_id, "test", "") };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path;
    const std::vector<KmerNodePtr> kmer_node_max_likelihood_path;

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(local_prg_ptr) };

    const auto actual {
            get_covgs_along_localnode_path(pangraph_node, local_node_max_likelihood_path, kmer_node_max_likelihood_path,
                                           0) };
    const std::vector<uint32_t> expected;

    EXPECT_EQ(actual, expected);
}

TEST(GetCovgsAlongLocalnodePathTest, get_covgs_along_localnode_path) {
    auto l3 { std::make_shared<LocalPRG>(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T") };
    auto l4 { std::make_shared<LocalPRG>(4, "much more complex", "TC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AG") };

    auto index = std::make_shared<Index>();
    KmerHash hash;

    l3->minimizer_sketch(index, 2, 3);
    vector<KmerNodePtr> kmp = {l3->kmer_prg.nodes[2], l3->kmer_prg.nodes[4]};
    vector<LocalNodePtr> lmp = l3->localnode_path_from_kmernode_path(kmp, 2);
    shared_ptr<pangenome::Node> pn3(make_shared<pangenome::Node>(l3));
    for (const auto &n : pn3->kmer_prg_with_coverage.kmer_prg->nodes) {
        pn3->kmer_prg_with_coverage.increment_covg(n->id, 0, 0);
    }
    vector<uint> covgs = get_covgs_along_localnode_path(pn3, lmp, kmp, 0);
    vector<uint> covgs_exp = {0, 1, 1, 1};
    EXPECT_ITERABLE_EQ(vector<uint>, covgs_exp, covgs);

    index->clear();
    l4->minimizer_sketch(index, 1, 3);
    kmp = {l4->kmer_prg.nodes[0], l4->kmer_prg.nodes[1], l4->kmer_prg.nodes[3], l4->kmer_prg.nodes[5], l4->kmer_prg.nodes[7],
           l4->kmer_prg.nodes[9], l4->kmer_prg.nodes[12], l4->kmer_prg.nodes[15], l4->kmer_prg.nodes[18],
           l4->kmer_prg.nodes[21], l4->kmer_prg.nodes[23], l4->kmer_prg.nodes[25], l4->kmer_prg.nodes[27],
           l4->kmer_prg.nodes[29]};
    lmp = l4->localnode_path_from_kmernode_path(kmp, 1);
    shared_ptr<pangenome::Node> pn4(make_shared<pangenome::Node>(l4));
    for (const auto &n : pn4->kmer_prg_with_coverage.kmer_prg->nodes) {
        pn4->kmer_prg_with_coverage.increment_covg(n->id, 0, 0);
    }
    covgs = get_covgs_along_localnode_path(pn4, lmp, kmp, 0);
    //covgs_exp = {1,2,3,3,3,3,3,3,3,3,3,3,2,1};
    covgs_exp = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    EXPECT_ITERABLE_EQ(vector<uint>, covgs_exp, covgs);

    kmp = {l4->kmer_prg.nodes[0], l4->kmer_prg.nodes[3], l4->kmer_prg.nodes[5], l4->kmer_prg.nodes[12],
           l4->kmer_prg.nodes[15], l4->kmer_prg.nodes[18], l4->kmer_prg.nodes[25]};
    lmp = l4->localnode_path_from_kmernode_path(kmp, 2);
    covgs = get_covgs_along_localnode_path(pn4, lmp, kmp, 0);
    //covgs_exp = {0,1,2,2,1,1,2,3,2,1,1,1,1,0};
    covgs_exp = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0};

    EXPECT_ITERABLE_EQ(vector<uint>, covgs_exp, covgs);
}


TEST(LocalPRGTest, write_covgs_to_file) {
    auto l3 { std::make_shared<LocalPRG>(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T") };

    auto index = std::make_shared<Index>();
    KmerHash hash;

    l3->minimizer_sketch(index, 2, 3);
    vector<KmerNodePtr> kmp = {l3->kmer_prg.nodes[2], l3->kmer_prg.nodes[4]};
    vector<LocalNodePtr> lmp = l3->localnode_path_from_kmernode_path(kmp, 2);
    shared_ptr<pangenome::Node> pn3(make_shared<pangenome::Node>(l3));
    for (const auto &n : pn3->kmer_prg_with_coverage.kmer_prg->nodes) {
        pn3->kmer_prg_with_coverage.increment_covg(n->id, 0, 0);
    }
    vector<uint> covgs = get_covgs_along_localnode_path(pn3, lmp, kmp, 0);
    vector<uint> covgs_exp = {0, 1, 1, 1};
    EXPECT_ITERABLE_EQ(vector<uint>, covgs_exp, covgs);

    l3->write_covgs_to_file("localPRG_test.covgs", covgs);
}

TEST(LocalPRGTest, write_path_to_fasta) {
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");

    auto index = std::make_shared<Index>();
    l3.minimizer_sketch(index, 1, 3);

    vector<LocalNodePtr> lmp3 = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};
    l3.write_path_to_fasta("localPRG_test.maxpath.fa", lmp3, 0.00);
}

TEST(LocalPRGTest, append_path_to_fasta) {
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");

    auto index = std::make_shared<Index>();
    l3.minimizer_sketch(index, 1, 3);

    vector<LocalNodePtr> lmp3 = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};
    l3.append_path_to_fasta("localPRG_test.maxpath.fa", lmp3, 0.00);


}

TEST(LocalPRGTest, write_aligned_path_to_fasta) {
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");

    auto index = std::make_shared<Index>();
    l3.minimizer_sketch(index, 1, 3);

    vector<LocalNodePtr> lmp3 = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};
    l3.write_aligned_path_to_fasta("localPRG_test.alignedpath.fa", lmp3, 0.00);


}

TEST(LocalPRGTest, build_vcf) {
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");
    LocalPRG l4(4, "small real PRG", "ATGACAAAACGAAGTGGAAGTAATACGCGCAGGCGGGCTATCAGTCGCCCTGTTCGTCTGACGGCAGAAGAAGACCAGG"
                                     "AAATCAGAAAAAGGGCTGCTGAATGCGGCAAGACCGTTTC 5 T 6 C 5 GGTTTTTTACGGGCGGCAGCTCTCGGTAAGAAAGTTAA 7 TTCACTGACTGA"
                                     "TGACCGAGTGCTGAAAGAAGTCATGCGACTGGGGGCGTTG 8 CTCACTGACTGATGATCGGGTACTGAAAGAAGTTATGAGACTGGGGGCGTTA 7 CAGAAA"
                                     "AAACTCTTTATCGACGGCAAGCGTGTCGGGGACAG 9 A 10 G 9 GAGTATGCGGAGGTGCTGAT 11 A 12 C 11 GCTATTACGGAGTATCACCG 13"
                                     " G 14 T 13 GCCCTGTTATCCAGGCTTATGGCAGATTAG");
    LocalPRG l5(5, "another real PRG",
                "ATGACAAAGGTTACACCGT 5 C 6 T 5 TGACGTGCTACGCCTGTCAGGCCTATTCGACTCCTGCAAT 7 G 8 A 7 TATTGAATTTGCATAGTTTT 9 G 10 A 9 TAGGTCGA 11 G 12 A 11 TAAGGCGTTCACGCCGCATCCGGCGTGAACAAA 13 G 14 T 13 TACTCTTTTT 15  17  19 C 20 T 19 GCACAATCCAA 18 CGCACAAACCAA 17  16  21 CGCACAATCCAA 22  23 CGT 24 CGC 23 ACAAACCA 25 A 26 T 25  21 TATGTGCAAATTATTACTTTTTCCAGAAATCATCGAAAACGG 15 ");

    VCF vcf;

    l1.build_vcf(vcf, l1.prg.top_path());
    uint j = 0;
    EXPECT_EQ(j, vcf.records.size());
    EXPECT_EQ(j, vcf.samples.size());

    vcf = VCF();
    l2.build_vcf(vcf, l2.prg.top_path());
    j = 1;
    EXPECT_EQ(j, vcf.records.size());
    EXPECT_EQ("varsite", vcf.records[0]->chrom);
    EXPECT_EQ((uint) 1, vcf.records[0]->pos);
    EXPECT_EQ("GC", vcf.records[0]->ref);
    EXPECT_EQ("G", vcf.records[0]->alts[0]);
    EXPECT_EQ("SVTYPE=INDEL;GRAPHTYPE=SIMPLE", vcf.records[0]->info);

    vcf = VCF();
    vector<LocalNodePtr> lmp = {l2.prg.nodes[0], l2.prg.nodes[2], l2.prg.nodes[3]};
    l2.build_vcf(vcf, lmp);
    j = 1;
    EXPECT_EQ(j, vcf.records.size());
    EXPECT_EQ("varsite", vcf.records[0]->chrom);
    EXPECT_EQ((uint) 1, vcf.records[0]->pos);
    EXPECT_EQ("G", vcf.records[0]->ref);
    EXPECT_EQ("GC", vcf.records[0]->alts[0]);
    EXPECT_EQ("SVTYPE=INDEL;GRAPHTYPE=SIMPLE", vcf.records[0]->info);

    vcf = VCF();
    l3.build_vcf(vcf, l3.prg.top_path());
    vcf.sort_records();
    j = 2;
    EXPECT_EQ(j, vcf.records.size());
    EXPECT_EQ("nested varsite", vcf.records[0]->chrom);
    EXPECT_EQ((uint) 1, vcf.records[0]->pos);
    EXPECT_EQ("GC", vcf.records[0]->ref);
    EXPECT_EQ("G", vcf.records[0]->alts[0]);
    EXPECT_EQ("SVTYPE=INDEL;GRAPHTYPE=NESTED", vcf.records[0]->info);
    EXPECT_EQ((uint) 2, vcf.records[1]->pos);
    EXPECT_EQ("C", vcf.records[1]->ref);
    EXPECT_EQ("T", vcf.records[1]->alts[0]);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=NESTED", vcf.records[1]->info);

    vcf = VCF();
    lmp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};
    l3.build_vcf(vcf, lmp);
    vcf.sort_records();
    EXPECT_EQ(j, vcf.records.size());
    EXPECT_EQ("nested varsite", vcf.records[0]->chrom);
    EXPECT_EQ((uint) 1, vcf.records[0]->pos);
    EXPECT_EQ("GT", vcf.records[0]->ref);
    EXPECT_EQ("G", vcf.records[0]->alts[0]);
    EXPECT_EQ("SVTYPE=INDEL;GRAPHTYPE=NESTED", vcf.records[0]->info);
    EXPECT_EQ((uint) 2, vcf.records[1]->pos);
    EXPECT_EQ("T", vcf.records[1]->ref);
    EXPECT_EQ("C", vcf.records[1]->alts[0]);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=NESTED", vcf.records[1]->info);

    vcf = VCF();
    lmp = {l3.prg.nodes[0], l3.prg.nodes[5], l3.prg.nodes[6]};
    l3.build_vcf(vcf, lmp);
    vcf.sort_records();
    EXPECT_EQ(j, vcf.records.size());
    EXPECT_EQ("nested varsite", vcf.records[0]->chrom);
    EXPECT_EQ((uint) 1, vcf.records[0]->pos);
    EXPECT_EQ("G", vcf.records[0]->ref);
    EXPECT_EQ("GC", vcf.records[0]->alts[0]);
    EXPECT_EQ("SVTYPE=INDEL;GRAPHTYPE=SIMPLE", vcf.records[0]->info);
    EXPECT_EQ((uint) 1, vcf.records[1]->pos);
    EXPECT_EQ("G", vcf.records[1]->ref);
    EXPECT_EQ("GT", vcf.records[1]->alts[0]);
    EXPECT_EQ("SVTYPE=INDEL;GRAPHTYPE=SIMPLE", vcf.records[1]->info);

    vcf = VCF();
    l4.build_vcf(vcf, l4.prg.top_path());
    vcf.sort_records();
    j = 5;
    EXPECT_EQ(j, vcf.records.size());
    EXPECT_EQ("small real PRG", vcf.records[0]->chrom);
    EXPECT_EQ((uint) 119, vcf.records[0]->pos);
    EXPECT_EQ("T", vcf.records[0]->ref);
    EXPECT_EQ("C", vcf.records[0]->alts[0]);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", vcf.records[0]->info);

    EXPECT_EQ((uint) 158, vcf.records[1]->pos);
    EXPECT_EQ("TTCACTGACTGATGACCGAGTGCTGAAAGAAGTCATGCGACTGGGGGCGTTG", vcf.records[1]->ref);
    EXPECT_EQ("CTCACTGACTGATGATCGGGTACTGAAAGAAGTTATGAGACTGGGGGCGTTA", vcf.records[1]->alts[0]);
    EXPECT_EQ("SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE", vcf.records[1]->info);

    EXPECT_EQ((uint) 251, vcf.records[2]->pos);
    EXPECT_EQ("A", vcf.records[2]->ref);
    EXPECT_EQ("G", vcf.records[2]->alts[0]);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", vcf.records[2]->info);

    EXPECT_EQ((uint) 272, vcf.records[3]->pos);
    EXPECT_EQ("A", vcf.records[3]->ref);
    EXPECT_EQ("C", vcf.records[3]->alts[0]);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", vcf.records[3]->info);

    EXPECT_EQ((uint) 293, vcf.records[4]->pos);
    EXPECT_EQ("G", vcf.records[4]->ref);
    EXPECT_EQ("T", vcf.records[4]->alts[0]);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", vcf.records[4]->info);

    vcf = VCF();
    lmp = {l4.prg.nodes[0], l4.prg.nodes[2], l4.prg.nodes[3], l4.prg.nodes[4], l4.prg.nodes[6], l4.prg.nodes[8],
           l4.prg.nodes[9], l4.prg.nodes[10], l4.prg.nodes[12], l4.prg.nodes[14], l4.prg.nodes[15]};
    l4.build_vcf(vcf, lmp);
    vcf.sort_records();
    j = 5;
    EXPECT_EQ(j, vcf.records.size());
    EXPECT_EQ("small real PRG", vcf.records[0]->chrom);
    EXPECT_EQ((uint) 119, vcf.records[0]->pos);
    EXPECT_EQ("C", vcf.records[0]->ref);
    EXPECT_EQ("T", vcf.records[0]->alts[0]);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", vcf.records[0]->info);

    EXPECT_EQ((uint) 158, vcf.records[1]->pos);
    EXPECT_EQ("TTCACTGACTGATGACCGAGTGCTGAAAGAAGTCATGCGACTGGGGGCGTTG", vcf.records[1]->ref);
    EXPECT_EQ("CTCACTGACTGATGATCGGGTACTGAAAGAAGTTATGAGACTGGGGGCGTTA", vcf.records[1]->alts[0]);
    EXPECT_EQ("SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE", vcf.records[1]->info);

    EXPECT_EQ((uint) 251, vcf.records[2]->pos);
    EXPECT_EQ("G", vcf.records[2]->ref);
    EXPECT_EQ("A", vcf.records[2]->alts[0]);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", vcf.records[2]->info);

    EXPECT_EQ((uint) 272, vcf.records[3]->pos);
    EXPECT_EQ("A", vcf.records[3]->ref);
    EXPECT_EQ("C", vcf.records[3]->alts[0]);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", vcf.records[3]->info);

    EXPECT_EQ((uint) 293, vcf.records[4]->pos);
    EXPECT_EQ("T", vcf.records[4]->ref);
    EXPECT_EQ("G", vcf.records[4]->alts[0]);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", vcf.records[4]->info);

    vcf = VCF();
    l5.build_vcf(vcf, l5.prg.top_path());
    vcf.sort_records();
}

TEST(LocalPRGTest, build_vcf_real) {
    LocalPRG l1(1, "GC00000008_13", "ATGAAATTAAAAATAGTT 5 G 6 A 5 CGGTGGTTGTAACTGGTTTGTTAGCTGCGAACGTAGC 7 ACACGCT 8 ACACGCC 8 GCACGCT 7 GCCGAAGTC 9 T 10 G 9 ATAACAAGGATGGTAATAAACT 11 C 12 A 11 GACCTTTATGGCAAGGTTACCGCTCTACGTTATTTTACTGATGATAAGCGTGA 13 C 14 A 13 GATGGTGATAAAA 15 C 16 T 16 A 15 TTATGCCCGTCTCGGCTTTAAAGGAGAAACG 17 CAAATCAATGATCAAATGA 19 C 20 T 19 TGGTTTTGGTCACTGGGAATATGATTTTAAAGGCTATAACGATGAAGCCAACGG 21 C 22 T 21 TCGCGCG 23 AC 24 GC 24 GT 23 AACAAGACCCGTCT 25 G 26 A 25 GCCTATGC 27 T 28 A 27 GGTTTAAAAATTAGTGAATTTGGCTCTCTGGACTATGG 29 C 30 T 29 CGTAACTACGGTGTCGGCTATGACATTGGTTCATGGAC 31 CT 32 CG 32 TG 31 ATATGTTGCCAGAATTTGGTGGCGATACCTGGAGTCAGAAAGATGTCTTCATGACATACCGTAC 33 C 34 T 34 A 33 ACCGGTGT 35 G 36 A 35 GCAACCTATCGCAACTACGATTTCTT 37 C 38 T 37 GGCTTAATTGAAGG 39 T 40 G 39 CTGAACTTTGCCGCGCAATATCAAGGCAAAAATGAACG 41 C 42 T 41 ACTGACAA 43 TGGT 44 CAGT 44 TGGC 43 CATCTTTATGGTGCTGACTA 45  47 C 48 T 47 ACGCGTGCCAA 49 C 50 T 49  46 CACGCGCGCCAAC 45 GGTGACGGTTTCGGTATCTCCTCAACTTATGTTTATGATGGCTTTGGTATCGG 51 TGCGGTA 52 AGCGGTG 52 TGCGGTG 51 TATACCAAATCCGATCGGACAA 53 ATGCA 54 TTGCA 54 ATGCG 53 CAGGAAAGAGCCGCTGCTAATCCTCTCAATGCCTCCGGTAAGAATGCAGAACTGTGGGC 55 C 56 T 55 ACAGGTATAAA 57 G 58 A 57 TATGATGCCAACAACAT 59 C 60 T 59 TACTTTGCAGCTAATTACGCTGAAACATTAAACATGACCACCTATGG 61 C 62 G 61 GATGGTTATAT 63 CTCT 64 CTCG 64 TTCT 63 AACAAAGCACAAAGTTTTGAAGT 65 AGTGACA 66 AGTGGCA 66 GGTGGCG 65 CAATATCAATTCGACTTCGGCTTGC 67 GCCCA 68 ACCCC 68 GCCCC 67 TCACTCGCTTACCTGAAATCGAAAGGCA 69 T 70 G 69 AGATCTGGGCCGCTACGGCGA 71  73 C 74 T 73 CAGGACATGATTGAGTATATCGACGTTGGTGCGACGTATTTCTTCAACAAAAATATGTCGACCTATGTTGATTATAAAATCAACCTGATTGATGAAAGCGACTTTACCCGTGCCGTAGATATTCGCAC 75 CGATAACATCGTCGC 77 AACGGGT 78 AACGGGC 78 TACGGGC 78 AACGGGA 77 ATTACCTATCAGTTC 76 GGCTTTGTTGAATAAATCGAACTTTTGC 75  72 AGTTAACGGCATCAACAATGATCCACTTGCCCAAATGCAGTACTGGACTGCAGTAAGAAATATAATTGATGACACTAATGAAGTGACCATTGAATTATCTTATAACCTGGCAATCACAAATATCGATACCAGCGATGAACATCTTGTAGAAGTAAGCGAGAATTCCGAAGGAAATCATATAAAAGACAATGACTCAATGTCTATTCGTTATAGATCAAAATATTATTCCAGAGAGTACGCTTTAATAGAAGAAGAAACAATATTTTCTGACGCAGAACTAAAAGCCATTCTGCCTATGCATCGCATGTACGGGGTTGGTGACTATAAGTCAAATTCCTCTTCTCTACCCTCACACTCGGGGCTAAAGGACCCAACGGGCACACCCGTCTGTTATTATATTCATAATGAGGATAAACCTTCCTTAGGTTTTGGTCCAATATCCAATAATTGGTTAAGCCAATCCTTTACAACAGAGTTA 71  18  18 CAAATCAATG 79 T 80 A 79 TCAAATG 81 CTTG 82 ATTA 82 ATTG 81 GTTTTGGTCAC 83 T 84 A 83 GGGAATATGATTTTAAAGGCTATAACGATGAAGCCAAC 85 GGCTCGC 86 GGCTCGT 86 TGCTCGC 86 GGTTCGC 85 GCGGCAA 87  89 CAATCTC 90 GAATCTC 89  88 GAATCTCTTATTGAGT 87  17 ");

    VCF vcf;
    auto ref_path = l1.prg.top_path();
    l1.build_vcf(vcf, ref_path);
    auto ref_seq = l1.string_along_path(ref_path);

    vcf.correct_dot_alleles(ref_seq, "GC00000008_13");
}

TEST(LocalPRGTest, add_sample_gt_to_vcf) {
    LocalPRG l1(1, "simple", "AGCT");
    LocalPRG l2(2, "varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");
    LocalPRG l4(4, "small real PRG", "ATGACAAAACGAAGTGGAAGTAATACGCGCAGGCGGGCTATCAGTCGCCCTGTTCGTCTGACGGCAGAAGAAGACCAGG"
                                     "AAATCAGAAAAAGGGCTGCTGAATGCGGCAAGACCGTTTC 5 T 6 C 5 GGTTTTTTACGGGCGGCAGCTCTCGGTAAGAAAGTTAA 7 TTCACTGACTGA"
                                     "TGACCGAGTGCTGAAAGAAGTCATGCGACTGGGGGCGTTG 8 CTCACTGACTGATGATCGGGTACTGAAAGAAGTTATGAGACTGGGGGCGTTA 7 CAGAAA"
                                     "AAACTCTTTATCGACGGCAAGCGTGTCGGGGACAG 9 A 10 G 9 GAGTATGCGGAGGTGCTGAT 11 A 12 C 11 GCTATTACGGAGTATCACCG 13"
                                     " G 14 T 13 GCCCTGTTATCCAGGCTTATGGCAGATTAG");
    LocalPRG l5(5, "another real PRG", " 5 ATGCTTATTGGCTATGT 7  9 ACGCGTA 10 TCGCGTA 10 ACGTGTG 9 TCAACAAATGACCAGAACA"
                                       "C 11 A 12 C 11  8 ACGCGTATCAACAAATGATCAGAACACA 7 GATCTACAACGTAATGCG 6 AAGT 5 ");

    VCF vcf;

    vector<LocalNodePtr> lmp1 = {l1.prg.nodes[0]};
    l1.build_vcf(vcf, l1.prg.top_path());
    l1.add_sample_gt_to_vcf(vcf, l1.prg.top_path(), lmp1, "sample");
    uint j = 1;
    EXPECT_EQ(j, vcf.samples.size());

    vcf = VCF();
    vector<LocalNodePtr> lmp2 = {l2.prg.nodes[0], l2.prg.nodes[2], l2.prg.nodes[3]};
    l2.build_vcf(vcf, l2.prg.top_path());
    l2.add_sample_gt_to_vcf(vcf, l2.prg.top_path(), lmp2, "sample");
    j = 1;
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);

    vcf = VCF();
    vector<LocalNodePtr> lmp3 = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};
    l3.build_vcf(vcf, l3.prg.top_path());
    vcf.sort_records();
    l3.add_sample_gt_to_vcf(vcf, l3.prg.top_path(), lmp3, "sample");
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);

    vcf = VCF();
    vector<LocalNodePtr> lmp4 = {l4.prg.nodes[0], l4.prg.nodes[1], l4.prg.nodes[3], l4.prg.nodes[5], l4.prg.nodes[6],
                                 l4.prg.nodes[8], l4.prg.nodes[9], l4.prg.nodes[10], l4.prg.nodes[12], l4.prg.nodes[13],
                                 l4.prg.nodes[15]};
    l4.build_vcf(vcf, l4.prg.top_path());
    vcf.sort_records();
    l4.add_sample_gt_to_vcf(vcf, l4.prg.top_path(), lmp4, "sample");
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[2]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 1, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[3]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[4]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);

    vcf = VCF();
    vector<LocalNodePtr> lmp5 = {l5.prg.nodes[0], l5.prg.nodes[1], l5.prg.nodes[10], l5.prg.nodes[11],
                                 l5.prg.nodes[13]};
    l5.build_vcf(vcf, l5.prg.top_path());
    vcf.sort_records();
    l5.add_sample_gt_to_vcf(vcf, l5.prg.top_path(), lmp5, "sample");
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ((uint) 5, vcf.records.size());
    EXPECT_EQ(j, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_TRUE(vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].end());
    EXPECT_EQ(j, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_TRUE(vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].end());
    EXPECT_EQ(j, vcf.records[2]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_TRUE(vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].end());
    EXPECT_EQ(j, vcf.records[3]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 1, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[4]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_TRUE(vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].end());

    // add the ref path
    l5.add_sample_gt_to_vcf(vcf, l5.prg.top_path(), l5.prg.top_path(), "sample2");
    EXPECT_EQ((uint) 2, vcf.samples.size());
    EXPECT_EQ((uint) 5, vcf.records.size());
    EXPECT_EQ((uint) 2, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
    EXPECT_EQ((uint) 2, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
    EXPECT_EQ((uint) 2, vcf.records[2]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
    EXPECT_EQ((uint) 2, vcf.records[3]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
    EXPECT_EQ((uint) 2, vcf.records[4]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);

}

TEST(LocalPRGTest, moreupdateVCF) {
    // load PRGs from file
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, "../../test/test_cases/updatevcf_test.fa");

    EXPECT_EQ((uint) 3, prgs.size());

    VCF vcf;

    prgs[0]->build_vcf(vcf, prgs[0]->prg.top_path());
    prgs[1]->build_vcf(vcf, prgs[1]->prg.top_path());
    prgs[2]->build_vcf(vcf, prgs[2]->prg.top_path());
    vcf.sort_records();

    vector<LocalNodePtr> lmp1 = {prgs[1]->prg.nodes[0], prgs[1]->prg.nodes[11], prgs[1]->prg.nodes[12],
                                 prgs[1]->prg.nodes[17], prgs[1]->prg.nodes[65], prgs[1]->prg.nodes[67]};
    prgs[1]->add_sample_gt_to_vcf(vcf, prgs[1]->prg.top_path(), lmp1, "sample");

    vector<LocalNodePtr> lmp2 = {prgs[2]->prg.nodes[0], prgs[2]->prg.nodes[1], prgs[2]->prg.nodes[3],
                                 prgs[2]->prg.nodes[4], prgs[2]->prg.nodes[6], prgs[2]->prg.nodes[7],
                                 prgs[2]->prg.nodes[9], prgs[2]->prg.nodes[10], prgs[2]->prg.nodes[11],
                                 prgs[2]->prg.nodes[13], prgs[2]->prg.nodes[14], prgs[2]->prg.nodes[16],
                                 prgs[2]->prg.nodes[17], prgs[2]->prg.nodes[19], prgs[2]->prg.nodes[44],
                                 prgs[2]->prg.nodes[45], prgs[2]->prg.nodes[47], prgs[2]->prg.nodes[118],
                                 prgs[2]->prg.nodes[119], prgs[2]->prg.nodes[121], prgs[2]->prg.nodes[123],
                                 prgs[2]->prg.nodes[125], prgs[2]->prg.nodes[126], prgs[2]->prg.nodes[130],
                                 prgs[2]->prg.nodes[131], prgs[2]->prg.nodes[133], prgs[2]->prg.nodes[135],
                                 prgs[2]->prg.nodes[141], prgs[2]->prg.nodes[142], prgs[2]->prg.nodes[144],
                                 prgs[2]->prg.nodes[145], prgs[2]->prg.nodes[160]};
    prgs[2]->add_sample_gt_to_vcf(vcf, prgs[2]->prg.top_path(), lmp2, "sample");
}

TEST(LocalPRGTest, find_alt_path) {
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT 9 T 10  9 ATG");

    vector<LocalNodePtr> top = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6]};
    vector<LocalNodePtr> middle = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};
    vector<LocalNodePtr> bottom = {l3.prg.nodes[0], l3.prg.nodes[5], l3.prg.nodes[6]};

    vector<LocalNodePtr> alt_path = l3.find_alt_path(top, 2, "C", "T");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, middle, alt_path);

    alt_path = l3.find_alt_path(top, 1, "GC", "G");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, bottom, alt_path);

    alt_path = l3.find_alt_path(middle, 2, "T", "C");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, top, alt_path);

    alt_path = l3.find_alt_path(top, 1, "GT", "G");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, bottom, alt_path);

    alt_path = l3.find_alt_path(bottom, 1, "G", "GT");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, middle, alt_path);

    alt_path = l3.find_alt_path(bottom, 1, "G", "GC");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, top, alt_path);

    // and now for the one where the alts or ref is "."
    top = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6], l3.prg.nodes[7],
           l3.prg.nodes[9]};
    bottom = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6], l3.prg.nodes[8],
              l3.prg.nodes[9]};
    alt_path = l3.find_alt_path(top, 6, "T", ".");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, bottom, alt_path);

    alt_path = l3.find_alt_path(bottom, 6, ".", "T");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, top, alt_path);

    // if the site is at the start and alts is "."
    LocalPRG l3_(3, "nested varsite", " 5 G 7 C 8 T 7  6  5 TAT 9 T 10  9 ");
    top = {l3_.prg.nodes[0], l3_.prg.nodes[1], l3_.prg.nodes[2], l3_.prg.nodes[4], l3_.prg.nodes[6]};
    bottom = {l3_.prg.nodes[0], l3_.prg.nodes[5], l3_.prg.nodes[6]};

    alt_path = l3_.find_alt_path(top, 0, "GC", ".");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, bottom, alt_path);

    alt_path = l3_.find_alt_path(bottom, 0, ".", "GC");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, top, alt_path);

    // if the site at the end has ref/alts as "."
    top = {l3_.prg.nodes[0], l3_.prg.nodes[1], l3_.prg.nodes[2], l3_.prg.nodes[4],
           l3_.prg.nodes[6], l3_.prg.nodes[7], l3_.prg.nodes[9]};
    bottom = {l3_.prg.nodes[0], l3_.prg.nodes[1], l3_.prg.nodes[2], l3_.prg.nodes[4],
              l3_.prg.nodes[6], l3_.prg.nodes[8], l3_.prg.nodes[9]};

    alt_path = l3_.find_alt_path(top, 5, "T", ".");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, bottom, alt_path);

    alt_path = l3_.find_alt_path(bottom, 5, ".", "T");
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, top, alt_path);

}

class LocalPRGMock : public LocalPRG {
public:
    LocalPRGMock() : LocalPRG(0, "", "") {}
};


class LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture : public ::testing::Test {
protected:
    class LocalPRGMockExposesTestedMethod : public LocalPRGMock {
    public:
        virtual uint32_t
        get_number_of_bases_in_local_path_before_a_given_position(const std::vector<LocalNodePtr> &local_path,
                                                                  uint32_t position) const {
            return LocalPRGMock::get_number_of_bases_in_local_path_before_a_given_position(local_path, position);
        }
    };


    void SetUp() override {
        local_path_with_two_intervals.push_back(std::make_shared<LocalNode>("", Interval(10, 20), 0));
        local_path_with_two_intervals.push_back(std::make_shared<LocalNode>("", Interval(35, 45), 1));
    }

    void TearDown() override {
    }

    LocalPRGMockExposesTestedMethod local_prg_mock;
    std::vector<LocalNodePtr> empty_local_path;
    std::vector<LocalNodePtr> local_path_with_two_intervals;
};


TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture,
       empty_local_path___returns_0) {
    uint32_t position = 30;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(empty_local_path,
                                                                                               position);

    uint32_t expected = 0;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture,
       local_path_with_two_intervals___position_just_before_first_interval___returns_0) {
    uint32_t position = 9;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(
            local_path_with_two_intervals, position);

    uint32_t expected = 0;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture,
       local_path_with_two_intervals___position_at_the_start_of_first_interval___returns_0) {
    uint32_t position = 10;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(
            local_path_with_two_intervals, position);

    uint32_t expected = 0;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture,
       local_path_with_two_intervals___position_just_after_the_start_of_first_interval___returns_1) {
    uint32_t position = 11;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(
            local_path_with_two_intervals, position);

    uint32_t expected = 1;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture,
       local_path_with_two_intervals___position_just_before_the_end_of_first_interval___returns_8) {
    uint32_t position = 18;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(
            local_path_with_two_intervals, position);

    uint32_t expected = 8;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture,
       local_path_with_two_intervals___position_at_the_end_of_first_interval___returns_9) {
    uint32_t position = 19;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(
            local_path_with_two_intervals, position);

    uint32_t expected = 9;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture,
       local_path_with_two_intervals___position_just_after_the_end_of_first_interval___returns_10) {
    uint32_t position = 20;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(
            local_path_with_two_intervals, position);

    uint32_t expected = 10;
    ASSERT_EQ(actual, expected);
}


TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture,
       local_path_with_two_intervals___position_between_intervals___returns_10) {
    uint32_t position = 30;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(
            local_path_with_two_intervals, position);

    uint32_t expected = 10;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture,
       local_path_with_two_intervals___position_at_the_start_of_second_interval___returns_10) {
    uint32_t position = 35;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(
            local_path_with_two_intervals, position);

    uint32_t expected = 10;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_in_local_path_before_a_given_position___Fixture,
       local_path_with_two_intervals___position_just_after_the_start_of_second_interval___returns_11) {
    uint32_t position = 36;

    uint32_t actual = local_prg_mock.get_number_of_bases_in_local_path_before_a_given_position(
            local_path_with_two_intervals, position);

    uint32_t expected = 11;
    ASSERT_EQ(actual, expected);
}


class LocalPRGTest___get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node___Fixture
        : public ::testing::Test {
protected:
    class LocalPRGMockExposesTestedMethod : public LocalPRGMock {
    public:
        virtual uint32_t
        get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node(const KmerNodePtr &previous_kmer_node,
                                                                           const KmerNodePtr &current_kmer_node) const {
            return LocalPRGMock::get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node(previous_kmer_node,
                                                                                                    current_kmer_node);
        }
    };


    void SetUp() override {
        {
            prg::Path path_from_3_to_30;
            path_from_3_to_30.push_back(Interval(3, 7));
            path_from_3_to_30.push_back(Interval(12, 20));
            path_from_3_to_30.push_back(Interval(20, 25));
            path_from_3_to_30.push_back(Interval(25, 30));
            kmer_node_from_3_to_30 = std::make_shared<KmerNode>(1, path_from_3_to_30);
        }

        {
            prg::Path path_from_3_to_50;
            path_from_3_to_50.push_back(Interval(3, 50));
            kmer_node_from_3_to_50 = std::make_shared<KmerNode>(2, path_from_3_to_50);
        }

        {
            prg::Path path_from_4_to_50;
            path_from_4_to_50.push_back(Interval(4, 50));
            kmer_node_from_4_to_50 = std::make_shared<KmerNode>(2, path_from_4_to_50);
        }

        {
            prg::Path path_from_7_to_50;
            path_from_7_to_50.push_back(Interval(7, 50));
            kmer_node_from_7_to_50 = std::make_shared<KmerNode>(2, path_from_7_to_50);
        }

        {
            prg::Path path_from_10_to_50;
            path_from_10_to_50.push_back(Interval(10, 50));
            kmer_node_from_10_to_50 = std::make_shared<KmerNode>(2, path_from_10_to_50);
        }

        {
            prg::Path path_from_15_to_50;
            path_from_15_to_50.push_back(Interval(15, 50));
            kmer_node_from_15_to_50 = std::make_shared<KmerNode>(2, path_from_15_to_50);
        }


        {
            prg::Path path_from_40_to_50;
            path_from_40_to_50.push_back(Interval(40, 50));
            kmer_node_from_40_to_50 = std::make_shared<KmerNode>(2, path_from_40_to_50);
        }
    }

    void TearDown() override {
    }

    LocalPRGMockExposesTestedMethod local_prg_mock;
    KmerNodePtr kmer_node_from_3_to_30;
    KmerNodePtr kmer_node_from_3_to_50;
    KmerNodePtr kmer_node_from_4_to_50;
    KmerNodePtr kmer_node_from_7_to_50;
    KmerNodePtr kmer_node_from_10_to_50;
    KmerNodePtr kmer_node_from_15_to_50;
    KmerNodePtr kmer_node_from_40_to_50;
};


TEST_F(LocalPRGTest___get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node___Fixture,
       previous_node_exactly_at_the_start_of_current_node___returns_0) {
    KmerNodePtr previous_node = kmer_node_from_3_to_30;
    KmerNodePtr current_node = kmer_node_from_3_to_50;

    uint32_t actual = local_prg_mock.get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node(previous_node,
                                                                                                        current_node);

    uint32_t expected = 0;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node___Fixture,
       previous_node_just_before_the_start_of_current_node___returns_1) {
    KmerNodePtr previous_node = kmer_node_from_3_to_30;
    KmerNodePtr current_node = kmer_node_from_4_to_50;

    uint32_t actual = local_prg_mock.get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node(previous_node,
                                                                                                        current_node);

    uint32_t expected = 1;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node___Fixture,
       current_node_right_at_the_end_of_first_interval_of_previous_node___returns_4) {
    KmerNodePtr previous_node = kmer_node_from_3_to_30;
    KmerNodePtr current_node = kmer_node_from_7_to_50;

    uint32_t actual = local_prg_mock.get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node(previous_node,
                                                                                                        current_node);

    uint32_t expected = 4;
    ASSERT_EQ(actual, expected);
}

TEST_F(LocalPRGTest___get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node___Fixture,
       current_node_between_the_first_and_second_interval_of_previous_node___returns_4) {
    KmerNodePtr previous_node = kmer_node_from_3_to_30;
    KmerNodePtr current_node = kmer_node_from_10_to_50;

    uint32_t actual = local_prg_mock.get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node(previous_node,
                                                                                                        current_node);

    uint32_t expected = 4;
    ASSERT_EQ(actual, expected);
}


TEST_F(LocalPRGTest___get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node___Fixture,
       current_node_inside_second_interval_of_previous_node___returns_7) {
    KmerNodePtr previous_node = kmer_node_from_3_to_30;
    KmerNodePtr current_node = kmer_node_from_15_to_50;

    uint32_t actual = local_prg_mock.get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node(previous_node,
                                                                                                        current_node);

    uint32_t expected = 7;
    ASSERT_EQ(actual, expected);
}


TEST_F(LocalPRGTest___get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node___Fixture,
       current_node_after_previous_node___returns_22) {
    KmerNodePtr previous_node = kmer_node_from_3_to_30;
    KmerNodePtr current_node = kmer_node_from_40_to_50;

    uint32_t actual = local_prg_mock.get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node(previous_node,
                                                                                                        current_node);

    uint32_t expected = 22;
    ASSERT_EQ(actual, expected);
}



TEST(LocalPRGTest, get_forward_and_reverse_kmer_coverages_in_range) {
    auto index = std::make_shared<Index>();
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");
    l3.minimizer_sketch(index, 1, 3);
    KmerGraphWithCoverage kg(&l3.kmer_prg);

    kg.set_covg(2, 4, 0, 0);
    kg.set_covg(2, 3, 1, 0);
    kg.set_covg(5, 4, 0, 0);
    kg.set_covg(5, 5, 1, 0);
    kg.set_covg(7, 2, 0, 0);
    kg.set_covg(7, 3, 1, 0);
    kg.set_covg(8, 4, 0, 0);
    kg.set_covg(8, 6, 1, 0);

    vector<LocalNodePtr> lmp = {};
    vector<KmerNodePtr> kmp = {
            l3.kmer_prg.nodes[0],
            l3.kmer_prg.nodes[2],
            l3.kmer_prg.nodes[5],
            l3.kmer_prg.nodes[8],
            l3.kmer_prg.nodes[10],
            l3.kmer_prg.nodes[11]
    };
    vector<uint32_t> fwd, rev, exp_fwd, exp_rev;

    std::tie(fwd, rev) = l3.get_forward_and_reverse_kmer_coverages_in_range(kg, kmp, lmp, 0, 0, 0);
    EXPECT_TRUE(fwd.empty());
    EXPECT_TRUE(rev.empty());

    std::tie(fwd, rev) = l3.get_forward_and_reverse_kmer_coverages_in_range(kg, kmp, lmp, 0, 1, 0);
    exp_fwd = {4};
    exp_rev = {3};
    EXPECT_ITERABLE_EQ(vector<uint32_t>, exp_fwd, fwd);
    EXPECT_ITERABLE_EQ(vector<uint32_t>, exp_rev, rev);

    std::tie(fwd, rev) = l3.get_forward_and_reverse_kmer_coverages_in_range(kg, kmp, lmp, 0, 2, 0);
    exp_fwd = {4, 4};
    exp_rev = {3, 5};
    EXPECT_ITERABLE_EQ(vector<uint32_t>, exp_fwd, fwd);
    EXPECT_ITERABLE_EQ(vector<uint32_t>, exp_rev, rev);

    std::tie(fwd, rev) = l3.get_forward_and_reverse_kmer_coverages_in_range(kg, kmp, lmp, 0, 3, 0);
    exp_fwd = {4, 4, 4};
    exp_rev = {3, 5, 6};
    EXPECT_ITERABLE_EQ(vector<uint32_t>, exp_fwd, fwd);
    EXPECT_ITERABLE_EQ(vector<uint32_t>, exp_rev, rev);

    std::tie(fwd, rev) = l3.get_forward_and_reverse_kmer_coverages_in_range(kg, kmp, lmp, 1, 2, 0);
    exp_fwd = {4, 4};
    exp_rev = {3, 5};
    EXPECT_ITERABLE_EQ(vector<uint32_t>, exp_fwd, fwd);
    EXPECT_ITERABLE_EQ(vector<uint32_t>, exp_rev, rev);
}

TEST(LocalPRGTest, add_sample_covgs_to_vcf) {
    uint32_t min_kmer_covgs = 0;
    auto index = std::make_shared<Index>();
    vector<string> short_formats = {"GT"};
    vector<string> formats = {"GT", "MEAN_FWD_COVG", "MEAN_REV_COVG",
                              "MED_FWD_COVG", "MED_REV_COVG",
                              "SUM_FWD_COVG", "SUM_REV_COVG", "GAPS"};

    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");
    l3.minimizer_sketch(index, 1, 3);

    VCF vcf;

    vector<LocalNodePtr> lmp3 = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};
    l3.build_vcf(vcf, l3.prg.top_path());
    vcf.sort_records();
    l3.add_sample_gt_to_vcf(vcf, l3.prg.top_path(), lmp3, "sample");
    EXPECT_EQ((uint) 1, vcf.samples.size());
    EXPECT_EQ((uint) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_ITERABLE_EQ(vector<string>, short_formats, vcf.records[0]->format);
    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);

    KmerGraphWithCoverage kg(&l3.kmer_prg);
    l3.add_sample_covgs_to_vcf(vcf, kg, l3.prg.top_path(), min_kmer_covgs, "sample",
                               0);
    EXPECT_EQ((uint) 1, vcf.samples.size());
    EXPECT_EQ((uint) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_ITERABLE_EQ(vector<string>, formats, vcf.records[0]->format);
    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"][1]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"][1]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MED_FWD_COVG"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MED_REV_COVG"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MED_FWD_COVG"][1]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MED_REV_COVG"][1]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["SUM_FWD_COVG"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["SUM_REV_COVG"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["SUM_FWD_COVG"][1]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["SUM_REV_COVG"][1]);

    // ref
    kg.set_covg(1, 1, 0, 0);
    kg.set_covg(1, 0, 1, 0);
    kg.set_covg(4, 1, 0, 0);
    kg.set_covg(4, 0, 1, 0);
    kg.set_covg(7, 1, 0, 0);
    kg.set_covg(7, 0, 1, 0);

    // alts
    kg.set_covg(2, 6, 0, 0);
    kg.set_covg(2, 8, 1, 0);
    kg.set_covg(5, 5, 0, 0);
    kg.set_covg(5, 5, 1, 0);
    kg.set_covg(8, 4, 0, 0);
    kg.set_covg(8, 5, 1, 0);

    l3.add_sample_covgs_to_vcf(vcf, kg, l3.prg.top_path(), min_kmer_covgs, "sample", 0);
    EXPECT_EQ((uint) 1, vcf.samples.size());
    EXPECT_EQ((uint) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_ITERABLE_EQ(vector<string>, formats, vcf.records[0]->format);
    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"][0]);
    EXPECT_EQ((uint16_t) 5, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"][1]);
    EXPECT_EQ((uint16_t) 6, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"][1]);
    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MED_FWD_COVG"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MED_REV_COVG"][0]);
    EXPECT_EQ((uint16_t) 5, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MED_FWD_COVG"][1]);
    EXPECT_EQ((uint16_t) 5, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MED_REV_COVG"][1]);
    EXPECT_EQ((uint16_t) 3, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["SUM_FWD_COVG"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["SUM_REV_COVG"][0]);
    EXPECT_EQ((uint16_t) 15, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["SUM_FWD_COVG"][1]);
    EXPECT_EQ((uint16_t) 18, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["SUM_REV_COVG"][1]);
}

TEST(LocalPRGTest, add_consensus_path_to_fastaq_bin) {
    auto index = std::make_shared<Index>();

    auto l3 { std::make_shared<LocalPRG>(3, "three", "A 5 G 7 C 8 T 7  6 G 5 TAT") };
    l3->minimizer_sketch(index, 1, 3);

    shared_ptr<pangenome::Node> pn3(make_shared<pangenome::Node>(l3));
    pn3->kmer_prg_with_coverage.set_covg(2, 4, 0, 0);
    pn3->kmer_prg_with_coverage.set_covg(2, 3, 1, 0);
    pn3->kmer_prg_with_coverage.set_covg(5, 4, 0, 0);
    pn3->kmer_prg_with_coverage.set_covg(5, 5, 0, 0);
    pn3->kmer_prg_with_coverage.set_covg(7, 2, 0, 0);
    pn3->kmer_prg_with_coverage.set_covg(7, 3, 1, 0);
    pn3->kmer_prg_with_coverage.set_covg(8, 4, 0, 0);
    pn3->kmer_prg_with_coverage.set_covg(8, 6, 0, 0);

    pn3->kmer_prg_with_coverage.set_num_reads(6);
    pn3->kmer_prg_with_coverage.set_binomial_parameter_p(0.0001);
    shared_ptr<pangenome::Read> pr(make_shared<pangenome::Read>(0));
    pn3->reads.insert(pr);

    Fastaq fq(false, true);
    vector<KmerNodePtr> kmp;
    vector<LocalNodePtr> lmp;

    uint32_t max_num_kmers_to_average = 100;
    l3->add_consensus_path_to_fastaq(fq, pn3, kmp, lmp, 1, true, 8, max_num_kmers_to_average, 0);
    EXPECT_EQ("AGTTAT", l3->string_along_path(lmp));
    bool added_to_fq = find(fq.names.begin(), fq.names.end(), "three") != fq.names.end();
    EXPECT_TRUE(added_to_fq);
    bool added_to_seqs = fq.sequences.find("three") != fq.sequences.end();
    EXPECT_TRUE(added_to_seqs);
    bool added_to_scores = fq.scores.find("three") != fq.scores.end();
    EXPECT_TRUE(added_to_scores);
    bool added_to_headers = fq.headers.find("three") != fq.headers.end();
    EXPECT_TRUE(added_to_headers);
    EXPECT_EQ("AGTTAT", fq.sequences["three"]);
    EXPECT_EQ(fq.scores["three"], "DDD\?\?!");

}

TEST(LocalPRGTest, add_consensus_path_to_fastaq_nbin) {
    auto index = std::make_shared<Index>();
    auto l3 { std::make_shared<LocalPRG>(3, "three", "A 5 G 7 C 8 T 7  6 G 5 TAT") };
    l3->minimizer_sketch(index, 1, 3);

    shared_ptr<pangenome::Node> pn3(make_shared<pangenome::Node>(l3));
    pn3->kmer_prg_with_coverage.set_covg(2 ,4, 0, 0);
    pn3->kmer_prg_with_coverage.set_covg(2, 3, 1, 0);
    pn3->kmer_prg_with_coverage.set_covg(5, 5, 0, 0);
    pn3->kmer_prg_with_coverage.set_covg(7, 2, 0, 0);
    pn3->kmer_prg_with_coverage.set_covg(7, 3, 1, 0);
    pn3->kmer_prg_with_coverage.set_covg(8, 6, 0, 0);
    pn3->kmer_prg_with_coverage.set_num_reads(6);
    pn3->kmer_prg_with_coverage.set_negative_binomial_parameters(0.05, 2.0);
    shared_ptr<pangenome::Read> pr(make_shared<pangenome::Read>(0));
    pn3->reads.insert(pr);

    Fastaq fq(false, true);
    vector<KmerNodePtr> kmp;
    vector<LocalNodePtr> lmp;

    uint32_t max_num_kmers_to_average = 100;
    l3->add_consensus_path_to_fastaq(fq, pn3, kmp, lmp, 1, false, 8, max_num_kmers_to_average, 0);

    EXPECT_NE(kmp.size(), 0);
    std::vector<uint32_t> expected = {2, 5, 8, 10};
    std::vector<uint32_t> result = {};
    for (const auto &kmer_node_ptr: kmp)
        result.push_back(kmer_node_ptr->id);

    EXPECT_ITERABLE_EQ(std::vector<uint32_t>, expected , result);

    EXPECT_EQ("AGTTAT", l3->string_along_path(lmp));
    bool added_to_fq = find(fq.names.begin(), fq.names.end(), "three") != fq.names.end();
    EXPECT_TRUE(added_to_fq);
    bool added_to_seqs = fq.sequences.find("three") != fq.sequences.end();
    EXPECT_TRUE(added_to_seqs);
    bool added_to_scores = fq.scores.find("three") != fq.scores.end();
    EXPECT_TRUE(added_to_scores);
    bool added_to_headers = fq.headers.find("three") != fq.headers.end();
    EXPECT_TRUE(added_to_headers);
    EXPECT_EQ("AGTTAT", fq.sequences["three"]);
    EXPECT_EQ(fq.scores["three"], "DDD\?\?!");
}

TEST(LocalPRGTest, get_valid_vcf_reference_real_example) {
    auto index = std::make_shared<Index>(Index());

    LocalPRG l(3, "GC00003042", " 5  7  9 ATGTTAGTTAGTAAAAGCAACGGATTTAACGCTAGCGCA 11 G 12 T 11 TTTTGGGTAGTGGAAGTTATAATGAAAATAAATCTTCTAAACAC 10 ATGTTAATTAATAAAAGCAACGGATTTAACGCTAGCGCAGTTTGGGGTAGTGGAAGTTATAATGAAAATAAATCTTCTAAACAC 10  9 ATGGAGCTACTAGCTCATAGTATT 13 T 14 G 13 TAAAATTAATTTGTAAGGAAGCTGCATCAGAGACGTATCGCGGTGCTCTTGAAA 15 CTTTACAAAAAATGAT 16 TTTTACAAAAAATAAT 16 CTTTACAAAAAATGAC 15 GTCTGAATGTATATATCA 17 AGAAGGCAACGCC 18 TGAAGGCAACGCC 17  8  19 ATGTTAGTTAGTAAAAACAACGAATTTAACACTAGCTCATTTATAGATAGTGGAAATTGTAATGAAAGAAAATCTTCTGAATCC 20  19 ATGGAGCTACTAGCTTATAGTATTATAAATTTAATTTGTAAGAAAGCTGCCTCAGGGACGTATCGCGGTGCTCTTGAAACTTTACAAAAAATGATGTCTGAATGTATATATCATGAAGGCAACGCC 8 ATGTTAGTTAGTAAAAGCAACGAAATTAACACTAGCACATTTATAAATAGTGGAA 21 A 22 G 21 TTGTAATGAAAGAAAATCTTCTAAATCCATGGAGTTACTAGCTTATAGTATTATAAAATTAATTTGTAAAGAAGCTGCCTCAGGGACGTATAACGGTGCTCTTGAAATTTTACAAAAAATGATGTCTGAATGTAAATATCATGAAGGCAATGCT 7 TTTGTCATTATGGGAG 23  25 CT 26 TT 26 CC 25 GGAGAACAATTAAAACGTATTAAATATGA 27 AGTTGGTGAAAATAACTTAAAGGTATTCAACGT 29 A 30 G 29 CACTTTAATAATAATCACGAGTTAGTTAGTTC 28 TGCTAATGAAAATAAATTAAAAGTATACAACGTACACTTTGATAATAATCAAGAGTTAGTTGCTGA 27 TGGTGAGCCTGACGTA 31  33 ATATGTTTAAGCAAGCAGGTCTGGGAAAATCTTCTCATTAAACTAAAGCTGGAAAA 35 C 36 T 35 AATGAAAATGTGTTTTCTGAAACTAAAAAATTATCGAATAAAAATAATG 37 CCGATCAGTTTTTTGAATGCGCTAA 39 AAGAAATGAA 40 GAGAAATGAA 39  38 ACGATCAGTTTTTTGAATGCGCTAAAAGAAATGAACAGAACCTTTTGATAATA 38 CCGATCAGTTTTTT 37  34 GTATTTTTAAAAAAACGGTGTGGGAAGACCTTCTCGTTAAAC 33  32 ATATGTTTAAGCAAGCAGGTCTGGGAAAATCTTCTCATTAAACTAAAGCTGGAAAACAATGAAAATGTGTTTTCTGAAACTAAAAAATTATCGAATAAAAATAATGACGATCAGTTTTTTGAATGCGCTAAAAGAAATGAACAGAACCTTTTCGATAATATAAGAAAAAGTGATTTTCATGTTGGTTTACTTAAGCCAAGTAGTACGCGTAGTGTTATTTTAGAAACGCCGCCAAATGTCTGTATGGAATCACGTAATTCATATGAAAAAAA 41 A 42 TAGA 41  31  24 CTGGGGAACAATTAAAACGTATTAAATATGATGTTGATGAAAATAACTTAAAGGTATTCAACGTACACTTTGATAATAATGAAGTGTTAGTTACTGATGGTGAGCCTGACGTAGTATGTTTAAGCAAGCAGGTCTGGGAAAATCTTCTCATTAAATTAAAACCGGAGATCAAGGAAAATGTGGCTTCTGAAGTTCATAAATCAGCGAATAAAGGTGAGATTGAGCAATTAGTTGAATGGTCTAAAAGAAATGAACAAACCCTTTTCGATAAT 43 A 44 G 43 TAATAAAAAGTGATTTTCATGTTGGTTCACTTAAGCCAGGTAGTATGAATGGTGTTATTTTAGAAATGCCACCAAATGTCTGTATGGAACCACGTAATTCATATGAAAACAAAATAGATGAGGTTTCATCTTTGTCAGAGTCAGAGGAACACCCCATAGATATTCAAAAAATAACAGATGCGTTTGTGAAGGAGTTCAAGGGGATATTATTTGATAAAAATGGAAG 45 A 46 G 45 TCTTCAGAGCTTCTGTTTAATTTTTATGAATGTTGCTATACGTTTTTACCAAGAGCGCAGCCTCAAGATAAAATCGATAGCTATAATTCAGCACTGCAAGCTTTTTCCATCTTTTGTTCATCTACGTTGACACATAATAATGTAGGCTTTGATTTCAAATTATTTCCAGAAGTCAAGCTGTCTGGAGAACATCTTGAAACGGTATTCAAATACAAAAATGGCGATGATGTCCGGGAGATAGCCAAAATTAACATTACTCTCCAAAAAGAAGAGGGTGGCTTATATAATCTACGTGGATTGGATTTTAAGGGATGCTTCTTTTCTGGACAGAACTTCAGTAACTATGATATTCAATATGTGAACTGGGGAATGTCATTGTTTGATGTTGATACTCCGTGTATTTTTAATACGCCTGCTAACCATGAGAGTTATGAAAAATCATTAAAACCTGTAAGCGAAAACGGTTTAAATGGAGTCTTATCTGATCGTAATAAAAAAATAAAAATGATCACGGGTGTGGCACCATTCGATGATATTTTATTTATGGATGATGACTTTGATGATAACTCCCCTGAGGATGCTCCCATTGAGAATAGTCCTGTTGTGAATAGTCCTCTTGTA 24 CTGGAAAACAA 23  6  47 ATGTTAGTTAGTAAAAGCAACGAAATTAACACTAGCACATTTATAAATAGTGGAA 49 G 50 A 49 TTGTAATGAAAGAAAATC 51 T 52 C 51 TCTAAATCCATGGAGTTACTAGCTTATAGTATTATAAAATTAATTTGTAAAGAAGCTGCCTCAGGGACGTATAACGGTGCTCTTGAAATTTTACAAAAAATGATGTCTGAATGTAAA 53 T 54 C 53 ATCATGAAGGCAATGC 55 T 56 C 55 TTTGTCATTATGGGAGCCGGAGAACAATTAAAACGTATTAAATATGATGCTAATGAAAATAAATTAAAAGTATACAACGTACACTTTGATAATAATCAAGAGTTAGTTGCTGATGGTGAGCCTGACGTAGTATTTTTAAAAAAAACGGTGTGGGAAGA 57 CCTTCTCG 58 TCTTCTCA 57  48 GTGTGGGAAGACCTTCTCG 47 TTAAACTAAAGCTGGAGAACAAAGAAAATGCGGTTTCTGAAA 59 TTAACC 60 CTGACA 59 TGTCATCTAATAAAAATAATGTTGATCAGTTTATTGAATGCGCTAAGAGAAATGAACAGACCCTATTCGGCAATATAAGAAAAAGTGATTTTCA 61 T 62 C 61 GTTGCTTCACTTCAGCCAGGTAGAACGCGTAGTGTTATTTCAGAAACGCCGCCAAATGACTGTATGGAATCACATAATTTATATGAAAACCACACAGATA 63 C 64 A 63 GGTTTCAACTGTAACTAAAAATTCTCAGCAAGTTAAAGGTCACTATGG 65 A 66 G 65 GATAAGTTGAAAGAAATGCAGTT 67 G 68 T 67 TTCCTCAACCAGATGAGCAATGCACTTCAACAGGATTCATCTTTGTTAGAGTCAAAGGAACACACCATAGATATTCAAGAAAAAACGAATAAGTTTGTGCAGCATTTTCAGCGGGTATTATTTGATAAAAATGGAAGGTCATCAGAGTTTCTACTTAATTTTTATGAGTGTTGCTATAAGTTTTTACCAAGAGCGCAGCCTCAAGATAAAATCGATAGCTATAATTCAGCACTGCAAGCTTTTTCCATCTTTTGTTCATCTACGTTGA 69 C 70 T 69 ACATAATAATGTAGGGTTTAATTTCAAATTATTTCCAGAAGTCAAG 71 C 72 T 71 TGTCTGGAGGAGAGCTTGAAACGGTATTCAAATACAAAAATGG 73 T 74 C 73 AATTTTGTCTGGGAGATAGCCAGAATTAAAATT 75 A 76 G 75 CTCTCCCAAAAGAAGAGGGT 77 G 78 A 77 GTTTATATAATTTACGTGGATTGGATTTTAAGGGATGCTTCTTTTCTGGACAGAACTTCAGTAACTATGATATTCAATATGTGAACTGGGGAACGTCATTGTTTGAT 79  81 CTTGATACTCCA 82 CTTGATACTCCG 81 TGTATTTTTAATGCGCCTGCTTACAACAAGAGTAATGAAAAATCATTA 83 G 84 A 83 AACGCGTCAGCGAAAACGGTTTAAGTGGAGTCTTGTCTGATCGTAATAA 85 A 86 T 85  80 GTTGATACTCCGTGTATTTTTAATGTGCCTGATGACAATAAGAGTTATGATAAATTATTAAAATCCGTCAGCGAAAATGGTTTAAATGGAGTCTTGACTGATCGTAATAAT 79 AAAATAAAACTAATCACGGG 87 T 88 C 87 GTGGCACCATTCGATG 89 ATATTTTATT 90 GTATTTCATC 89 TATGGATGATGACTTTGATGATA 91  93 ACTCCCCTGATGATGG 94 GTTCTTCTGAGGATGA 93 TCCCGTTGAGAATAGTCCTGTTGTGAATAGTCCCCTTGTA 92 GTTCCTCTGAGGAAAATTCCCCTGAGGATAGTCCCATTGAGAATCGTCCCCTTGTA 91  6  95  97 ATGGAATCACGTAATTCATATGAAAACAAAATAGATGAGATTTCATCTTTGTCAGAGTCAAAGGAACACCCCATAGATATTCAAGAAAAAAAAGATG 99 CGTTT 100 TGTTT 99  98  98 ATGGAATCACGTAATTCATATGAAAACAAAATAGATGAGATTTCATCTTTGTCAGAGTTAAAGGAACACCCCATAGATATTCAAGAAAAAAAAGATGCGTTT 98 ATGGAATCACGTAATTCATATAAAAACAAAATAGATGAGATTTCATCTTTGTCAGAGTCAAAGGAACACCCCATAGATATTCAAGAAAAAAAAGATGCGTTT 97 GTGAATGAGTTCAAGGGGGTATTATTTGATAA 101 AAATAC 102 GAATAC 102 AAATAT 101  96  103 ATGTTAGTTAGTAAAAACAACGAATTTAACACTAGCTCATTTATAGATAGTGGAAATTGTAATGAAAGAAAATCTTCTGAATCCATGGAGCTACTAGCTTATAGTATTATAAATTTAATTTGTAAGAAAGCTGCCTCAGGGACGTATCGCGGTGCTCTTGAAACTTTACAAAAAATG 104  103 ATGTCTGAATGTATATATCA 105 T 106 A 105 GAAGGCAACGCCTTTGTCATTATGGGAG 107 T 108 C 107 TGGAGAACAATTAAAACGTATTAAATATGAAGTTGGTGAAAATAACTTAAAGGTATTCAACGTACACTTTAATAATAATCACGAGTTAGTTAGTTCTGGTGAGCCTGACGTAATATGTTTAAGCAAGCAGGTCTGGGAAAATCTTCTCATTAAACTAAAGCTGGAAAACAATGAAAATGTGTTTTCTGAAACTAAAAAATTATCGAATAAAAATAATG 109 A 110 C 109 CGATCAGTTTTTTGAATGCGCTAAAAGAAATGAACAGAACCTTTTCGATAATATAAGAAAAAGTGATTTTCATGTTGGTTTACTTAAGCCAAGTAGTACGCGTAGTGTTATTTTAGAAACGCCGCCAAATGTCTGTATGGAATCACGTAATTCATATGAAAACAAAATAGATGAGATTTCATCTTTGTCAGAGTCAAAGGAACACCCCATAGATATTCAAGAAAAAAAAGATGCGTTTGTGAATGAGTTCAAGGGGGTATTATTTGATAAAAATAC 96 ATGTTAGTTAGTAAAAGCAACGAAATTAACACTAGCACATTTATAAATAGTGGAAGTTGTAATGAAAGAAAATCTTCTAAATCCATGGAGTTACTAGCTTATAGTATTATAAAATTAATTTGTAAGGAAGCTGCCTCAGGGACGTATAACGGTGCTCTTGAAATTTTACAAAAAATGATGTCTGAATGTAAATATCATGAAAGCAATGCTTTTGTCATTATGGGAGCCGGAGAACAATTAAAACGTATTAAATATGATGTTAATGAAGATAAATTAAAAGTATACAACGTACACTTTGATAATAATCAAGAGTTAGTTGCTGATGGTGAGCCTGACGTAGTATTTTTAAAAAAAACGGTGTGGGAAGACCTTCTCGTTAAACTAAAGCTGGAGAACAAAGAAAATGCGGTTTCTGAAGTTCATAAATCAGCGAATAAAGGTGAGGTCGAGCAATTAGTTGAATGTTCTGAAAGAAATGAAAAGAGGCTTCTTGATAATATAAATAAAAACACTACTATCTATAATATTTACAACAACCAACGAGCAACAAACATTACTACATCTGTATCTGAACCACCTGCACAAGGTAAAAGTAATGATGTAGATAAATTAAAACAAACGCAGTTTATCAACGCTAGTCAATTTGATATTAATGAGCTTCAGCAGAGTTCTGCTTTTTTAGCGTCAAAATATCACACCATAGATATTCAAGAAAAAACAGATGCGTTTGTGAAGAAGTTCAAGGGGATATTATTTGATAAAAATGG 95 AAGGTCTTCAGAGCTTCT 111 G 112 T 111 TTTAATTTTTATGAGTGTTGCTATAAGTTTTTACCAAGAGC 113  115 T 116 G 115 CAGCCTCAAGATAAAAT 117 CGAT 118 TGAA 117  114 TCAGCCACAAGATAAAATCGAT 113 AGCTATAATTCAGCACTGCAAGCTTTTTCCAT 119 TTTTTGTTCATCTACGTTGACACATAATAATATAGGCTTTGATTTCAAATTATTTCCT 120 CTTTCGTTCATCTACTTTGAATAATAATGATGTAGGGTTTAATTTCAAATTATTCCCA 119 GAAGTCAAGCTGTCTGGA 121 G 122 A 121 AACATCTTGAAACGGTATTCAAATACAAAAATGGCGATGATGTCCGGGAGATAGCCAAAATTAACATTACTCTCCAAAAAGAAGAGGGTGGTTTATATAATTTACGTGGATTGGATTTTAAGGG 123 GT 124 AT 124 GG 123 GCTTCTTTTCTGGACAGAACTTCAGTA 125 A 126 T 125 CTATGATATTCAATATGTGAACTGGGGAACGTCATTGTTTGAT 127 GTTGATAC 128 CTTGATAC 128 GTTGATAT 127 TCCGTGTATTTTTAATGCGCCTGCTTACAACAAGAGTAATGAAAAATCATTA 129 AAACCTGTG 130 GAACGCGTC 129 AGCGAAAACGGTTTAAGTGGAGT 131 ATTGA 132 CTTGT 132 CTTGA 131 CTGATCGTAATAATAAAATAAAACTCATCACGGGCGTGGCACCATTCGATGATATTTTATTTATGGATGATGACTTTGATGATAGTTC 133 CTCTGAGGATG 134 CTCTGAGGATA 134 TTCTGAGGATG 133 ATCCCGTTGAGAATAGTCCTGTTGTGACTAGTCCC 135 G 136 C 135 TTGTATCAAGTTCTAAAAGCAGTTTTCAA 6  137  139  140 ATGTTAGTTAGTAAAAGCAACGAACTTAACACTAGCGCATTTTTTGCTAGCAGGAATTATAATGGAAACAATTCTTCCAACCCCATGGAGCTACTAGCTCATAGCATTATAAAATTAATTTGTAAGGAAGCTGCCTCAGCGACGTATTGCGGTGCTCTTGAAACTTTACAAAAAATGATGTCTGAATGTATATATCAAGAAGGCAACGCCTTTGTCATTATGGGAGCTGGAGAACAATTAAAACGTATTAAATATGATGTTGATGAAAATAACTTAAAGGTATTCAACGTACACTTTGATAATAATGAAGTGTTAGTTACTGATGGTGAGCCTGACGTAGTATGTTTAAGCAAGCAGGTCTGGGAAAATCTTCTCATTAAATTAAAACCGGAGATTAAAGAAAATGTGTTTTCTGAAAATAACAAATTATCGAGAGAAAATAATGTTGATCAGTATGTTCAATCCGCTAAAAGAAATGAACAGGCCCTTTTCGATAATATAATAAAAAGTGATTTTCATGTTGCTTCACTTCAGCCAGGTAGAACACGTAGTGTTATTTCAGAAACGCCGCCAAATGACTGTATGGAAGCACGTAATTCATATGAAAACCAAATAGATGAGATTTCATCT 139 TTGTCAGAGTCAAAGGAACACCCCATAGATATTCAAGAAAAAA 141 A 142 C 141 AGATGCGTTTGTGAATGAGTTCA 143 A 144 G 143  138 ATGTTAGTTAGTAAAAG 145 C 146 A 145 AACGGATTTAACGCTAGCGCAGTTTTGGGTAGTGGAAGTTATAATGAAAATAAATCTTCTAAACACATGGAGCTACTAGCTCATAGTATTTTAAAATTAATTTGTAAG 147 G 148 A 147 AAGCTGCATCAGAGACGTATCGCGGTGCTCTTGAAACTTTACAAAAAATGATGTCTGAATGT 149 A 150 G 149 TATATCAAGAAGGCAACGCCTTTGTCATTATGGGAGCTGGAGAACAATTAAAACGTATTAAATATGAAGTTGGTGAAAATAACTTAAAGGTATTCAACGTACACTTTAATAATAATCACGAGTTAGTTAGTTCTGGTGAGCCTGACGTAATATGTTTAAGCAAGCAGGTCTGGGAAAATCTTCTCATTAAACTAAAGCTGGAAAACAATGAAAATGTGTTTTCTGAAACTAAAAAA 151 T 152 C 151 TATCGAATAAAAATAATGCCGATCAGTTTTTTGAATGCGCTAA 153 GAGAAATGAACAGACCCTA 154 AAGAAATGAACAGAACCTT 153 TTCGATAATATAAGAAAAAGTGATTTTCATGTTGGTT 155 CACTTAAGCCAG 156 TACTTAAGCCAA 155 GTAGTACGCGTAGTGTTATTTTAGAAACGCCGCCAAATGTCTGTATGGAATCACGTAATTCATATGAAAACAAAATAGATGAGATTTCATCTTTGTCAGAGTCAAAGGAACACCCCATAGATATTCAAGAAAAAA 157 CAGATGCATTTGTGAAGAAGTTCAA 158 AAGATGCGTTTGTGAATGAGTTCAA 157  138 GTGAATGAGTTCAA 138 ATGTTAGTTAGTAAAAGCAACGGATTTAACGCTAGCGCAGTTTTGGGTAGTGGAAGTTATAATGAAAATAAATCTTCTAAACACATGGAGCTACTAGCTCATAGTATTGTAAAATTAATTTGTAAGGAAGCTGCATCAGAGACGTATCGCGGTGCTCTTGAAATTTTACAAAAAATAATGTCTGAATGTATATATCATGAAGGCAACGCCTTTGTCATTATGGGAGCTGGAGAACAATTAAAACGTATTAAATATGATGTTGATGAAAATAACTTAAAGGTATTCAACGTACACTTTGATAATAATGAAGTGTTAGTTACTGATGGTGAGCCTGACGTAGTATGTTTAAGCAAGCAGGTCTGGGAAAATCTTCTCATTAAATTAAAACCGGAGATCAAGGAAAATGCGGCTTCTGAAGTTCATAAATCAGCGAATAAAGGTGAGATTGAGCAATTAGTTGAATGCTCTAAAAGAAATGAACAGACCCTTTTCGATAATATAAGAAAAAGTGATTTTCATGTTGGTTCACTTAAGCCAGGTAGTATGAATAGTGTTATTTTAGAAATGCCACCAAATGTCTGTATGGAACCACGTAATCCATATGAAAACAAAATAGATGAGGTTTCATCTTTGTCAGAGCCAAAGGAACACACCATAGATATTCAAGAAAAAACAGATGCGTTTGTGAATGAGTTCAA 137 GGGGATATTATTTGAT 159 AAAAATGG 160 AAAAATAC 160 CAAAATGG 160 AAAAATAG 159 AAGGTCTTCAGAGTTTCT 161  163 AC 164 GT 163 TTAATTTTTATGAATGTTGCTAT 165 GA 166 GT 166 AA 165  162 ATTTAATTTTTATGAGTGTTGCTATAA 161 GTTTTTACCAAGAGC 167 GCAGCCTCAGGATAAAAT 169 CGAA 170 TGAA 169  168 TCAGCCTCAAGATAAAATCGAT 167 AGCTATAATTCAGCACTGCAAGCTTT 171 CTCCATC 172 TTCCATT 172 TTCCATC 171 TTTTGTTCATCTACGTTGA 173  175 T 176 C 175 ACATAATAATATAGGCTTTG 174 CACATAATGGTGTAGGGTTTA 173 ATTTCAAATTATTTCC 177 AGAAGTCAAACTGTGTGGGG 179 AAA 180 GAA 180 AAC 179  178  181 AGAAGTCAAGCTGTCTGGAGAAC 182 TGAAGTCAAACTGTCTGGAGAAC 181  177 ATCTTGAAACGGTATTCAAATA 183 T 184 C 183 AAAAATGGCGATGATGTCCGGGAGATAGCCAAAATTAACATT 185  187 A 188 G 187 CTCTCCAAAAAGAAGAGGATGGTTTATATAATTTAG 186 ATCCTACCAAAAGGAGAGGGTGATTTATATAATTTGG 186 ACTCTCCAAAAAGAAGA 189 GGGTGGATTATATAATCTAC 190  191 A 192 G 191 GGTGGTTTATATAATTTAC 189  185 GTGGATTGGATTTTAAGGG 193 A 194 G 193 TGCTTCTTTTCTGGACAGAACTTCAGTAACTATGATATTCAATATGTGAACTGGGG 195 A 196 T 195 ACGTCATTGTTTGAT 197  199 CTTGATACTCCG 200 GTTGATACTCCG 200 CTTGATACTCCA 199 TGTATTTTTAATGCGCCTGCTTACAACAAGAGTAATG 201  203 A 204 G 203 AAAATCATTAAAACCTGTG 202 AAAAATCATTAGAACGCGTC 201  198 GTTGATACTCCGTGTATTTTTAATACGCCTGCTAACCATGAGAGTTATGAAAAATCATTAAAACCTGTA 197 AGCGAAAACGGTTTAA 205  207 GTGGAGTC 208 GTGGAGTA 207 TTGACTGATCGTAATAATAAAATAAAACTCATCACGGGC 206  209 ATGGAGTCTTA 210 GTGGAGTCTTG 209 TCTGATCGTAATAAAAAAATAAAA 211 ATGATCACGGGT 212 CTAATCACGGGT 211  205 GTGGCACCATTCGATGATATTTTATTTATGGATGATGACTTTGATGATA 213 GTTCCTCTGAGGATGATCCCG 214 ACTCCCCTGAGGATG 215 CTCCCA 216 GTCCCG 215  213 TTGAGAATAGTCCTGTTGTGA 217 CTAGTCCCGTTGTATCAAGTTCTAAAAGCAGTTTTCAA 218 ATAGTCCCCTTGTA 217  6 ATGTTAGTTAGTAAAAG 219 TAATGGCCATAACGTTAATACAGTTTTGGGTAGTAGAAGTTGTAATGAAAACAATTCTTCCAACCCCATGGG 220 CAACGGATTTAACGCTAGCGCAGTTTTGGGTAGTGGAAGTTATAATGAAAATAAATCTTCTAAACACATGGA 219 GCTACTAGCTCATAGTATT 221 G 222 T 221 TAAAATTAATTTGTAAGGAAGCTGCATCAGAGACGTATCGCGGTGCTCTTGAAACTTTACAAAAAATGATGTCTGAATGTATATATCA 223 T 224 A 223 GAAGGCAACGCCTTTGTCATTATGGGA 225 TCTGGAGTG 226 GCTGGAGAA 226 TCTGGAGAG 225 CAATTAAAACGTATTAAATATGATGTTGATGAAAATAACTTAAAAGTATTCAACGTATACTTTGATAATAATGAAGAGTTAGTTACTGATGGTGAGCCTGACGTAGTATGTTTAAGCAAGCAGGTCTGGGAAGATCTTCTCATTAAACTAAAGCTGGAGAACAAAGAAAATGCGGTTTCTGAAACTGACCTGTCATCGAATGAAAATAATGTTGATCATTTTTTTCAATCCGCTAAGAGAGATGAACAGACCCTATTCGGCAATATAAGAAAAAGTGAGTTTCATGTTGATTCATTAAAGCCAGGTAGTACGCGTAGTGTTATTTTAGAAACGCAACCAAATGTCTCTATGGAACCACATAATTTATATGACAACCAAATAGATAAGGTTTCACCTGTAACTAAAAACTCTCAGCAAGTTAAAGGCCACTATGGAGATAAATTGAAAGAAATGCAGATGTTCCTTAACCAGATGAGCAATGCACTTCAACAGGCTCCATCTTTGTCAGAGTCAAAGGAACACACCATAGATATTCAAAAAATAACGGATGCCTTTGTGAAGGAGTTCCGGGGGATATTATTTGATAAAAATGGAAATTCATCAGAACGTTTATTCAATTTTTATGAGTGCTGCTATATATTTTTACCAAGAGCGCAGCCTCAAGATAAAATCGAAAGCTATAATTCAGCGCTGCAAGCTTTTTCCATCTTCCGTTCATCTACTTTGAATAATAATGATGTAGGGTTTAATTTCAAATTATTCCCAGAAGTCAAGCTGTCTGGTGAAAATCTTGAAACGGTATTCAAATACAAAAAAGGAAGTTTTGTCCGGGAGATAGCCAGAATTAACATTACTCTCCAAAAAGAAGAGGATGGTTTATATAATTTAGGTGGATTGGATTTTAAGGGATGCTTCTTTTCTGGACAGAACTTCAGTAACTATGATATTCAATATGTGAACTGGGGAACGTCATTGTTTGATCTTGATACTCCGTGTATTTTTAATGCGCCTGCTTACAACAAGAGTAATGAAAAATCATTAAAACCTGTGAGCGAAAACGGTTTAAGTGGAGTCTTGACTGATCGTAATAATAAAATAAAACTCATCACGGGCGTGGCACCATTCGATGATATTTTATTTATGGATGATGACTTTGATGATAGTTCCTCTGAGGATGATCCCGTTGAGAATAGTCCTGTTGTGACTAGTCCCGTTGTATCAAGTTCTAAAAGCAGTTTTCAA 5 ");
    string bad = "ATGTTAGTTAGTAAAAGCAACGGATTTAACGCTAGCGCAGTTTTGGGTAGTGGAAGTTATAATGAAAATAAATCTTCTAAACACATGGAGCTACTAGCTCATAGTATTTTAAAATTAATTTGTAAGGAAGCTGCATCAGAGACGTATCGCGGTGCTCTTGAAACTTTACAAAAAATGATGTCTGAATGTATATATCAAGAAGGCAACGCCTTTGTCATTATGGGAGCTGGAGAACAATTAAAACGTATTAAATATGAAGTTGGTGAAAATAACTTAAAGGTATTCAACGTACACTTTAATAATAATCACGAGTTAGTTAGTTCTGGTGAGCCTGACGTAATATGTTTAAGCAAGCAGGTCTGGGAAAATCTTCTCATTAAACTAAAGCTGGAAAACAATGAAAATGTGTTTTCTGAAACTAAAAAATTATCGAATAAAAATAATGCCGATCAGTTTTTTGAATGCGCTAAAAGAAATGAA";
    EXPECT_TRUE(l.get_valid_vcf_reference(bad).empty());
}

TEST(LocalPRGTest, get_valid_vcf_reference_valid_simple) {
    LocalPRG test_prg(3, "long_enough", "AGTATA 5 GCC 7 CCC 8 TATG 7  6 GGAGCG 5 TATTTACGTTCGAGGTCCAGACGCTCTA");
    std::string valid1 = "AGTATAGCCCCCTATTTACGTTCGAGGTCCAGACGCTCTA";
    std::string valid2 = "AGTATAGCCTATGTATTTACGTTCGAGGTCCAGACGCTCTA";
    std::string valid3 = "AGTATAGGAGCGTATTTACGTTCGAGGTCCAGACGCTCTA";
    std::vector<LocalNodePtr> exp_nodes1 = {test_prg.prg.nodes[0], test_prg.prg.nodes[1], test_prg.prg.nodes[2],
                                            test_prg.prg.nodes[4], test_prg.prg.nodes[6]};
    std::vector<LocalNodePtr> exp_nodes2 = {test_prg.prg.nodes[0], test_prg.prg.nodes[1], test_prg.prg.nodes[3],
                                            test_prg.prg.nodes[4], test_prg.prg.nodes[6]};
    std::vector<LocalNodePtr> exp_nodes3 = {test_prg.prg.nodes[0], test_prg.prg.nodes[5], test_prg.prg.nodes[6]};

    auto result1 = test_prg.get_valid_vcf_reference(valid1);
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, exp_nodes1, result1);
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, exp_nodes2, test_prg.get_valid_vcf_reference(valid2));
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, exp_nodes3, test_prg.get_valid_vcf_reference(valid3));
}

TEST(LocalPRGTest, get_valid_vcf_reference_valid_rev) {
    LocalPRG test_prg(3, "long_enough", "AGTATA 5 GCC 7 CCC 8 TATG 7  6 GGAGCG 5 TATTTACGTTCGAGGTCCAGACGCTCTA");
    std::string valid1 = rev_complement("AGTATAGCCCCCTATTTACGTTCGAGGTCCAGACGCTCTA");
    std::string valid2 = rev_complement("AGTATAGCCTATGTATTTACGTTCGAGGTCCAGACGCTCTA");
    std::string valid3 = rev_complement("AGTATAGGAGCGTATTTACGTTCGAGGTCCAGACGCTCTA");
    std::vector<LocalNodePtr> exp_nodes1 = {test_prg.prg.nodes[0], test_prg.prg.nodes[1], test_prg.prg.nodes[2],
                                            test_prg.prg.nodes[4], test_prg.prg.nodes[6]};
    std::vector<LocalNodePtr> exp_nodes2 = {test_prg.prg.nodes[0], test_prg.prg.nodes[1], test_prg.prg.nodes[3],
                                            test_prg.prg.nodes[4], test_prg.prg.nodes[6]};
    std::vector<LocalNodePtr> exp_nodes3 = {test_prg.prg.nodes[0], test_prg.prg.nodes[5], test_prg.prg.nodes[6]};

    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, exp_nodes1, test_prg.get_valid_vcf_reference(valid1));
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, exp_nodes2, test_prg.get_valid_vcf_reference(valid2));
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, exp_nodes3, test_prg.get_valid_vcf_reference(valid3));
}

TEST(LocalPRGTest, get_valid_vcf_reference_invalid) {
    LocalPRG test_prg(3, "long_enough", "AGTATA 5 GCC 7 CCC 8 TATG 7  6 GGAGCGTCGAGGTCCAGTCGAGGTCCAG 6  5 TATTTACGTTCGAGGTCCAGACG");
    std::string null = "";
    std::string snp_away_from_graph = "AGTATAGCCCCCTAGTTACGTTCGAGGTCCAGACG";
    std::string too_short = "AGTATATATTTACGTTCGAGGTCCAGACG";
    std::string starts_late = "TATAGCCCCCTATTTACGTTCGAGGTCCAGACG";
    std::string ends_less_than_a_node_early = "AGTATAGCCCCCTATTTACGTTCGAGGTCCAGAC";
    std::string ends_a_node_early = "AGTATAGGAGCGTCGAGGTCCAGTCGAGGTCCAG";

    EXPECT_TRUE(test_prg.get_valid_vcf_reference(null).empty());
    EXPECT_TRUE(test_prg.get_valid_vcf_reference(snp_away_from_graph).empty());
    EXPECT_TRUE(test_prg.get_valid_vcf_reference(too_short).empty());
    EXPECT_TRUE(test_prg.get_valid_vcf_reference(starts_late).empty());
    std::vector<LocalNodePtr> exp_nodes1 = {test_prg.prg.nodes[0], test_prg.prg.nodes[1], test_prg.prg.nodes[2],
                                            test_prg.prg.nodes[4], test_prg.prg.nodes[7]};
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, exp_nodes1, test_prg.get_valid_vcf_reference(ends_less_than_a_node_early));
    EXPECT_TRUE(test_prg.get_valid_vcf_reference(ends_a_node_early).empty());

}

TEST(LocalPRGTest, random_path) {
    LocalPRG test_prg(3, "long_enough", "AGTATA 5 GCC 7 CCC 8 TATG 7  6 GGACCAG 6  5 TATTTACG");
    std::set<std::string> random_paths;
    while (random_paths.size() < 4) {
        random_paths.insert(test_prg.random_path());
    }
    std::set<std::string> exp_random_paths = {"AGTATAGCCCCCTATTTACG",
                                              "AGTATAGCCTATGTATTTACG",
                                              "AGTATAGGACCAGTATTTACG",
                                              "AGTATATATTTACG"};
    EXPECT_ITERABLE_EQ(std::set<std::string>, exp_random_paths, random_paths);
}