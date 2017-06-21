#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "localPRG.h"
#include "minimizer.h"
#include "minirecord.h"
#include "minihit.h"
#include "interval.h"
#include "path.h"
#include "localgraph.h"
#include "localnode.h"
#include "index.h"
#include "inthash.h"
#include "pannode.h"
#include "utils.h"
#include "seq.h"
#include "kmergraph.h"
#include "kmernode.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class LocalPRGTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(LocalPRGTest, create){
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

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

TEST_F(LocalPRGTest, isalphaEmptyString){
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    bool a0 = l0.isalpha_string("");
    bool a1 = l1.isalpha_string("");
    bool a2 = l2.isalpha_string("");
    bool a3 = l3.isalpha_string("");
    EXPECT_EQ(a0, 1) << "isalpha_string thinks the empty string is not alphabetic for l0";
    EXPECT_EQ(a1, 1) << "isalpha_string thinks the empty string is not alphabetic for l1";
    EXPECT_EQ(a2, 1) << "isalpha_string thinks the empty string is not alphabetic for l2";
    EXPECT_EQ(a3, 1) << "isalpha_string thinks the empty string is not alphabetic for l3";
}

TEST_F(LocalPRGTest, isalphaSpaceString){
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    bool a0 = l0.isalpha_string("AGCT T");
    bool a1 = l1.isalpha_string("AGCT T");
    bool a2 = l2.isalpha_string("AGCT T");
    bool a3 = l3.isalpha_string("AGCT T");         
    EXPECT_EQ(a0, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l0";
    EXPECT_EQ(a1, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l1";
    EXPECT_EQ(a2, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l2";
    EXPECT_EQ(a3, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l3";
}

TEST_F(LocalPRGTest, isalphaNumberString){
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    bool a0 = l0.isalpha_string("AGCT 8 T");
    bool a1 = l1.isalpha_string("AGCT 8 T");
    bool a2 = l2.isalpha_string("AGCT 8 T");
    bool a3 = l3.isalpha_string("AGCT 8 T");
    EXPECT_EQ(a0, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l0";
    EXPECT_EQ(a1, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l1";
    EXPECT_EQ(a2, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l2";
    EXPECT_EQ(a3, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l3";
}

TEST_F(LocalPRGTest, stringAlongPath){
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    // empty interval
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    EXPECT_EQ("", l0.string_along_path(p));
    EXPECT_EQ("", l1.string_along_path(p));
    EXPECT_EQ("", l2.string_along_path(p));
    EXPECT_EQ("", l3.string_along_path(p));

    // positive length interval
    d = {Interval(1,3)};
    p.initialize(d);
    EXPECT_EQ("GC", l1.string_along_path(p));
    EXPECT_EQ(" 5", l2.string_along_path(p));
    EXPECT_EQ(" 5", l3.string_along_path(p));

    // multiple intervals
    d = {Interval(0,1), Interval(2,3)};
    p.initialize(d);
    EXPECT_EQ("AC", l1.string_along_path(p));
    EXPECT_EQ("A5", l2.string_along_path(p));
    EXPECT_EQ("A5", l3.string_along_path(p));
    
    // including empty interval
    d = {Interval(0,1), Interval(2,2)};
    p.initialize(d);
    EXPECT_EQ("A", l1.string_along_path(p));
    EXPECT_EQ("A", l2.string_along_path(p));
    EXPECT_EQ("A", l3.string_along_path(p));

    // forbidden paths
    d = {Interval(2,3), Interval(13,25)};
    p.initialize(d);
    EXPECT_DEATH(l1.string_along_path(p),"");
    EXPECT_DEATH(l1.string_along_path(p),"");
    EXPECT_DEATH(l2.string_along_path(p),"");
    EXPECT_DEATH(l3.string_along_path(p),"");
}

TEST_F(LocalPRGTest, nodesAlongPath)
{
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");
    
    // empty interval expects no nodes along
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    vector<LocalNode*> v;

    //EXPECT_EQ(v, l0.nodes_along_path(p));
    EXPECT_EQ(v, l1.nodes_along_path(p));
    EXPECT_EQ(v, l2.nodes_along_path(p));
    EXPECT_EQ(v, l3.nodes_along_path(p));

    // positive length interval
    d = {Interval(1,3)};
    p.initialize(d);
    uint32_t j = 1;
    EXPECT_EQ(j, l1.nodes_along_path(p).size());
    j = 0;
    EXPECT_EQ(j, l1.nodes_along_path(p)[0]->id);
    EXPECT_EQ(j, l2.nodes_along_path(p).size()); // no nodes in this interval
    EXPECT_EQ(j, l3.nodes_along_path(p).size());
    // different interval
    d = {Interval(4,5)};
    p.initialize(d);
    j = 1;
    EXPECT_EQ(j, l2.nodes_along_path(p).size());
    EXPECT_EQ(j, l3.nodes_along_path(p).size());
    EXPECT_EQ(j, l2.nodes_along_path(p)[0]->id);
    EXPECT_EQ(j, l3.nodes_along_path(p)[0]->id);

    // multiple intervals
    d = {Interval(0,1), Interval(4,5)};
    p.initialize(d);
    j = 1;
    EXPECT_EQ(j, l1.nodes_along_path(p).size());
    j = 2;
    EXPECT_EQ(j, l2.nodes_along_path(p).size());
    EXPECT_EQ(j, l3.nodes_along_path(p).size());
    j = 0;
    EXPECT_EQ(j, l1.nodes_along_path(p)[0]->id);
    EXPECT_EQ(j, l2.nodes_along_path(p)[0]->id);
    EXPECT_EQ(j, l3.nodes_along_path(p)[0]->id);
    j = 1;
    EXPECT_EQ(j, l2.nodes_along_path(p)[1]->id);
    EXPECT_EQ(j, l3.nodes_along_path(p)[1]->id);

    // including empty interval
    d = {Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    j = 3;
    vector<LocalNode*> w = l3.nodes_along_path(p);
    EXPECT_EQ(j, w.size());
    EXPECT_EQ(j, l3.nodes_along_path(p)[0]->id);
    j = 4;
    EXPECT_EQ(j, l3.nodes_along_path(p)[1]->id);
    j = 6;
    EXPECT_EQ(j, l3.nodes_along_path(p)[2]->id);

    // a path with an empty node at end
    d = {Interval(12,13), Interval(16,16), Interval(23,23)};
    p.initialize(d);
    j = 3;
    w = l3.nodes_along_path(p);
    EXPECT_EQ(j, w.size());
    EXPECT_EQ(j, l3.nodes_along_path(p)[0]->id);
    j = 4;
    EXPECT_EQ(j, l3.nodes_along_path(p)[1]->id);
    j = 6;
    EXPECT_EQ(j, l3.nodes_along_path(p)[2]->id);

    // and a path which ends on a null node
    d = {Interval(12,13), Interval(16,16)};
    p.initialize(d);
    j = 2;
    w = l3.nodes_along_path(p);
    EXPECT_EQ(j, w.size());
    j = 3;
    EXPECT_EQ(j, l3.nodes_along_path(p)[0]->id);
    j = 4;
    EXPECT_EQ(j, l3.nodes_along_path(p)[1]->id);

    // and a path that can't really exist still works
    d = {Interval(12,13),Interval(19,20)};
    p.initialize(d);
    j = 2;
    EXPECT_EQ(j, l3.nodes_along_path(p).size());
    j = 3;
    EXPECT_EQ(j, l3.nodes_along_path(p)[0]->id);
    j = 5;
    EXPECT_EQ(j, l3.nodes_along_path(p)[1]->id);
}

TEST_F(LocalPRGTest, split_by_siteNoSites){
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    vector<Interval> v0, v1;
    v0.push_back(Interval(0,0));
    EXPECT_ITERABLE_EQ( vector< Interval >, v0, l0.split_by_site(Interval(0,0)));// << "Failed to split empty string with input Interval";
    v1.push_back(Interval(0,4));
    EXPECT_ITERABLE_EQ( vector< Interval >, v1, l1.split_by_site(Interval(0,4)));// << "Failed to split string with input Interval";
    v1.clear();
    v1.push_back(Interval(0,2));
    EXPECT_ITERABLE_EQ( vector< Interval >, v1, l1.split_by_site(Interval(0,2)));// << "Failed to split string with short input Interval";
    v1.clear();
    v1.push_back(Interval(1,3));
    EXPECT_ITERABLE_EQ( vector< Interval >, v1, l1.split_by_site(Interval(1,3)));// << "Failed to split string with middle input Interval";
}

TEST_F(LocalPRGTest, split_by_siteSite){
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");

    vector<Interval> v2;
    v2.push_back(Interval(0,1));
    l2.next_site = 5;
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,1)));// << "Failed to split string in Interval" << Interval(0,1);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,2)));// << "Failed to split string in Interval" << Interval(0,2);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,3)));// << "Failed to split string in Interval" << Interval(0,3);
    //EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,4)));// << "Failed to split string in Interval" << Interval(0,4);
    v2.push_back(Interval(4,6));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,6)));// << "Failed to split string in Interval" << Interval(0,6);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,7)));// << "Failed to split string in Interval" << Interval(0,7);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,8)));// << "Failed to split string in Interval" << Interval(0,8);
    v2.push_back(Interval(9,10));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,10)));// << "Failed to split string in Interval" << Interval(0,10);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,11)));// << "Failed to split string in Interval" << Interval(0,11);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,12)));// << "Failed to split string in Interval" << Interval(0,12);
    v2.push_back(Interval(13,14));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(0,14)));// << "Failed to split string in Interval" << Interval(0,14);
    v2.clear();
    v2.push_back(Interval(5,6));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.split_by_site(Interval(5,8)));// << "Failed to split string in mid Interval" << Interval(5,8);
}

TEST_F(LocalPRGTest, split_by_siteNestedSite){
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");
    LocalPRG l4(4,"nested varsite start immediately", " 5 G 7 C 8 T 7  6 G 5 ");

    vector<Interval> v3;
    v3.push_back(Interval(0,1));
    l3.next_site = 5;
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,1)));// << "Failed to split string in Interval" << Interval(0,1);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,2)));// << "Failed to split string in Interval" << Interval(0,2);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,3)));// << "Failed to split string in Interval" << Interval(0,3);
    //EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,4)));// << "Failed to split string in Interval" << Interval(0,4);
    v3.push_back(Interval(4,16));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,16)));// << "Failed to split string in Interval" << Interval(0,6);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,17)));// << "Failed to split string in Interval" << Interval(0,7);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,18)));// << "Failed to split string in Interval" << Interval(0,8);
    v3.push_back(Interval(19,20));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,20)));// << "Failed to split string in Interval" << Interval(0,10);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,21)));// << "Failed to split string in Interval" << Interval(0,11);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,22)));// << "Failed to split string in Interval" << Interval(0,12);
    v3.push_back(Interval(23,24));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(0,24)));// << "Failed to split string in Interval" << Interval(0,14);
    l3.next_site = 7;    
    v3.clear();
    v3.push_back(Interval(4,5));
    v3.push_back(Interval(8,9));
    v3.push_back(Interval(12,13));
    v3.push_back(Interval(16,16));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.split_by_site(Interval(4,16)));// << "Failed to split string in mid Interval" << Interval(5,8);

    vector<Interval> v4;
    v4.push_back(Interval(0,0));
    l4.next_site = 5;
    EXPECT_ITERABLE_EQ( vector< Interval >,v4, l4.split_by_site(Interval(0,1)));// << "Failed to split string in Interval" << Interval(0,1);
    EXPECT_ITERABLE_EQ( vector< Interval >,v4, l4.split_by_site(Interval(0,2)));// << "Failed to split string in Interval" << Interval(0,2);
    v4.push_back(Interval(3,15));
    EXPECT_ITERABLE_EQ( vector< Interval >,v4, l4.split_by_site(Interval(0,15)));// << "Failed to split string in Interval" << Interval(0,6);
    EXPECT_ITERABLE_EQ( vector< Interval >,v4, l4.split_by_site(Interval(0,16)));// << "Failed to split string in Interval" << Interval(0,7);
    EXPECT_ITERABLE_EQ( vector< Interval >,v4, l4.split_by_site(Interval(0,17)));// << "Failed to split string in Interval" << Interval(0,8);
    v4.push_back(Interval(18,19));
    EXPECT_ITERABLE_EQ( vector< Interval >,v4, l4.split_by_site(Interval(0,19)));// << "Failed to split string in Interval" << Interval(0,10);
    EXPECT_ITERABLE_EQ( vector< Interval >,v4, l4.split_by_site(Interval(0,20)));// << "Failed to split string in Interval" << Interval(0,11);
    l4.next_site = 7;
    v4.clear();
    v4.push_back(Interval(0,4));
    v4.push_back(Interval(7,8));
    v4.push_back(Interval(11,12));
    v4.push_back(Interval(15,22));
    EXPECT_ITERABLE_EQ( vector< Interval >,v4, l4.split_by_site(Interval(0,22)));// << "Failed to split string in mid Interval" << Interval(5,8);
    
}

TEST_F(LocalPRGTest, buildGraph)
{
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    LocalGraph lg0;
    lg0.add_node(0,"",Interval(0,0));
    EXPECT_EQ(lg0, l0.prg);

    LocalGraph lg1;
    lg1.add_node(0,"AGCT", Interval(0,4));
    EXPECT_EQ(lg1, l1.prg);

    LocalGraph lg2;
    lg2.add_node(0,"A", Interval(0,1));
    lg2.add_node(1,"GC", Interval(4,6));
    lg2.add_node(2,"G", Interval(9,10));
    lg2.add_node(3,"T", Interval(13,14));
    lg2.add_edge(0,1);
    lg2.add_edge(0,2);
    lg2.add_edge(1,3);
    lg2.add_edge(2,3);
    EXPECT_EQ(lg2, l2.prg);

    LocalGraph lg3;
    lg3.add_node(0,"A", Interval(0,1));
    lg3.add_node(1,"G", Interval(4,5));
    lg3.add_node(2,"C", Interval(8,9));
    lg3.add_node(3,"T", Interval(12,13));
    lg3.add_node(4,"", Interval(16,16));
    lg3.add_node(5,"G", Interval(19,20));
    lg3.add_node(6,"T", Interval(23,24));
    lg3.add_edge(0,1);
    lg3.add_edge(0,5);
    lg3.add_edge(1,2);
    lg3.add_edge(1,3);
    lg3.add_edge(2,4);
    lg3.add_edge(3,4);
    lg3.add_edge(4,6);
    lg3.add_edge(5,6);
    EXPECT_EQ(lg3, l3.prg);
}

TEST_F(LocalPRGTest, shift){
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "AT 5 G 7 C 8 T 7  6 G 5 T");
    //LocalPRG l4(4, "much more complex", "TCATTC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AGCTG");
    LocalPRG l5(5, "one with lots of null at start and end, and a long stretch in between", " 5  7  9  11 AGTTCTGAAACATTGCGCGTGAGATCTCTG 12 T 11  10 A 9  8 C 7  6 G 5 ");
    LocalPRG l6(6, "one representing a possible deletion at end", "GATCTCTAG 5 TTATG 6  5 ");

    deque<Interval> d = {Interval(0,3)};
    Path p, q;
    p.initialize(d);
    d = {Interval(1,4)};
    q.initialize(d);
    vector<Path> v_exp = {q};
    EXPECT_ITERABLE_EQ(vector<Path>, v_exp, l1.shift(p));
    v_exp.clear();
    EXPECT_ITERABLE_EQ(vector<Path>, v_exp, l1.shift(q)); // there are no shifts over end of prg

    d = {Interval(0,1), Interval(4,6)};
    p.initialize(d);
    d = {Interval(4,6), Interval(13,14)};
    q.initialize(d);
    v_exp = {q};
    EXPECT_ITERABLE_EQ(vector<Path>, v_exp, l2.shift(p));
    v_exp.clear();
    EXPECT_ITERABLE_EQ(vector<Path>, v_exp, l2.shift(q));

    v_exp.clear();
    d = {Interval(0,2)};
    p.initialize(d);
    d = {Interval(1,2), Interval(5,6)};
    q.initialize(d);
    v_exp.push_back(q);
    d = {Interval(1,2), Interval(20,21)};
    q.initialize(d);
    v_exp.push_back(q);
    EXPECT_ITERABLE_EQ(vector<Path>, v_exp, l3.shift(p));

    v_exp.clear();
    d = {Interval(1,2), Interval(5,6)};
    p.initialize(d);
    d = {Interval(5,6), Interval(9,10)};
    q.initialize(d);
    v_exp.push_back(q);
    d = {Interval(5,6), Interval(13,14)};
    q.initialize(d);
    v_exp.push_back(q);
    EXPECT_ITERABLE_EQ(vector<Path>, v_exp, l3.shift(p));

    v_exp.clear();
    d = {Interval(0, 0), Interval(3, 3), Interval(6, 6), Interval(9, 9), Interval(13, 18)};
    p.initialize(d);
    d = {Interval(14, 19)};
    q.initialize(d);
    v_exp.push_back(q);
    EXPECT_ITERABLE_EQ(vector<Path>, v_exp, l5.shift(p));

    v_exp.clear();
    d = {Interval(3, 8)};
    p.initialize(d);
    d = {Interval(4, 9), Interval(20,20), Interval(23,23)};
    q.initialize(d);
    v_exp.push_back(q);
    d = {Interval(4, 9)};
    q.initialize(d);
    v_exp.push_back(q);
    EXPECT_ITERABLE_EQ(vector<Path>, v_exp, l6.shift(p));

    v_exp.clear();
    d = {Interval(4, 9)};
    p.initialize(d);
    d = {Interval(5, 9), Interval(12,13)};
    q.initialize(d);
    v_exp.push_back(q);
    EXPECT_ITERABLE_EQ(vector<Path>, v_exp, l6.shift(p));
}

TEST_F(LocalPRGTest, minimizerSketch){
    // note this is a bad test
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");
    LocalPRG l4(4, "much more complex", "TCATTC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AGCTG");
    LocalPRG l5(5, "one with lots of null at start and end, and a long stretch in between", " 5  7  9  11 AGTTCTGAAACATTGCGCGTGAGATCTCTG 12 T 11  10 A 9  8 C 7  6 G 5 ");

    Index* idx;
    idx = new Index();

    KmerHash hash;

    l0.minimizer_sketch(idx, 1, 3);
    uint32_t j = 0;
    EXPECT_EQ(j, idx->minhash.size());

    l1.minimizer_sketch(idx, 2, 3);
    j = 1;
    EXPECT_EQ(j, idx->minhash.size());
    l1.minimizer_sketch(idx, 1, 3);
    EXPECT_EQ(j, idx->minhash.size());
    j = 2;
    pair<uint64_t,uint64_t> kh = hash.kmerhash("AGC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    
    idx->clear();
    l2.minimizer_sketch(idx, 2, 3);
    j = 1;
    EXPECT_EQ(j, idx->minhash.size());
    l2.minimizer_sketch(idx, 1, 3);
    j = 2;
    EXPECT_EQ(j, idx->minhash.size());
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    j = 1;
    kh = hash.kmerhash("AGT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());

    idx->clear();
    l3.minimizer_sketch(idx, 2, 3);
    j = 2;
    EXPECT_EQ(j, idx->minhash.size());
    l3.minimizer_sketch(idx, 1, 3);
    j = 3;
    EXPECT_EQ(j, idx->minhash.size());
    j = 2;
    kh = hash.kmerhash("AGC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size()); //AGC
    kh = hash.kmerhash("AGT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size()); //AGTx2
    j = 1;
    kh = hash.kmerhash("GTT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    
    idx->clear();
    l4.minimizer_sketch(idx, 1, 3);
    j = 16;
    EXPECT_EQ(j, idx->minhash.size());
    j = 5;
    kh = hash.kmerhash("TCA",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    j = 4;
    kh = hash.kmerhash("CTA",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    j = 3;
    kh = hash.kmerhash("ACT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size()); 
    kh = hash.kmerhash("CAA",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("AAG",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("TCT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("AGC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    j = 2;
    kh = hash.kmerhash("TTC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size()); 
    kh = hash.kmerhash("CAC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("CTC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    j = 1;
    kh = hash.kmerhash("CAT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("ATT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("GTC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("GTT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("TGT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("CTG",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());

    idx->clear();
    l4.minimizer_sketch(idx, 3, 3);
    j = 10;
    EXPECT_EQ(j, idx->minhash.size());
    j = 4;
    kh = hash.kmerhash("CTA",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    j = 3;
    kh = hash.kmerhash("CTT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    j = 2;
    kh = hash.kmerhash("CAC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    j = 1;
    kh = hash.kmerhash("ATT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("ACT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("TCA",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("AAC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("GTC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("GAG",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    kh = hash.kmerhash("CTG",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());

    idx->clear();
    l5.minimizer_sketch(idx, 4, 5);
    EXPECT_EQ((idx->minhash.size()>2), true);

    idx->clear();
    delete idx;
}

TEST_F(LocalPRGTest, minimizerSketchSameAsSeqw1){
    string st = "ATGGCAATCCGAATCTTCGCGATACTTTTCTCCATTTTTTCTCTTGCCACTTTCGCGCATGCGCAAGAAGGCACGCTAGAACGTTCTGACTGGAGGAAGTTTTTCAGCGAATTTCAAGCCAAAGGCACGATAGTTGTGGCAGACGAACGCCAAGCGGATCGTGCCATGTTGGTTTTTGATCCTGTGCGATCGAAGAAACGCTACTCGCCTGCATCGACATTCAAGATACCTCATACACTTTTTGCACTTGATGCAGGCGCTGTTCGTGATGAGTTCCAGATTTTTCGATGGGACGGCGTTAACAGGGGCTTTGCAGGCCACAATCAAGACCAAGATTTGCGATCAGCAATGCGGAATTCTACTGTTTGGGTGTATGAGCTATTTGCAAAGGAAATTGGTGATGACAAAGCTCGGCGCTATTTGAAGAAAATCGACTATGGCAACGCCGATCCTTCGACAAGTAATGGCGATTACTGTATAGAAGGCAGCCTTGCAATCTCGGCGCAGGAGCAAATTGCATTTCTCAGGAAGCTCTATCGTAACGAGCTGCCCTTTCGGGTAGAACATCAGCGCTTGGTCAAGGATCTCATGATTGTGGAAGCCGGTCGCAACTGGATACTGCGTGCAAAGACGGGCTGGGAAGGCCGTATGGGTTGGTGGGTAGGATGGGTTGAGTGGCCGACTGGCTCCGTATTCTTCGCACTGAATATTGATACGCCAAACAGAATGGATGATCTTTTCAAGAGGGAGGCAATCGTGCGGGCAATCCTT";

    Index* idx;
    idx = new Index();
    LocalPRG l(0,"prg", st);
    l.minimizer_sketch(idx, 1, 15);

    Seq s = Seq(0, "read", st, 1, 15);

    //cout << l.kmer_prg.nodes.size() << " " << s.sketch.size() << endl;
    EXPECT_EQ(l.kmer_prg.nodes.size(), s.sketch.size()+2);

    set<Minimizer*, MiniPos> sketch(s.sketch.begin(), s.sketch.end());
    vector<KmerNode*>::iterator lit = l.kmer_prg.nodes.begin();
    lit++;

    for (set<Minimizer*, MiniPos>::iterator sit = sketch.begin(); sit != sketch.end(); ++sit)
    {
        EXPECT_EQ((*sit)->pos, (*lit)->path.path[0]);
	++lit;
    }
}

TEST_F(LocalPRGTest, minimizerSketchSameAsSeqw5){
    string st = "ATGGCAATCCGAATCTTCGCGATACTTTTCTCCATTTTTTCTCTTGCCACTTTCGCGCATGCGCAAGAAGGCACGCTAGAACGTTCTGACTGGAGGAAGTTTTTCAGCGAATTTCAAGCCAAAGGCACGATAGTTGTGGCAGACGAACGCCAAGCGGATCGTGCCATGTTGGTTTTTGATCCTGTGCGATCGAAGAAACGCTACTCGCCTGCATCGACATTCAAGATACCTCATACACTTTTTGCACTTGATGCAGGCGCTGTTCGTGATGAGTTCCAGATTTTTCGATGGGACGGCGTTAACAGGGGCTTTGCAGGCCACAATCAAGACCAAGATTTGCGATCAGCAATGCGGAATTCTACTGTTTGGGTGTATGAGCTATTTGCAAAGGAAATTGGTGATGACAAAGCTCGGCGCTATTTGAAGAAAATCGACTATGGCAACGCCGATCCTTCGACAAGTAATGGCGATTACTGTATAGAAGGCAGCCTTGCAATCTCGGCGCAGGAGCAAATTGCATTTCTCAGGAAGCTCTATCGTAACGAGCTGCCCTTTCGGGTAGAACATCAGCGCTTGGTCAAGGATCTCATGATTGTGGAAGCCGGTCGCAACTGGATACTGCGTGCAAAGACGGGCTGGGAAGGCCGTATGGGTTGGTGGGTAGGATGGGTTGAGTGGCCGACTGGCTCCGTATTCTTCGCACTGAATATTGATACGCCAAACAGAATGGATGATCTTTTCAAGAGGGAGGCAATCGTGCGGGCAATCCTT";

    Index* idx;
    idx = new Index();
    LocalPRG l(0,"prg", st);
    l.minimizer_sketch(idx, 5, 15);

    Seq s = Seq(0, "read", st, 5, 15);

    //cout << l.kmer_prg.nodes.size() << " " << s.sketch.size() << endl;
    EXPECT_EQ(l.kmer_prg.nodes.size(), s.sketch.size()+2);

    set<Minimizer*, MiniPos> sketch(s.sketch.begin(), s.sketch.end());
    vector<KmerNode*>::iterator lit = l.kmer_prg.nodes.begin();
    lit++;

    for (set<Minimizer*, MiniPos>::iterator sit = sketch.begin(); sit != sketch.end(); ++sit)
    {
        EXPECT_EQ((*sit)->pos, (*lit)->path.path[0]);
        ++lit;
    }
}

TEST_F(LocalPRGTest, minimizerSketchSameAsSeqw10){
    string st = "ATGGCAATCCGAATCTTCGCGATACTTTTCTCCATTTTTTCTCTTGCCACTTTCGCGCATGCGCAAGAAGGCACGCTAGAACGTTCTGACTGGAGGAAGTTTTTCAGCGAATTTCAAGCCAAAGGCACGATAGTTGTGGCAGACGAACGCCAAGCGGATCGTGCCATGTTGGTTTTTGATCCTGTGCGATCGAAGAAACGCTACTCGCCTGCATCGACATTCAAGATACCTCATACACTTTTTGCACTTGATGCAGGCGCTGTTCGTGATGAGTTCCAGATTTTTCGATGGGACGGCGTTAACAGGGGCTTTGCAGGCCACAATCAAGACCAAGATTTGCGATCAGCAATGCGGAATTCTACTGTTTGGGTGTATGAGCTATTTGCAAAGGAAATTGGTGATGACAAAGCTCGGCGCTATTTGAAGAAAATCGACTATGGCAACGCCGATCCTTCGACAAGTAATGGCGATTACTGTATAGAAGGCAGCCTTGCAATCTCGGCGCAGGAGCAAATTGCATTTCTCAGGAAGCTCTATCGTAACGAGCTGCCCTTTCGGGTAGAACATCAGCGCTTGGTCAAGGATCTCATGATTGTGGAAGCCGGTCGCAACTGGATACTGCGTGCAAAGACGGGCTGGGAAGGCCGTATGGGTTGGTGGGTAGGATGGGTTGAGTGGCCGACTGGCTCCGTATTCTTCGCACTGAATATTGATACGCCAAACAGAATGGATGATCTTTTCAAGAGGGAGGCAATCGTGCGGGCAATCCTT";

    Index* idx;
    idx = new Index();
    LocalPRG l(0,"prg", st);
    l.minimizer_sketch(idx, 10, 15);

    Seq s = Seq(0, "read", st, 10, 15);

    //cout << l.kmer_prg.nodes.size() << " " << s.sketch.size() << endl;
    EXPECT_EQ(l.kmer_prg.nodes.size(), s.sketch.size()+2);

    set<Minimizer*, MiniPos> sketch(s.sketch.begin(), s.sketch.end());
    vector<KmerNode*>::iterator lit = l.kmer_prg.nodes.begin();
    lit++;

    for (set<Minimizer*, MiniPos>::iterator sit = sketch.begin(); sit != sketch.end(); ++sit)
    {
        EXPECT_EQ((*sit)->pos, (*lit)->path.path[0]);
        ++lit;
    }
}

TEST_F(LocalPRGTest, minimizerSketchSameAsSeqw15){
    string st = "ATGGCAATCCGAATCTTCGCGATACTTTTCTCCATTTTTTCTCTTGCCACTTTCGCGCATGCGCAAGAAGGCACGCTAGAACGTTCTGACTGGAGGAAGTTTTTCAGCGAATTTCAAGCCAAAGGCACGATAGTTGTGGCAGACGAACGCCAAGCGGATCGTGCCATGTTGGTTTTTGATCCTGTGCGATCGAAGAAACGCTACTCGCCTGCATCGACATTCAAGATACCTCATACACTTTTTGCACTTGATGCAGGCGCTGTTCGTGATGAGTTCCAGATTTTTCGATGGGACGGCGTTAACAGGGGCTTTGCAGGCCACAATCAAGACCAAGATTTGCGATCAGCAATGCGGAATTCTACTGTTTGGGTGTATGAGCTATTTGCAAAGGAAATTGGTGATGACAAAGCTCGGCGCTATTTGAAGAAAATCGACTATGGCAACGCCGATCCTTCGACAAGTAATGGCGATTACTGTATAGAAGGCAGCCTTGCAATCTCGGCGCAGGAGCAAATTGCATTTCTCAGGAAGCTCTATCGTAACGAGCTGCCCTTTCGGGTAGAACATCAGCGCTTGGTCAAGGATCTCATGATTGTGGAAGCCGGTCGCAACTGGATACTGCGTGCAAAGACGGGCTGGGAAGGCCGTATGGGTTGGTGGGTAGGATGGGTTGAGTGGCCGACTGGCTCCGTATTCTTCGCACTGAATATTGATACGCCAAACAGAATGGATGATCTTTTCAAGAGGGAGGCAATCGTGCGGGCAATCCTT";

    Index* idx;
    idx = new Index();
    LocalPRG l(0,"prg", st);
    l.minimizer_sketch(idx, 15, 15);

    Seq s = Seq(0, "read", st, 15, 15);

    //cout << l.kmer_prg.nodes.size() << " " << s.sketch.size() << endl;
    EXPECT_EQ(l.kmer_prg.nodes.size(), s.sketch.size()+2);

    set<Minimizer*, MiniPos> sketch(s.sketch.begin(), s.sketch.end());
    vector<KmerNode*>::iterator lit = l.kmer_prg.nodes.begin();
    lit++;

    for (set<Minimizer*, MiniPos>::iterator sit = sketch.begin(); sit != sketch.end(); ++sit)
    {
        EXPECT_EQ((*sit)->pos, (*lit)->path.path[0]);
        ++lit;
    }
}

TEST_F(LocalPRGTest, localnode_path_from_kmernode_path)
{   
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");
    LocalPRG l4(4, "much more complex", "TC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AG");
    //LocalPRG l5(5, "one with lots of null at start and end, and a long stretch in between", " 5  7  9  11 AGTTCTGAAACATTGCGCGTGAGATCTCTG 12 T 11  10 A 9  8 C 7  6 G 5 ");
    
    Index* idx;
    idx = new Index();
    
    KmerHash hash;
    
    l3.minimizer_sketch(idx, 2, 3);
    vector<KmerNode*> kmp = {l3.kmer_prg.nodes[2], l3.kmer_prg.nodes[4]};
    vector<LocalNode*> lmp = l3.localnode_path_from_kmernode_path(kmp);
    vector<LocalNode*> lmp_exp = {l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6]};
    EXPECT_ITERABLE_EQ( vector<LocalNode*>,lmp_exp, lmp);
    lmp = l3.localnode_path_from_kmernode_path(kmp, 2);
    lmp_exp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6]};
    EXPECT_ITERABLE_EQ( vector<LocalNode*>,lmp_exp, lmp);

    idx->clear();
    l4.minimizer_sketch(idx, 3, 3);
    kmp = {l4.kmer_prg.nodes[2], l4.kmer_prg.nodes[4], l4.kmer_prg.nodes[8]};
    lmp = l4.localnode_path_from_kmernode_path(kmp, 2);
    lmp_exp = {l4.prg.nodes[1], l4.prg.nodes[3]};
    EXPECT_ITERABLE_EQ( vector<LocalNode*>,lmp_exp, lmp);
    lmp = l4.localnode_path_from_kmernode_path(kmp, 3);
    lmp_exp = {l4.prg.nodes[0], l4.prg.nodes[1], l4.prg.nodes[3], l4.prg.nodes[4], l4.prg.nodes[6]};
    EXPECT_ITERABLE_EQ( vector<LocalNode*>,lmp_exp, lmp);

    /*idx->clear();
    l5.minimizer_sketch(idx, 4, 3);
    cout << l5.kmer_prg;*/

    delete idx;
}

TEST_F(LocalPRGTest, updateCovgWithHit)
{
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");

    Index* idx;
    idx = new Index();

    l3.minimizer_sketch(idx, 1, 3);

    KmerHash hash;
    Minimizer* m;
    pair<uint64_t,uint64_t> kh = hash.kmerhash("AGC", 3);
    m = new Minimizer(min(kh.first,kh.second), 1,4,1);
    deque<Interval> d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    Path p;
    p.initialize(d);
    MiniRecord* mr;
    mr = new MiniRecord(3,p,0);
    MinimizerHit* mh;
    mh = new MinimizerHit(1, m, mr);

    l3.update_covg_with_hit(mh);
    uint j = 1;
    EXPECT_EQ(j, l3.kmer_prg.nodes[2]->covg[0]);
    j = 0;
    EXPECT_EQ(j, l3.kmer_prg.nodes[2]->covg[1]);
    for (uint i=3; i!=l3.kmer_prg.nodes.size(); ++i)
    {
        EXPECT_EQ(j, l3.kmer_prg.nodes[i]->covg[0]);
	EXPECT_EQ(j, l3.kmer_prg.nodes[i]->covg[1]);
    }
     
    l3.update_covg_with_hit(mh);
    j = 2;
    EXPECT_EQ(j, l3.kmer_prg.nodes[2]->covg[0]);
    j = 0;
    EXPECT_EQ(j, l3.kmer_prg.nodes[2]->covg[1]);
    for (uint i=3; i!=l3.kmer_prg.nodes.size(); ++i)
    {
        EXPECT_EQ(j, l3.kmer_prg.nodes[i]->covg[0]);
        EXPECT_EQ(j, l3.kmer_prg.nodes[i]->covg[1]);
    }

    delete m;
    delete mr;
    delete mh;

    kh = hash.kmerhash("TAT", 3);
    m = new Minimizer(min(kh.first,kh.second), 1,4,1);
    d = {Interval(16,16), Interval(23, 26)};
    p.initialize(d);
    mr = new MiniRecord(3,p,1);
    mh = new MinimizerHit(1, m, mr);

    l3.update_covg_with_hit(mh);
    j = 1;
    EXPECT_EQ(j, l3.kmer_prg.nodes[8]->covg[1]);
    j = 2;
    EXPECT_EQ(j, l3.kmer_prg.nodes[2]->covg[0]);
    j = 0;
    EXPECT_EQ(j, l3.kmer_prg.nodes[8]->covg[0]);
    EXPECT_EQ(j, l3.kmer_prg.nodes[2]->covg[1]);
    for (uint i=3; i<8; ++i)
    {
        EXPECT_EQ(j, l3.kmer_prg.nodes[i]->covg[0]);
        EXPECT_EQ(j, l3.kmer_prg.nodes[i]->covg[1]);
    }
 
    l3.update_covg_with_hit(mh);
    j = 2;
    EXPECT_EQ(j, l3.kmer_prg.nodes[8]->covg[1]);
    EXPECT_EQ(j, l3.kmer_prg.nodes[2]->covg[0]);
    j = 0;
    EXPECT_EQ(j, l3.kmer_prg.nodes[8]->covg[0]);
    EXPECT_EQ(j, l3.kmer_prg.nodes[2]->covg[1]);
    for (uint i=3; i<8; ++i)
    {
        EXPECT_EQ(j, l3.kmer_prg.nodes[i]->covg[0]);
        EXPECT_EQ(j, l3.kmer_prg.nodes[i]->covg[1]);
    }   

    delete m;
    delete mr;
    delete mh;

    // could add futher examples for inner kmers?

    delete idx;
}

TEST_F(LocalPRGTest, write_max_path_to_fasta)
{
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");

    Index* idx;
    idx = new Index();

    l3.minimizer_sketch(idx, 1, 3);

    vector<LocalNode*> lmp3 = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};
    l3.write_max_path_to_fasta("../test/test_cases/localPRG_test.maxpath.fa", lmp3, 0.00);
    
    delete idx;
}

TEST_F(LocalPRGTest, buildVcf)
{
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");
    LocalPRG l4(4, "small real PRG", "ATGACAAAACGAAGTGGAAGTAATACGCGCAGGCGGGCTATCAGTCGCCCTGTTCGTCTGACGGCAGAAGAAGACCAGGAAATCAGAAAAAGGGCTGCTGAATGCGGCAAGACCGTTTC 5 T 6 C 5 GGTTTTTTACGGGCGGCAGCTCTCGGTAAGAAAGTTAA 7 TTCACTGACTGATGACCGAGTGCTGAAAGAAGTCATGCGACTGGGGGCGTTG 8 CTCACTGACTGATGATCGGGTACTGAAAGAAGTTATGAGACTGGGGGCGTTA 7 CAGAAAAAACTCTTTATCGACGGCAAGCGTGTCGGGGACAG 9 A 10 G 9 GAGTATGCGGAGGTGCTGAT 11 A 12 C 11 GCTATTACGGAGTATCACCG 13 G 14 T 13 GCCCTGTTATCCAGGCTTATGGCAGATTAG");

    l1.build_vcf();
    uint j = 0;
    EXPECT_EQ(j, l1.vcf.records.size());
    EXPECT_EQ(j, l1.vcf.samples.size());

    l2.build_vcf();
    j = 1;
    EXPECT_EQ(j, l2.vcf.records.size());
    EXPECT_EQ("varsite", l2.vcf.records[0].chrom);
    EXPECT_EQ((uint)1, l2.vcf.records[0].pos);
    EXPECT_EQ("GC", l2.vcf.records[0].ref);
    EXPECT_EQ("G", l2.vcf.records[0].alt);
    EXPECT_EQ("SVTYPE=INDEL;GRAPHTYPE=SIMPLE", l2.vcf.records[0].info);
    
    l3.build_vcf();
    j = 2;
    EXPECT_EQ(j, l3.vcf.records.size());
    EXPECT_EQ("nested varsite", l3.vcf.records[0].chrom);
    EXPECT_EQ((uint)1, l3.vcf.records[0].pos);
    EXPECT_EQ("GC", l3.vcf.records[0].ref);
    EXPECT_EQ("G", l3.vcf.records[0].alt);
    EXPECT_EQ("SVTYPE=INDEL;GRAPHTYPE=COMPLEX", l3.vcf.records[0].info);
    EXPECT_EQ((uint)1, l3.vcf.records[1].pos);
    EXPECT_EQ("GC", l3.vcf.records[1].ref);
    EXPECT_EQ("GT", l3.vcf.records[1].alt);
    EXPECT_EQ("SVTYPE=PH_SNPs;GRAPHTYPE=COMPLEX", l3.vcf.records[1].info);    

    l4.build_vcf();
    j = 5;
    EXPECT_EQ(j, l4.vcf.records.size());
    EXPECT_EQ("small real PRG", l4.vcf.records[0].chrom);
    EXPECT_EQ((uint)119, l4.vcf.records[0].pos);
    EXPECT_EQ("T", l4.vcf.records[0].ref);
    EXPECT_EQ("C", l4.vcf.records[0].alt);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", l4.vcf.records[0].info);

    EXPECT_EQ((uint)158, l4.vcf.records[1].pos);
    EXPECT_EQ("TTCACTGACTGATGACCGAGTGCTGAAAGAAGTCATGCGACTGGGGGCGTTG", l4.vcf.records[1].ref);
    EXPECT_EQ("CTCACTGACTGATGATCGGGTACTGAAAGAAGTTATGAGACTGGGGGCGTTA", l4.vcf.records[1].alt);
    EXPECT_EQ("SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE", l4.vcf.records[1].info);
    
    EXPECT_EQ((uint)251, l4.vcf.records[2].pos);
    EXPECT_EQ("A", l4.vcf.records[2].ref);
    EXPECT_EQ("G", l4.vcf.records[2].alt);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", l4.vcf.records[2].info);

    EXPECT_EQ((uint)272, l4.vcf.records[3].pos);
    EXPECT_EQ("A", l4.vcf.records[3].ref);
    EXPECT_EQ("C", l4.vcf.records[3].alt);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", l4.vcf.records[3].info);

    EXPECT_EQ((uint)293, l4.vcf.records[4].pos);
    EXPECT_EQ("G", l4.vcf.records[4].ref);
    EXPECT_EQ("T", l4.vcf.records[4].alt);
    EXPECT_EQ("SVTYPE=SNP;GRAPHTYPE=SIMPLE", l4.vcf.records[4].info);
}

TEST_F(LocalPRGTest, add_sample_to_vcf)
{
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");
    LocalPRG l4(4, "small real PRG", "ATGACAAAACGAAGTGGAAGTAATACGCGCAGGCGGGCTATCAGTCGCCCTGTTCGTCTGACGGCAGAAGAAGACCAGGAAATCAGAAAAAGGGCTGCTGAATGCGGCAAGACCGTTTC 5 T 6 C 5 GGTTTTTTACGGGCGGCAGCTCTCGGTAAGAAAGTTAA 7 TTCACTGACTGATGACCGAGTGCTGAAAGAAGTCATGCGACTGGGGGCGTTG 8 CTCACTGACTGATGATCGGGTACTGAAAGAAGTTATGAGACTGGGGGCGTTA 7 CAGAAAAAACTCTTTATCGACGGCAAGCGTGTCGGGGACAG 9 A 10 G 9 GAGTATGCGGAGGTGCTGAT 11 A 12 C 11 GCTATTACGGAGTATCACCG 13 G 14 T 13 GCCCTGTTATCCAGGCTTATGGCAGATTAG");
    LocalPRG l5(5, "another real PRG", " 5 ATGCTTATTGGCTATGT 7  9 ACGCGTA 10 TCGCGTA 10 ACGTGTG 9 TCAACAAATGACCAGAACAC 11 A 12 C 11  8 ACGCGTATCAACAAATGATCAGAACACA 7 GATCTACAACGTAATGCG 6 AAGT 5 ");

    vector<LocalNode*> lmp1 = {l1.prg.nodes[0]};
    l1.build_vcf();
    l1.add_sample_to_vcf(lmp1);
    uint j = 1;
    EXPECT_EQ(j, l1.vcf.samples.size());

    vector<LocalNode*> lmp2 = {l2.prg.nodes[0], l2.prg.nodes[2], l2.prg.nodes[3]};
    l2.build_vcf();
    l2.add_sample_to_vcf(lmp2);
    j = 1;
    EXPECT_EQ(j, l2.vcf.samples.size());
    EXPECT_EQ(j, l2.vcf.records[0].samples.size());
    EXPECT_EQ("1", l2.vcf.records[0].samples[0]);

    vector<LocalNode*> lmp3 = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};
    l3.build_vcf();
    l3.add_sample_to_vcf(lmp3);
    EXPECT_EQ(j, l3.vcf.samples.size());
    EXPECT_EQ(j, l3.vcf.records[0].samples.size());
    EXPECT_EQ("1", l3.vcf.records[1].samples[0]);

    vector<LocalNode*> lmp4 = {l4.prg.nodes[0], l4.prg.nodes[1], l4.prg.nodes[3], l4.prg.nodes[5], l4.prg.nodes[6], l4.prg.nodes[8], l4.prg.nodes[9], l4.prg.nodes[10], l4.prg.nodes[12], l4.prg.nodes[13], l4.prg.nodes[15]};
    l4.build_vcf();
    l4.add_sample_to_vcf(lmp4);
    EXPECT_EQ(j, l4.vcf.samples.size());
    EXPECT_EQ(j, l4.vcf.records[0].samples.size());
    EXPECT_EQ("0", l4.vcf.records[0].samples[0]);
    EXPECT_EQ(j, l4.vcf.records[1].samples.size());
    EXPECT_EQ("1", l4.vcf.records[1].samples[0]);
    EXPECT_EQ(j, l4.vcf.records[2].samples.size());
    EXPECT_EQ("1", l4.vcf.records[2].samples[0]);
    EXPECT_EQ(j, l4.vcf.records[3].samples.size());
    EXPECT_EQ("0", l4.vcf.records[3].samples[0]);
    EXPECT_EQ(j, l4.vcf.records[4].samples.size());
    EXPECT_EQ("0", l4.vcf.records[4].samples[0]);

    vector<LocalNode*> lmp5 = {l5.prg.nodes[0], l5.prg.nodes[1], l5.prg.nodes[9], l5.prg.nodes[10], l5.prg.nodes[11], l5.prg.nodes[13]};
    l5.build_vcf(); 
    l5.add_sample_to_vcf(lmp5);
    EXPECT_EQ(j, l5.vcf.samples.size());
    EXPECT_EQ((uint)7, l5.vcf.records.size());
    EXPECT_EQ(j, l5.vcf.records[0].samples.size());
    EXPECT_EQ(".", l5.vcf.records[0].samples[0]);
    EXPECT_EQ(j, l5.vcf.records[1].samples.size());
    EXPECT_EQ(".", l5.vcf.records[1].samples[0]);
    EXPECT_EQ(j, l5.vcf.records[2].samples.size());
    EXPECT_EQ("1", l5.vcf.records[2].samples[0]);
    EXPECT_EQ(j, l5.vcf.records[3].samples.size());
    EXPECT_EQ(".", l5.vcf.records[3].samples[0]);
    EXPECT_EQ(j, l5.vcf.records[4].samples.size());
    EXPECT_EQ(".", l5.vcf.records[4].samples[0]);
    EXPECT_EQ(j, l5.vcf.records[5].samples.size());
    EXPECT_EQ(".", l5.vcf.records[5].samples[0]);
    EXPECT_EQ(j, l5.vcf.records[6].samples.size());
    EXPECT_EQ(".", l5.vcf.records[6].samples[0]);
}

TEST_F(LocalPRGTest, moreupdateVCF)
{
    // load PRGs from file
    vector<LocalPRG*> prgs;
    read_prg_file(prgs, "../test/test_cases/updatevcf_test.fa");

    EXPECT_EQ((uint)3, prgs.size());
    prgs[0]->build_vcf();
    prgs[1]->build_vcf();
    prgs[2]->build_vcf();

    for (uint i=0; i!=prgs[2]->vcf.records.size(); ++i)
    {
	cout << prgs[2]->vcf.records[i];
    }

    vector<LocalNode*> lmp1 = {prgs[1]->prg.nodes[0], prgs[1]->prg.nodes[11], prgs[1]->prg.nodes[12], prgs[1]->prg.nodes[17], prgs[1]->prg.nodes[65], prgs[1]->prg.nodes[67]};
    cout << "PRG 1 has " << prgs[1]->prg.nodes.size() << " nodes" << endl;
    prgs[1]->add_sample_to_vcf(lmp1);

    vector<LocalNode*> lmp2 = {prgs[2]->prg.nodes[0], prgs[2]->prg.nodes[1], prgs[2]->prg.nodes[3], prgs[2]->prg.nodes[4], prgs[2]->prg.nodes[6], prgs[2]->prg.nodes[7],
			       prgs[2]->prg.nodes[9], prgs[2]->prg.nodes[10], prgs[2]->prg.nodes[11], prgs[2]->prg.nodes[13], prgs[2]->prg.nodes[14], prgs[2]->prg.nodes[16],
			       prgs[2]->prg.nodes[17], prgs[2]->prg.nodes[19], prgs[2]->prg.nodes[44], prgs[2]->prg.nodes[45], prgs[2]->prg.nodes[47], prgs[2]->prg.nodes[118],
			       prgs[2]->prg.nodes[119], prgs[2]->prg.nodes[121], prgs[2]->prg.nodes[123], prgs[2]->prg.nodes[125], prgs[2]->prg.nodes[126], 
			       prgs[2]->prg.nodes[130], prgs[2]->prg.nodes[131], prgs[2]->prg.nodes[133], prgs[2]->prg.nodes[135], prgs[2]->prg.nodes[141], 
			       prgs[2]->prg.nodes[142], prgs[2]->prg.nodes[144], prgs[2]->prg.nodes[145], prgs[2]->prg.nodes[160]};
    cout << "PRG 2 has " << prgs[2]->prg.nodes.size() << " nodes" << endl;
    prgs[2]->add_sample_to_vcf(lmp2);
}

TEST_F(LocalPRGTest, find_path_and_variants)
{

    Index* idx;
    idx = new Index();

    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");
    l3.minimizer_sketch(idx, 1, 3);

    cout << l3.kmer_prg;  
    delete idx;
}
