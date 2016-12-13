#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "localPRG.h"
#include "minimizer.h"
#include "minirecord.h"
#include "interval.h"
#include "path.h"
#include "localgraph.h"
#include "localnode.h"
#include "index.h"
#include "inthash.h"
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
    cout << "define" << endl;
    //LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");
    
    // empty interval expects no nodes along
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    vector<LocalNode*> v;
    cout << "test" << endl;
    //EXPECT_EQ(v, l0.nodes_along_path(p));
    EXPECT_EQ(v, l1.nodes_along_path(p));
    /*EXPECT_EQ(v, l2.nodes_along_path(p));
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
    EXPECT_EQ(j, l3.nodes_along_path(p)[2]->id);

    // including empty interval
    d = {Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    j = 3;
    EXPECT_EQ(j, l3.string_along_path(p).size());
    j = 3;
    EXPECT_EQ(j, l3.nodes_along_path(p)[0]->id);
    j = 4;
    EXPECT_EQ(j, l3.nodes_along_path(p)[1]->id);
    j = 6;
    EXPECT_EQ(j, l3.nodes_along_path(p)[2]->id);

    // and a path that can't really exist still works
    d = {Interval(12,13),Interval(19,20)};
    p.initialize(d);
    j = 2;
    EXPECT_EQ(j, l3.nodes_along_path(p).size());
    j = 3;
    EXPECT_EQ(j, l3.nodes_along_path(p)[0]->id);
    j = 5;
    EXPECT_EQ(j, l3.nodes_along_path(p)[2]->id);*/
}

TEST_F(LocalPRGTest, splitBySiteNoSites){
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    vector<Interval> v0, v1;
    v0.push_back(Interval(0,0));
    EXPECT_ITERABLE_EQ( vector< Interval >, v0, l0.splitBySite(Interval(0,0)));// << "Failed to split empty string with input Interval";
    v1.push_back(Interval(0,4));
    EXPECT_ITERABLE_EQ( vector< Interval >, v1, l1.splitBySite(Interval(0,4)));// << "Failed to split string with input Interval";
    v1.clear();
    v1.push_back(Interval(0,2));
    EXPECT_ITERABLE_EQ( vector< Interval >, v1, l1.splitBySite(Interval(0,2)));// << "Failed to split string with short input Interval";
    v1.clear();
    v1.push_back(Interval(1,3));
    EXPECT_ITERABLE_EQ( vector< Interval >, v1, l1.splitBySite(Interval(1,3)));// << "Failed to split string with middle input Interval";
}

TEST_F(LocalPRGTest, splitBySiteSite){
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    vector<Interval> v2;
    v2.push_back(Interval(0,1));
    l2.next_site = 5;
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,1)));// << "Failed to split string in Interval" << Interval(0,1);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,2)));// << "Failed to split string in Interval" << Interval(0,2);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,3)));// << "Failed to split string in Interval" << Interval(0,3);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,4)));// << "Failed to split string in Interval" << Interval(0,4);
    v2.push_back(Interval(4,6));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,6)));// << "Failed to split string in Interval" << Interval(0,6);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,7)));// << "Failed to split string in Interval" << Interval(0,7);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,8)));// << "Failed to split string in Interval" << Interval(0,8);
    v2.push_back(Interval(9,10));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,10)));// << "Failed to split string in Interval" << Interval(0,10);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,11)));// << "Failed to split string in Interval" << Interval(0,11);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,12)));// << "Failed to split string in Interval" << Interval(0,12);
    v2.push_back(Interval(13,14));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,14)));// << "Failed to split string in Interval" << Interval(0,14);
    v2.clear();
    v2.push_back(Interval(5,6));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(5,8)));// << "Failed to split string in mid Interval" << Interval(5,8);
}

TEST_F(LocalPRGTest, splitBySiteNestedSite){
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    vector<Interval> v3;
    v3.push_back(Interval(0,1));
    l3.next_site = 5;
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,1)));// << "Failed to split string in Interval" << Interval(0,1);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,2)));// << "Failed to split string in Interval" << Interval(0,2);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,3)));// << "Failed to split string in Interval" << Interval(0,3);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,4)));// << "Failed to split string in Interval" << Interval(0,4);
    v3.push_back(Interval(4,16));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,16)));// << "Failed to split string in Interval" << Interval(0,6);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,17)));// << "Failed to split string in Interval" << Interval(0,7);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,18)));// << "Failed to split string in Interval" << Interval(0,8);
    v3.push_back(Interval(19,20));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,20)));// << "Failed to split string in Interval" << Interval(0,10);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,21)));// << "Failed to split string in Interval" << Interval(0,11);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,22)));// << "Failed to split string in Interval" << Interval(0,12);
    v3.push_back(Interval(23,24));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,24)));// << "Failed to split string in Interval" << Interval(0,14);
    l3.next_site = 7;    
    v3.clear();
    v3.push_back(Interval(4,5));
    v3.push_back(Interval(8,9));
    v3.push_back(Interval(12,13));
    v3.push_back(Interval(16,16));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(4,16)));// << "Failed to split string in mid Interval" << Interval(5,8);
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
    uint32_t j = 0;
    EXPECT_EQ(j, l0.max_path_index.size());

    LocalGraph lg1;
    lg1.add_node(0,"AGCT", Interval(0,4));
    EXPECT_EQ(lg1, l1.prg);
    j = 1;
    EXPECT_EQ(j, l1.prg.index.size());
    EXPECT_EQ(j, l1.prg.index[0].size());

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
    j = 2;
    EXPECT_EQ(j, l2.prg.index.size());
    j = 1;
    EXPECT_EQ(j, l2.prg.index[0].size());
    EXPECT_EQ(j, l2.prg.index[1].size());

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
    j = 3;
    EXPECT_EQ(j, l3.prg.index.size());
    j = 1;
    EXPECT_EQ(j, l3.prg.index[0].size());
    EXPECT_EQ(j, l3.prg.index[1].size());
    EXPECT_EQ(j, l3.prg.index[2].size());
    j = 0;
    EXPECT_EQ(j, l3.prg.index[0][0].first);
    EXPECT_EQ(j, l3.prg.index[1][0].first);
    j = 6;
    EXPECT_EQ(j, l3.prg.index[0][0].second);
    EXPECT_EQ(j, l3.prg.index[1][0].second);
    j = 1;
    EXPECT_EQ(j, l3.prg.index[2][0].first);
    j = 4;
    EXPECT_EQ(j, l3.prg.index[2][0].second);
}

//could test that no int is included in more than one node
TEST_F(LocalPRGTest, minimizerSketch){
    // note this is a bad test
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    Index* idx;
    idx = new Index();

    l0.minimizer_sketch(idx, 1, 3);
    uint32_t j = 0;
    EXPECT_EQ(j, idx->minhash.size());

    l1.minimizer_sketch(idx, 2, 3);
    j = 1;
    EXPECT_EQ(j, idx->minhash.size());
    l1.minimizer_sketch(idx, 1, 3);
    EXPECT_EQ(j, idx->minhash.size());
    j = 2;
    pair<uint64_t,uint64_t> kh = kmerhash("AGC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)].size());
    
    idx->clear();
    l2.minimizer_sketch(idx, 2, 3);
    j = 1;
    EXPECT_EQ(j, idx->minhash.size());
    l2.minimizer_sketch(idx, 1, 3);
    //j = 2;
    EXPECT_EQ(j, idx->minhash.size());
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)].size());

    idx->clear();
    l3.minimizer_sketch(idx, 2, 3);
    j = 2;
    EXPECT_EQ(j, idx->minhash.size());
    l3.minimizer_sketch(idx, 1, 3);
    //j = 3;
    EXPECT_EQ(j, idx->minhash.size());
    j = 1;
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)].size()); //AGC, GCT
    kh = kmerhash("AGT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)].size()); //AGTx2
    j = 0;
    kh = kmerhash("GTT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)].size());
    
    delete idx;
}

TEST_F(LocalPRGTest, updateCovgWithHit)
{
}

TEST_F(LocalPRGTest, inferMostLikelyPrgPathsForCorrespondingPannode)
{
}

TEST_F(LocalPRGTest, writeFasta)
{
}

// this version does not apply to forcing minis to be disjoint (on a path)
/*TEST_F(LocalPRGTest, minimizerSketch){
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    Index* idx;
    idx = new Index();

    l0.minimizer_sketch(idx, 1, 3);
    uint32_t j = 0;
    EXPECT_EQ(j, idx->minhash.size());

    l1.minimizer_sketch(idx, 2, 3);
    j = 1;
    EXPECT_EQ(j, idx->minhash.size());
    l1.minimizer_sketch(idx, 1, 3);
    EXPECT_EQ(j, idx->minhash.size());
    j = 2;
    pair<uint64_t,uint64_t> kh = kmerhash("AGC",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)].size());
    
    idx->clear();
    l2.minimizer_sketch(idx, 2, 3);
    j = 1;
    EXPECT_EQ(j, idx->minhash.size());
    l2.minimizer_sketch(idx, 1, 3);
    j = 2;
    EXPECT_EQ(j, idx->minhash.size());
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)].size());

    idx->clear();
    l3.minimizer_sketch(idx, 2, 3);
    j = 2;
    EXPECT_EQ(j, idx->minhash.size());
    l3.minimizer_sketch(idx, 1, 3);
    j = 3;
    EXPECT_EQ(j, idx->minhash.size());
    j = 2;
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)].size()); //AGC, GCT
    kh = kmerhash("AGT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)].size()); //AGTx2
    j = 1;
    kh = kmerhash("GTT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)].size());
    
    delete idx;
}*/

/*TEST_F(LocalPRGTest, getCovgs)
{
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T"); 

    MinimizerHits* mh;
    mh = new MinimizerHits();
    Minimizer* m;
    pair<uint64_t,uint64_t> kh = kmerhash("AGC",3);
    m = new Minimizer(min(kh.first,kh.second),0,3,(min(kh.first,kh.second)==kh.first));
    MiniRecord* r;
    Path p;
    deque<Interval> d = {Interval(0,3)};
    p.initialize(d);
    r = new MiniRecord(1, p, 0);
    mh->add_hit(0, m, r);
    l1.get_covgs(mh);
    uint32_t j = 3;
    EXPECT_EQ(j, l1.prg.nodes[0]->covg);
    // add a second identical read
    mh->add_hit(1, m, r);
    l1.get_covgs(mh);
    j = 6;
    EXPECT_EQ(j, l1.prg.nodes[0]->covg);

    delete r;
    d = {Interval(0,1), Interval(4,6)};
    p.initialize(d);
    r = new MiniRecord(2, p, 0);
    mh->add_hit(0, m, r);
    l2.get_covgs(mh);
    j = 1;
    EXPECT_EQ(j, l2.prg.nodes[0]->covg);
    j = 2;
    EXPECT_EQ(j, l2.prg.nodes[1]->covg);
    j = 0;
    EXPECT_EQ(j, l2.prg.nodes[2]->covg);
    EXPECT_EQ(j, l2.prg.nodes[3]->covg);


    delete r;
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    r = new MiniRecord(3, p, 0);
    mh->add_hit(0, m, r);
    l3.get_covgs(mh);
    j = 1;
    EXPECT_EQ(j, l3.prg.nodes[0]->covg);
    EXPECT_EQ(j, l3.prg.nodes[1]->covg);
    EXPECT_EQ(j, l3.prg.nodes[2]->covg);
    j = 0;
    EXPECT_EQ(j, l3.prg.nodes[3]->covg);
    EXPECT_EQ(j, l3.prg.nodes[4]->covg);
    EXPECT_EQ(j, l3.prg.nodes[5]->covg);
    EXPECT_EQ(j, l3.prg.nodes[6]->covg);
    // add same hit again
    mh->add_hit(1, m, r);
    l3.get_covgs(mh);
    j = 2;
    EXPECT_EQ(j, l3.prg.nodes[0]->covg);
    EXPECT_EQ(j, l3.prg.nodes[1]->covg);
    EXPECT_EQ(j, l3.prg.nodes[2]->covg);
    j = 0;
    EXPECT_EQ(j, l3.prg.nodes[3]->covg);
    EXPECT_EQ(j, l3.prg.nodes[4]->covg);
    EXPECT_EQ(j, l3.prg.nodes[5]->covg);
    EXPECT_EQ(j, l3.prg.nodes[6]->covg);
    // add a different hit
    delete r;
    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    r = new MiniRecord(3, p, 0);
    mh->add_hit(1, m, r);
    l3.get_covgs(mh);
    j = 3;
    EXPECT_EQ(j, l3.prg.nodes[0]->covg);
    EXPECT_EQ(j, l3.prg.nodes[1]->covg);
    j = 2;
    EXPECT_EQ(j, l3.prg.nodes[2]->covg);
    j = 1;
    EXPECT_EQ(j, l3.prg.nodes[3]->covg);
    j = 0;
    EXPECT_EQ(j, l3.prg.nodes[4]->covg);
    EXPECT_EQ(j, l3.prg.nodes[5]->covg);
    EXPECT_EQ(j, l3.prg.nodes[6]->covg);
    // and a hit crossing the empty string node
    delete r;
    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    r = new MiniRecord(3, p, 0);
    mh->add_hit(2, m, r);
    l3.get_covgs(mh);
    j = 4;
    EXPECT_EQ(j, l3.prg.nodes[1]->covg);
    j = 3;
    EXPECT_EQ(j, l3.prg.nodes[0]->covg);
    j = 2;
    EXPECT_EQ(j, l3.prg.nodes[2]->covg);
    EXPECT_EQ(j, l3.prg.nodes[3]->covg);
    j = 1;
    EXPECT_EQ(j, l3.prg.nodes[6]->covg);
    j = 0;
    EXPECT_EQ(j, l3.prg.nodes[4]->covg);
    EXPECT_EQ(j, l3.prg.nodes[5]->covg);
    delete r;
    delete m;
    delete mh;
}*/

/*TEST_F(LocalPRGTest,unpackLinearString){
    Seq s1 = Seq(0,"0", "AGCTAATGCGTT", 11, 3);
    Seq s2 = Seq(0,"0", "AGCTAATGCGTT", 10, 3);
    Seq s3 = Seq(0,"0", "AGCTAATGCGTT", 9, 3);
    Seq s4 = Seq(0,"0", "AGCTAGTGCGTT", 9, 3);
    uint32_t j = 0;
    EXPECT_EQ(s1.sketch.size(),j) << "Have " << s1.sketch.size() << " minimizer when string is too short";
    ++j;
    EXPECT_EQ(s2.sketch.size(),j) << "Have " << s2.sketch.size() << " minimizers when should have 1";
    EXPECT_EQ(s3.sketch.size(),j) << "Have " << s3.sketch.size() << " minimizers when should have 1";
    ++j;
    EXPECT_EQ(s4.sketch.size(),j) << "Have " << s4.sketch.size() << " minimizers when should have 2";
}

TEST_F(SeqTest,sketchIncludesEveryLetter){
    Seq s1 = Seq(0,"0", "AGCTAATGTGTT", 3, 3);
    Seq s2 = Seq(0,"0", "AGCTAATGTGTT", 2, 3);
    Seq s3 = Seq(0,"0", "AGCTAATGTGTT", 1, 3);
    Seq s4 = Seq(0,"0", "AGCTAATGTGAT", 3, 3);

    set<int> pos_inc;
    for(set<Minimizer*>::iterator it=s4.sketch.begin(); it != s4.sketch.end(); ++it)
    {
        for (std::deque<Interval>::iterator it2=((*it)->path).path.begin(); it2!=((*it)->path).path.end(); ++it2)
        {
            for (uint32_t j = it2->start; j<it2->end; ++j)
            {
                pos_inc.insert(j);
            }
        }
    } 
    set<int> expected = {0,1,2,3,4,5,6,7,8,9,10,11};
    EXPECT_EQ(pos_inc, expected) << "sketch misses a letter";
 
    uint32_t j = 10;
    EXPECT_EQ(s3.sketch.size(), j) << "sketch with w=1 has incorrect size " << s3.sketch.size();
    
    pos_inc.clear();
    for(set<Minimizer*>::iterator it=s2.sketch.begin(); it != s2.sketch.end(); ++it)
    {
        for (std::deque<Interval>::iterator it2=((*it)->path).path.begin(); it2!=((*it)->path).path.end(); ++it2)
        {
            for (uint32_t j = it2->start; j<it2->end; ++j)
            {
                pos_inc.insert(j);
            }
        }
    }
    EXPECT_EQ(pos_inc, expected) << "sketch for s2 includes/misses wrong letter";

    pos_inc.clear();
    for(set<Minimizer*>::iterator it=s1.sketch.begin(); it != s1.sketch.end(); ++it)
    {
        for (std::deque<Interval>::iterator it2=((*it)->path).path.begin(); it2!=((*it)->path).path.end(); ++it2)
        {
            for (uint32_t j = it2->start; j<it2->end; ++j)
            {
                pos_inc.insert(j);
            }
        }
    }
    expected = {0,1,2,3,4,5,6,7,8,9};
    EXPECT_EQ(pos_inc, expected) << "sketch for s1 includes/misses wrong letter";

}*/

//TEST_F(SeqTest,sketchCorrect){
//}
