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
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    deque<Interval> d = {Interval(0,3)};
    Path p, q;
    p.initialize(d);
    d = {Interval(1,4)};
    q.initialize(d);
    vector<Path> v_exp = {q};
    EXPECT_ITERABLE_EQ(vector<Path>, v_exp, l1.shift(p));

}

TEST_F(LocalPRGTest, minimizerSketch){
    // note this is a bad test
    LocalPRG l0(0,"empty", "");
    LocalPRG l1(1,"simple", "AGCT");
    LocalPRG l2(2,"varsite", "A 5 GC 6 G 5 T");
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

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
    j = 2;
    EXPECT_EQ(j, idx->minhash.size());
    l2.minimizer_sketch(idx, 1, 3);
    EXPECT_EQ(j, idx->minhash.size());
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());
    j = 1;
    kh = hash.kmerhash("AGT",3);
    EXPECT_EQ(j, idx->minhash[min(kh.first,kh.second)]->size());

    idx->clear();
    l3.minimizer_sketch(idx, 2, 3);
    j = 3;
    EXPECT_EQ(j, idx->minhash.size());
    l3.minimizer_sketch(idx, 1, 3);
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
    delete idx;
}

TEST_F(LocalPRGTest, updateCovgWithHit)
{
// do need a test for this, but function currently altered to only update on kmers, not for localnodes
/*    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 TAT");

    KmerHash hash;

    Minimizer* m;
    pair<uint64_t,uint64_t> kh = hash.kmerhash("AGC", 3);
    m = new Minimizer(min(kh.first,kh.second), 1,4,0);
    deque<Interval> d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    Path p;
    p.initialize(d);
    MiniRecord* mr;
    mr = new MiniRecord(3,p,0);
    MinimizerHit* mh;
    mh = new MinimizerHit(1, m, mr);

    l3.update_covg_with_hit(mh);
    uint j = 2;
    EXPECT_EQ(j, l3.prg.nodes[0]->covg);
    EXPECT_EQ(j, l3.prg.nodes[1]->covg);
    EXPECT_EQ(j, l3.prg.nodes[2]->covg);
    j = 0;
    EXPECT_EQ(j, l3.prg.nodes[4]->covg);
    j = 1;
    EXPECT_EQ(j, l3.prg.nodes[3]->covg);
    EXPECT_EQ(j, l3.prg.nodes[5]->covg);
    j = 3;
    EXPECT_EQ(j, l3.prg.nodes[6]->covg);
    
     
    delete m;
    delete mr;
    delete mh;
    kh = hash.kmerhash("CTA", 3); 
    m = new Minimizer(min(kh.first,kh.second), 3,6,0);
    d = {Interval(8, 9), Interval(16, 16), Interval(23, 25)};
    p.initialize(d);
    mr = new MiniRecord(3,p,0);
    mh = new MinimizerHit(1, m, mr);

    l3.update_covg_with_hit(mh);
    j = 2;
    EXPECT_EQ(j, l3.prg.nodes[0]->covg);
    EXPECT_EQ(j, l3.prg.nodes[1]->covg);
    j = 3;
    EXPECT_EQ(j, l3.prg.nodes[2]->covg);
    j = 1;
    EXPECT_EQ(j, l3.prg.nodes[4]->covg);
    EXPECT_EQ(j, l3.prg.nodes[3]->covg);
    EXPECT_EQ(j, l3.prg.nodes[5]->covg);
    j = 5;
    EXPECT_EQ(j, l3.prg.nodes[6]->covg);

    delete m;
    delete mr;
    delete mh;*/
}

/*TEST_F(LocalPRGTest, inferMostLikelyPrgPathsForCorrespondingPannode)
{
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();

    KmerHash hash;

    // initialize a prgs object
    vector<LocalPRG*> prgs;
    LocalPRG* lp3;
    lp3 = new LocalPRG(3, "3", "T 5 G 7 C 8 T 7  6 G 5 TATG");
    prgs.push_back(lp3);

    // initialize index as we would expect with example prgs
    Index *idx;
    idx = new Index();

    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    pair<uint64_t,uint64_t> kh = hash.kmerhash("TGC",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    kh = hash.kmerhash("TGT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kh = hash.kmerhash("TTT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = hash.kmerhash("GCT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = hash.kmerhash("GTT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(8,9), Interval(16,16), Interval(23,25)};
    p.initialize(d);
    kh = hash.kmerhash("CTA",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(12,13), Interval(16,16), Interval(23,25)};
    p.initialize(d);
    kh = hash.kmerhash("TTA",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(19,20), Interval(23,25)};
    p.initialize(d);
    kh = hash.kmerhash("TTA",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(23,26)};
    p.initialize(d);
    kh = hash.kmerhash("TAT",3); //inconsistent, i don't care
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(24,27)};
    p.initialize(d);
    kh = hash.kmerhash("ATG",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(27,27)};
    p.initialize(d);
    lp3->kmer_prg.add_node(p);
    
    PanNode* pn;
    pn = new PanNode(3);
    pn->add_read(0);
    add_read_hits(0, "read0", "TGTTATG", mhs, idx, 1, 3); //AGTTAAGCTAGCTACTTACGGTA
    pn->add_hits(mhs->hits);
    
    lp3->infer_most_likely_prg_paths_for_corresponding_pannode(pn, 3, 0.0015);
    vector<LocalNode*> v_exp = {lp3->prg.nodes[0], lp3->prg.nodes[1], lp3->prg.nodes[3], lp3->prg.nodes[4], lp3->prg.nodes[6]};
    EXPECT_EQ(v_exp.size(), lp3->max_path_index[0][0][0].npath.size());
    for (uint j = 0; j!= v_exp.size(); ++j)
    {
	EXPECT_EQ(*(v_exp[j]), *(lp3->max_path_index[0][0][0].npath[j]));
    }

    delete mhs;
    idx->clear();
    delete idx;
    delete pn;
}*/

TEST_F(LocalPRGTest, writeFasta)
{
}
