#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "utils.h"
#include "localPRG.h"
#include "pangraph.h"
#include "interval.h"
#include "path.h"
#include "minihit.h"
#include "minihits.h"
#include "index.h"
#include "inthash.h"
#include <stdint.h>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

class UtilsTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }

};

TEST_F(UtilsTest, split)
{
    vector<string> v = {"abc", "def", "ghi"};
    EXPECT_EQ(v, split("abc, def, ghi", ", "));
    EXPECT_EQ(v, split("abc, def, ghi, ", ", "));
    EXPECT_EQ(v, split(", abc, def, ghi", ", "));    
}

TEST_F(UtilsTest, revComplement)
{
    string s = "ACCTGATTGCGTA";
    EXPECT_EQ(s, rev_complement(rev_complement(s)));
    
    string t = "TACGCAATCAGGT";
    EXPECT_EQ(t, rev_complement(s));

    s = "ACCTGATTgCGTA";
    EXPECT_EQ(t, rev_complement(s));
    
    s = "ACCTGATTYCGTA";
    t = "TACGNAATCAGGT";
    EXPECT_EQ(t, rev_complement(s));
}

// don't bother with nchoosek test as will remove function

TEST_F(UtilsTest, readPrgFile){
    vector<LocalPRG*> prgs;

    // simple case first, single prg with empty string sequence
    // doesn't get added to prgs 
    read_prg_file(prgs, "../test/test_cases/prg0.fa");
    uint32_t j = 0;
    EXPECT_EQ(prgs.size(), j);
 
    // single prg with simple sequence
    read_prg_file(prgs, "../test/test_cases/prg1.fa");
    LocalPRG l1(1,"prg1", "AGCT");
    j = 1;
    EXPECT_EQ(prgs.size(), j);
    j = 0;
    EXPECT_EQ(prgs[0]->id, j);
    EXPECT_EQ(prgs[0]->name, "prg1");
    EXPECT_EQ(prgs[0]->seq, "AGCT");
    EXPECT_EQ(prgs[0]->prg, l1.prg);

    // single prg with a variant site
    read_prg_file(prgs, "../test/test_cases/prg2.fa");
    LocalPRG l2(2,"prg2", "A 5 GC 6 G 5 T");
    j = 2;
    EXPECT_EQ(prgs.size(), j);
    j = 0;
    EXPECT_EQ(prgs[1]->id, j);
    EXPECT_EQ(prgs[1]->name, "prg2");
    EXPECT_EQ(prgs[1]->seq, "A 5 GC 6 G 5 T");
    EXPECT_EQ(prgs[1]->prg, l2.prg);

    // single prg with a nested variant site
    read_prg_file(prgs, "../test/test_cases/prg3.fa");
    LocalPRG l3 = LocalPRG(3,"prg3", "A 5 G 7 C 8 T 7  6 G 5 T");
    j = 3;
    EXPECT_EQ(prgs.size(), j);
    j = 0;
    EXPECT_EQ(prgs[2]->id, j);
    EXPECT_EQ(prgs[2]->name, "prg3");
    EXPECT_EQ(prgs[2]->seq, "A 5 G 7 C 8 T 7  6 G 5 T");
    EXPECT_EQ(prgs[2]->prg, l3.prg);    

    // now a prg input file with all 4 in
    prgs.clear();
    EXPECT_EQ(prgs.size(), j);
    read_prg_file(prgs, "../test/test_cases/prg0123.fa");
    j = 3;
    EXPECT_EQ(prgs.size(), j);
}

TEST_F(UtilsTest, addReadHits){
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    MinimizerHits expected1;
    MinimizerHits expected2;
    MinimizerHits expected3;
    MinimizerHits expected4;
    MinimizerHit *m1, *m2, *m3, *m4, *m5, *m6, *m7, *m8, *m9, *m10, *m11;

    // initialize index as we would expect with example prgs 1 and 3 from above
    KmerHash hash;
    Index *idx;
    idx = new Index();
    deque<Interval> d = {Interval(0,3)};
    Path p;
    p.initialize(d);
    pair<uint64_t,uint64_t> kh = hash.kmerhash("AGC",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first<kh.second));
    m1 = new MinimizerHit(0, Interval(0,3), 1, p, 1);
    m2 = new MinimizerHit(0, Interval(1,4), 1, p, 0);
    expected1.hits.insert(m1);
    expected2.hits.insert(m2);
    d = {Interval(1,4)};
    p.initialize(d);
    kh = hash.kmerhash("GCT",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first<kh.second));
    m3 = new MinimizerHit(0, Interval(1,4), 1, p, 1);
    m4 = new MinimizerHit(0, Interval(0,3), 1, p, 0);
    expected2.hits.insert(m3);
    expected1.hits.insert(m4);
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    kh = hash.kmerhash("AGC",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first<kh.second));
    m5 = new MinimizerHit(0, Interval(0,3), 3, p, 1);
    m6 = new MinimizerHit(0, Interval(1,4), 3, p, 0);
    expected1.hits.insert(m5);
    expected2.hits.insert(m6);
    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    kh = hash.kmerhash("AGT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first<kh.second));
    m9 = new MinimizerHit(0, Interval(0,3), 3, p, 1);
    expected3.hits.insert(m9);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first<kh.second));
    m10 = new MinimizerHit(0, Interval(0,3), 3, p, 1);
    expected3.hits.insert(m10);
    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = hash.kmerhash("GCT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first<kh.second));
    m7 = new MinimizerHit(0, Interval(1,4), 3, p, 1);
    m8 = new MinimizerHit(0, Interval(0,3), 3, p, 0);
    expected2.hits.insert(m7);
    expected1.hits.insert(m8);
    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = hash.kmerhash("GTT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first<kh.second));
    m11 = new MinimizerHit(0, Interval(1,4), 3, p, 1);
    expected4.hits.insert(m11);

    Seq *s;
    s = new Seq(0, "read1", "AGC", 1, 3);
    add_read_hits(s, mhs, idx);
    EXPECT_EQ(expected1.hits.size(), mhs->hits.size());
    set<MinimizerHit*, pComp>::const_iterator it2 = expected1.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }
    
    // if take w=2 as sketch of read AGTT should miss AGT, which occurs twice in PRG and contain GTT
    delete mhs;
    delete s;
    mhs = new MinimizerHits();
    uint32_t j = 0;
    EXPECT_EQ(j, mhs->hits.size());
    s = new Seq(0, "read2", "AGTT", 2, 3);
    add_read_hits(s, mhs, idx);
    EXPECT_EQ(expected4.hits.size(), mhs->hits.size());
    it2 = expected4.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    // but for w=1, only add one more hit, for GTT
    expected3.hits.insert(m11);
    delete mhs;
    delete s;
    mhs = new MinimizerHits();
    EXPECT_EQ(j, mhs->hits.size());
    s = new Seq(0, "read2", "AGTT", 1, 3);
    add_read_hits(s, mhs, idx);
    EXPECT_EQ(expected3.hits.size(), mhs->hits.size());
    it2 = expected3.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    // now back to w = 1, add expected2 to expected1 as will get hits against both AGC and GCT
    delete mhs;
    delete s;
    mhs = new MinimizerHits();
    j = 0;
    EXPECT_EQ(j, mhs->hits.size());
    s = new Seq(0, "read3", "AGCT", 1, 3);
    add_read_hits(s, mhs, idx);
    expected1.hits.insert(expected2.hits.begin(), expected2.hits.end());
    EXPECT_EQ(expected1.hits.size(), mhs->hits.size());
    it2 = expected1.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    // same for w = 2, add expected2 to expected1 as will get hits against both because AGC and GCT are joint minimums
    delete mhs;
    delete s;
    mhs = new MinimizerHits();
    j = 0;
    EXPECT_EQ(j, mhs->hits.size());
    s = new Seq(0, "read3", "AGCT", 2, 3);
    add_read_hits(s, mhs, idx);
    EXPECT_EQ(expected1.hits.size(), mhs->hits.size());
    it2 = expected1.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }    

    expected1.hits.clear();
    expected2.hits.clear();
    expected3.hits.clear();
    expected4.hits.clear();
    idx->clear();
    delete idx;
    delete m1;
    delete m2;
    delete m3;
    delete m4;
    delete m5;
    delete m6;
    delete m7;
    delete m8;
    delete m9;
    delete m10;
    delete m11;
    delete mhs;
    delete s;
}

TEST_F(UtilsTest, simpleInferLocalPRGOrderForRead){    
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    KmerHash hash;

    // initialize a prgs object
    vector<LocalPRG*> prgs;
    LocalPRG* lp1;
    LocalPRG* lp3;
    lp1 = new LocalPRG(1, "1", "");
    lp3 = new LocalPRG(0, "0", "");
    prgs.push_back(lp3);
    prgs.push_back(lp1);

    // initialize index as we would expect with example prgs (variant of) 1 and 3 from above
    Index *idx;
    idx = new Index();

    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    lp1->kmer_prg.add_node(p);

    d = {Interval(0,3)};
    p.initialize(d);
    pair<uint64_t,uint64_t> kh = hash.kmerhash("TAC",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(1,4)};
    p.initialize(d);
    kh = hash.kmerhash("ACG",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(4,4)};
    p.initialize(d);
    lp1->kmer_prg.add_node(p);

    d = {Interval(0,0)};
    p.initialize(d);
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    kh = hash.kmerhash("AGC",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    kh = hash.kmerhash("AGT",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kh = hash.kmerhash("ATT",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = hash.kmerhash("GCT",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = hash.kmerhash("GTT",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(12,13), Interval(16,16), Interval(23,25)};
    p.initialize(d);
    kh = hash.kmerhash("TTA",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(23,26)};
    p.initialize(d);
    kh = hash.kmerhash("TAA",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(24,27)};
    p.initialize(d);
    kh = hash.kmerhash("AAG",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(27,27)};
    p.initialize(d);
    lp3->kmer_prg.add_node(p);

    // add read hits to mhs
    Seq *s;
    s = new Seq(0, "read1", "AGTTAAGTACG", 1, 3);
    add_read_hits(s, mhs, idx);
    delete s;
    //add_read_hits(0, "read1", "AGTTAAGTACG", mhs, idx, 1, 3);

    // initialize pangraph;
    PanGraph *pg;
    pg = new PanGraph();
    infer_localPRG_order_for_reads(prgs, mhs, pg, 1, 3);

    // create a pangraph object representing the truth we expect (prg 3 then 1)
    PanGraph pg_exp;
    MinimizerHits mhs_dummy;
    pg_exp.add_node(1,0, mhs_dummy.hits);
    pg_exp.add_node(0,0, mhs_dummy.hits);
    pg_exp.add_edge(0,1);

    EXPECT_EQ(pg_exp, *pg);
    idx->clear();
    delete idx;
    delete pg;
    delete lp1;
    delete lp3;
    delete mhs;
}

TEST_F(UtilsTest, biggerInferLocalPRGOrderForRead){
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    KmerHash hash;

    // initialize a prgs object
    vector<LocalPRG*> prgs;
    LocalPRG* lp1;
    LocalPRG* lp2;
    LocalPRG* lp3;
    LocalPRG* lp4;
    lp1 = new LocalPRG(1, "1", "");
    lp3 = new LocalPRG(3, "3", "");
    lp4 = new LocalPRG(0, "", "");
    lp2 = new LocalPRG(2, "2", "");
    prgs.push_back(lp4);
    prgs.push_back(lp1);
    prgs.push_back(lp2);
    prgs.push_back(lp3);

    // initialize index as we would expect with example prgs
    Index *idx;
    idx = new Index();

    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    lp1->kmer_prg.add_node(p);

    d = {Interval(0,3)};
    p.initialize(d);
    pair<uint64_t,uint64_t> kh = hash.kmerhash("TAC",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(1,4)};
    p.initialize(d);
    kh = hash.kmerhash("ACG",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(2,5)};
    p.initialize(d);
    kh = hash.kmerhash("CGG",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(3,6)};
    p.initialize(d);
    kh = hash.kmerhash("GGT",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(4,7)};
    p.initialize(d);
    kh = hash.kmerhash("GTA",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(7,7)};
    p.initialize(d);
    lp1->kmer_prg.add_node(p);

    d = {Interval(0,0)};
    p.initialize(d);
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    kh = hash.kmerhash("AGC",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    kh = hash.kmerhash("AGT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kh = hash.kmerhash("ATT",3);
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

    d = {Interval(12,13), Interval(16,16), Interval(23,25)};
    p.initialize(d);
    kh = hash.kmerhash("TTA",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(23,26)};
    p.initialize(d);
    kh = hash.kmerhash("TAT",3);//inconsistent but I don't care
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

    d = {Interval(8,8)};
    p.initialize(d);
    lp4->kmer_prg.add_node(p);

    d = {Interval(8,11)};
    p.initialize(d);
    kh = hash.kmerhash("CTA",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp4->kmer_prg.add_node(p);

    d = {Interval(9,12)};
    p.initialize(d);
    kh = hash.kmerhash("TAG",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp4->kmer_prg.add_node(p);

    d = {Interval(12,12)};
    p.initialize(d);
    lp4->kmer_prg.add_node(p);

    d = {Interval(0,0)};
    p.initialize(d);
    lp2->kmer_prg.add_node(p);

    d = {Interval(0,3)};
    p.initialize(d);
    kh = hash.kmerhash("CTA",3);
    idx->add_record(min(kh.first,kh.second), 2, p, (kh.first < kh.second));
    lp2->kmer_prg.add_node(p);

    d = {Interval(1,4)};
    p.initialize(d);
    kh = hash.kmerhash("TAC",3);
    idx->add_record(min(kh.first,kh.second), 2, p, (kh.first < kh.second));
    lp2->kmer_prg.add_node(p);

    d = {Interval(2,5)};
    p.initialize(d);
    kh = hash.kmerhash("ACT",3);
    idx->add_record(min(kh.first,kh.second), 2, p, (kh.first < kh.second));
    lp2->kmer_prg.add_node(p);

    d = {Interval(5,5)};
    p.initialize(d);
    lp2->kmer_prg.add_node(p);

    // add read hits to mhs
    Seq *s;
    s = new Seq(0, "read2", "AGTTATGCTAGCTACTTACGGTA", 1, 3);
    add_read_hits(s, mhs, idx);
    delete s;
    //add_read_hits(0, "read2", "AGTTATGCTAGCTACTTACGGTA", mhs, idx, 1, 3);

    // initialize pangraph;
    PanGraph *pg;
    pg = new PanGraph();
    infer_localPRG_order_for_reads(prgs, mhs, pg, 1, 3);

    // create a pangraph object representing the truth we expect (prg 3 4 2 1)
    // note that prgs 1, 3, 4 share no 3mer, but 2 shares a 3mer with each of 2 other prgs
    PanGraph pg_exp;
    MinimizerHits mhs_dummy;
    pg_exp.add_node(1,0, mhs_dummy.hits);
    pg_exp.add_node(2,0, mhs_dummy.hits);
    pg_exp.add_node(3,0, mhs_dummy.hits);
    pg_exp.add_node(0,0, mhs_dummy.hits);
    pg_exp.add_edge(3,0);
    pg_exp.add_edge(0,2);
    pg_exp.add_edge(2,1);

    EXPECT_EQ(pg_exp, *pg);
    delete pg;
    delete lp1;
    delete lp2;
    delete lp3;
    delete lp4;
    delete mhs;
    idx->clear();
    delete idx;
}

TEST_F(UtilsTest, pangraphFromReadFile)
{
    MinimizerHits* mhs;
    mhs = new MinimizerHits();
    KmerHash hash;

    // initialize a prgs object
    vector<LocalPRG*> prgs;
    LocalPRG* lp1;
    LocalPRG* lp2;
    LocalPRG* lp3;
    LocalPRG* lp4;
    lp1 = new LocalPRG(1, "1", "");
    lp3 = new LocalPRG(3, "3", "");
    lp4 = new LocalPRG(0, "0", "");
    lp2 = new LocalPRG(2, "2", "");
    prgs.push_back(lp4);
    prgs.push_back(lp1);
    prgs.push_back(lp2);
    prgs.push_back(lp3);

    // initialize index as we would expect with example prgs
    Index *idx;
    idx = new Index();

    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    lp1->kmer_prg.add_node(p);
    d = {Interval(0,3)};
    p.initialize(d);
    pair<uint64_t,uint64_t> kh = hash.kmerhash("TAC",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(1,4)};
    p.initialize(d);
    kh = hash.kmerhash("ACG",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(2,5)};
    p.initialize(d);
    kh = hash.kmerhash("CGG",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(3,6)};
    p.initialize(d);
    kh = hash.kmerhash("GGT",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(4,7)};
    p.initialize(d);
    kh = hash.kmerhash("GTA",3);
    idx->add_record(min(kh.first,kh.second), 1, p, (kh.first < kh.second));
    lp1->kmer_prg.add_node(p);

    d = {Interval(7,7)};
    p.initialize(d);
    lp1->kmer_prg.add_node(p);

    d = {Interval(0,0)};
    p.initialize(d);
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    kh = hash.kmerhash("AGC",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    kh = hash.kmerhash("AGT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kh = hash.kmerhash("ATT",3);
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

    d = {Interval(12,13), Interval(16,16), Interval(23,25)};
    p.initialize(d);
    kh = hash.kmerhash("TTA",3);
    idx->add_record(min(kh.first,kh.second), 3, p, (kh.first < kh.second));
    lp3->kmer_prg.add_node(p);

    d = {Interval(23,26)};
    p.initialize(d);
    kh = hash.kmerhash("TAT",3);//inconsistent but I don't care
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

    d = {Interval(8,8)};
    p.initialize(d);
    lp4->kmer_prg.add_node(p);

    d = {Interval(8,11)};
    p.initialize(d);
    kh = hash.kmerhash("CTA",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp4->kmer_prg.add_node(p);

    d = {Interval(9,12)};
    p.initialize(d);
    kh = hash.kmerhash("TAG",3);
    idx->add_record(min(kh.first,kh.second), 0, p, (kh.first < kh.second));
    lp4->kmer_prg.add_node(p);

    d = {Interval(12,12)};
    p.initialize(d);
    lp4->kmer_prg.add_node(p);

    d = {Interval(0,0)};
    p.initialize(d);
    lp2->kmer_prg.add_node(p);

    d = {Interval(0,3)};
    p.initialize(d);
    kh = hash.kmerhash("CTA",3);
    idx->add_record(min(kh.first,kh.second), 2, p, (kh.first < kh.second));
    lp2->kmer_prg.add_node(p);

    d = {Interval(1,4)};
    p.initialize(d);
    kh = hash.kmerhash("TAC",3);
    idx->add_record(min(kh.first,kh.second), 2, p, (kh.first < kh.second));
    lp2->kmer_prg.add_node(p);

    d = {Interval(2,5)};
    p.initialize(d);
    kh = hash.kmerhash("ACT",3);
    idx->add_record(min(kh.first,kh.second), 2, p, (kh.first < kh.second));
    lp2->kmer_prg.add_node(p);

    d = {Interval(5,5)};
    p.initialize(d);
    lp2->kmer_prg.add_node(p);

    // initialize pangraph;
    PanGraph *pg;
    pg = new PanGraph();
    pangraph_from_read_file("../test/test_cases/read2.fa", mhs, pg, idx, prgs, 1, 3, 1);

    // create a pangraph object representing the truth we expect (prg 3 4 2 1)
    // note that prgs 1, 3, 4 share no 3mer, but 2 shares a 3mer with each of 2 other prgs
    PanGraph pg_exp;
    MinimizerHits mhs_dummy;
    pg_exp.add_node(1,0, mhs_dummy.hits);
    pg_exp.add_node(2,0, mhs_dummy.hits);
    pg_exp.add_node(3,0, mhs_dummy.hits);
    pg_exp.add_node(0,0, mhs_dummy.hits);
    pg_exp.add_edge(3,0);
    pg_exp.add_edge(0,2);
    pg_exp.add_edge(2,1);

    EXPECT_EQ(pg_exp, *pg);
    delete pg;
    delete mhs;
    delete lp1;
    delete lp2;
    delete lp3;
    delete lp4;
    idx->clear();
    delete idx;
}

//update_covgs_from_hits
//p_null

