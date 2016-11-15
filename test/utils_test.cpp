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

TEST_F(UtilsTest, indexPrgFile){
    vector<LocalPRG*> prgs;
    Index *idx;
    idx = new Index;

    // simple case first, single prg with empty string sequence
    // doesn't get added to prgs 
    index_prg_file(prgs, "../test/test_cases/prg0.fa", idx, 1,3);
    uint32_t j = 0;
    EXPECT_EQ(idx->minhash.size(), j);
    EXPECT_EQ(prgs.size(), j);
 
    // single prg with simple sequence
    index_prg_file(prgs, "../test/test_cases/prg1.fa", idx, 1,3);
    LocalPRG l1(1,"prg1", "AGCT");
    j = 2;
    pair<uint64_t,uint64_t> kh = kmerhash("AGC",3);
    EXPECT_EQ(idx->minhash[min(kh.first, kh.second)].size(), j);
    j = 1;
    EXPECT_EQ(idx->minhash.size(), j);
    EXPECT_EQ(prgs.size(), j);
    j = 0;
    EXPECT_EQ(prgs[0]->id, j);
    EXPECT_EQ(prgs[0]->name, "prg1");
    EXPECT_EQ(prgs[0]->seq, "AGCT");
    EXPECT_EQ(prgs[0]->prg, l1.prg);

    // single prg with a variant site
    idx->clear();
    j = 0;
    EXPECT_EQ(idx->minhash.size(), j);
    index_prg_file(prgs, "../test/test_cases/prg2.fa", idx, 1,3);
    LocalPRG l2(2,"prg2", "A 5 GC 6 G 5 T");
    j = 2;
    EXPECT_EQ(idx->minhash.size(), j);
    EXPECT_EQ(idx->minhash[min(kh.first,kh.second)].size(), j);
    EXPECT_EQ(prgs.size(), j);
    j = 0;
    EXPECT_EQ(prgs[1]->id, j);
    EXPECT_EQ(prgs[1]->name, "prg2");
    EXPECT_EQ(prgs[1]->seq, "A 5 GC 6 G 5 T");
    EXPECT_EQ(prgs[1]->prg, l2.prg);

    // single prg with a nested variant site
    idx->clear();
    j = 0;
    EXPECT_EQ(idx->minhash.size(), j);
    index_prg_file(prgs, "../test/test_cases/prg3.fa", idx, 1,3);
    LocalPRG l3 = LocalPRG(3,"prg3", "A 5 G 7 C 8 T 7  6 G 5 T");
    j = 3;
    EXPECT_EQ(idx->minhash.size(), j);
    EXPECT_EQ(prgs.size(), j);
    j = 2;
    EXPECT_EQ(idx->minhash[min(kh.first,kh.second)].size(), j);
    kh = kmerhash("AGT",3);
    EXPECT_EQ(idx->minhash[min(kh.first,kh.second)].size(), j);
    j = 0;
    EXPECT_EQ(prgs[2]->id, j);
    EXPECT_EQ(prgs[2]->name, "prg3");
    EXPECT_EQ(prgs[2]->seq, "A 5 G 7 C 8 T 7  6 G 5 T");
    EXPECT_EQ(prgs[2]->prg, l3.prg);    

    // now a prg input file with all 4 in
    idx->clear();
    j = 0;
    EXPECT_EQ(idx->minhash.size(), j);
    prgs.clear();
    EXPECT_EQ(prgs.size(), j);
    index_prg_file(prgs, "../test/test_cases/prg0123.fa", idx, 1,3);
    j = 3;
    EXPECT_EQ(idx->minhash.size(), j);
    EXPECT_EQ(prgs.size(), j);
    j = 6;
    kh = kmerhash("AGC",3);
    EXPECT_EQ(idx->minhash[min(kh.first,kh.second)].size(), j);
    kh = kmerhash("GCT",3);
    EXPECT_EQ(idx->minhash[min(kh.first,kh.second)].size(), j); // same
    j = 3;
    kh = kmerhash("AGT",3);
    EXPECT_EQ(idx->minhash[min(kh.first,kh.second)].size(), j);
    kh = kmerhash("ACT",3);
    EXPECT_EQ(idx->minhash[min(kh.first,kh.second)].size(), j); // same
    j = 1;
    kh = kmerhash("GTT",3);
    EXPECT_EQ(idx->minhash[min(kh.first,kh.second)].size(), j); //same
    kh = kmerhash("AAC",3);
    EXPECT_EQ(idx->minhash[min(kh.first,kh.second)].size(), j);
    delete idx;
}

TEST_F(UtilsTest, addReadHits){
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    MinimizerHits expected1;
    MinimizerHits expected2;
    MinimizerHit *m1, *m2, *m3, *m4, *m5, *m6, *m7, *m8;

    // initialize index as we would expect with example prgs 1 and 3 from above
    Index *idx;
    idx = new Index();
    deque<Interval> d = {Interval(0,3)};
    Path p;
    p.initialize(d);
    pair<uint64_t,uint64_t> kh = kmerhash("AGC",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 0);
    m1 = new MinimizerHit(0, Interval(0,3), 1, p, 0);
    m2 = new MinimizerHit(0, Interval(1,4), 1, p, 1);
    expected1.hits.insert(m1);
    expected2.hits.insert(m2);
    d = {Interval(1,4)};
    p.initialize(d);
    kh = kmerhash("GCT",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 1);
    m3 = new MinimizerHit(0, Interval(1,4), 1, p, 0);
    m4 = new MinimizerHit(0, Interval(0,3), 1, p, 1);
    expected2.hits.insert(m3);
    expected1.hits.insert(m4);
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    kh = kmerhash("AGC",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);
    m5 = new MinimizerHit(0, Interval(0,3), 3, p, 0);
    m6 = new MinimizerHit(0, Interval(1,4), 3, p, 1);
    expected1.hits.insert(m5);
    expected2.hits.insert(m6);
    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    kh = kmerhash("AGT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);
    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = kmerhash("GCT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);
    m7 = new MinimizerHit(0, Interval(1,4), 3, p, 0);
    m8 = new MinimizerHit(0, Interval(0,3), 3, p, 1);
    expected2.hits.insert(m7);
    expected1.hits.insert(m8);
    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = kmerhash("GTT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    add_read_hits(0, "read1", "AGC", mhs, idx, 1, 3);
    set<MinimizerHit*, pComp>::const_iterator it2 = expected1.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }
    
    // same if take w=2 as sketch of read AGCT should just contain AGC
    delete mhs;
    mhs = new MinimizerHits();
    add_read_hits(0, "read2", "AGCT", mhs, idx, 2, 3);
    it2 = expected1.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    // now back to w = 1, add expected2 to expected1 as will get hits against both AGC and GCT
    delete mhs;
    mhs = new MinimizerHits();
    add_read_hits(0, "read3", "AGCT", mhs, idx, 1, 3);
    expected1.hits.insert(expected2.hits.begin(), expected2.hits.end());
    EXPECT_EQ(expected1.hits.size(), mhs->hits.size());
    it2 = expected1.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }
    expected1.hits.clear();
    expected2.hits.clear();
    delete idx;
    delete m1;
    delete m2;
    delete m3;
    delete m4;
    delete m5;
    delete m6;
    delete m7;
    delete m8;
    delete mhs;
}

TEST_F(UtilsTest, simpleInferLocalPRGOrderForRead){    
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();

    // initialize index as we would expect with example prgs (variant of) 1 and 3 from above
    Index *idx;
    idx = new Index();

    deque<Interval> d = {Interval(0,3)};
    Path p;
    p.initialize(d);
    pair<uint64_t,uint64_t> kh = kmerhash("TAC",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 1);

    d = {Interval(1,4)};
    p.initialize(d);
    kh = kmerhash("ACG",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 1);

    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    kh = kmerhash("AGC",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    kh = kmerhash("AGT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = kmerhash("GCT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = kmerhash("GTT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    // add read hits to mhs
    add_read_hits(0, "read1", "AGTTTACG", mhs, idx, 1, 3);

    // initialize pangraph;
    PanGraph *pg;
    pg = new PanGraph();
    infer_localPRG_order_for_read(mhs, pg, 1, 1, 3);

    // create a pangraph object representing the truth we expect (prg 3 then 1)
    PanGraph pg_exp;
    pg_exp.add_node(1,0);
    pg_exp.add_node(3,0);
    pg_exp.add_edge(3,1);

    EXPECT_EQ(pg_exp, *pg);
    delete idx;
    delete pg;
    delete mhs;
}

TEST_F(UtilsTest, biggerInferLocalPRGOrderForRead){
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();

    // initialize index as we would expect with example prgs
    Index *idx;
    idx = new Index();

    deque<Interval> d = {Interval(0,3)};
    Path p;
    p.initialize(d);
    pair<uint64_t,uint64_t> kh = kmerhash("TAC",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 0);

    d = {Interval(1,4)};
    p.initialize(d);
    kh = kmerhash("ACG",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 1);

    d = {Interval(2,5)};
    p.initialize(d);
    kh = kmerhash("CGG",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 1);

    d = {Interval(3,6)};
    p.initialize(d);
    kh = kmerhash("GGT",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 1);

    d = {Interval(4,7)};
    p.initialize(d);
    kh = kmerhash("GTA",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 0);

    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    kh = kmerhash("AGC",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    kh = kmerhash("AGT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = kmerhash("GCT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = kmerhash("GTT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(12,13), Interval(16,16), Interval(23,25)};
    p.initialize(d);
    kh = kmerhash("TTA",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(23,26)};
    p.initialize(d);
    kh = kmerhash("TAA",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    d = {Interval(24,27)};
    p.initialize(d);
    kh = kmerhash("AAG",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    d = {Interval(8,11)};
    p.initialize(d);
    kh = kmerhash("CTA",3);
    idx->add_record(min(kh.first,kh.second), 4, p, 0);

    d = {Interval(9,12)};
    p.initialize(d);
    kh = kmerhash("TAG",3);
    idx->add_record(min(kh.first,kh.second), 4, p, 0);

    d = {Interval(0,3)};
    p.initialize(d);
    kh = kmerhash("CTA",3);
    idx->add_record(min(kh.first,kh.second), 2, p, 0);

    d = {Interval(1,4)};
    p.initialize(d);
    kh = kmerhash("TAC",3);
    idx->add_record(min(kh.first,kh.second), 2, p, 0);

    d = {Interval(2,5)};
    p.initialize(d);
    kh = kmerhash("ACT",3);
    idx->add_record(min(kh.first,kh.second), 2, p, 0);

    // add read hits to mhs
    add_read_hits(0, "read2", "AGTTAAGCTAGCTACTTACGGTA", mhs, idx, 1, 3);

    // initialize pangraph;
    PanGraph *pg;
    pg = new PanGraph();
    infer_localPRG_order_for_read(mhs, pg, 1, 1, 3);

    // create a pangraph object representing the truth we expect (prg 3 4 2 1)
    // note that prgs 1, 3, 4 share no 3mer, but 2 shares a 3mer with each of 2 other prgs
    PanGraph pg_exp;
    pg_exp.add_node(1,0);
    pg_exp.add_node(2,0);
    pg_exp.add_node(3,0);
    pg_exp.add_node(4,0);
    pg_exp.add_edge(3,4);
    pg_exp.add_edge(4,2);
    pg_exp.add_edge(2,1);

    EXPECT_EQ(pg_exp, *pg);
    delete pg;
    delete mhs;
    delete idx;
}

TEST_F(UtilsTest, pangraphFromReadFile){
    vector<LocalPRG*> prgs;
    LocalPRG *lp0, *lp1, *lp2, *lp3, *lp4;
    lp0 = new LocalPRG(0,"0","");
    prgs.push_back(lp0);
    lp1 = new LocalPRG(1,"1","");
    prgs.push_back(lp1);
    lp2 = new LocalPRG(2,"2","");
    prgs.push_back(lp2);
    lp3 = new LocalPRG(3,"3", "");
    prgs.push_back(lp3);
    lp4 = new LocalPRG(4,"4", "");
    prgs.push_back(lp3);
    // should give exactly the same results, but read the read from a file

    // initialize index as we would expect with example prgs
    Index *idx;
    idx = new Index();

    deque<Interval> d = {Interval(0,3)};
    Path p;
    p.initialize(d);
    pair<uint64_t,uint64_t> kh = kmerhash("ACG",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 1);

    d = {Interval(2,5)};
    p.initialize(d);
    kh = kmerhash("CGG",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 1);

    d = {Interval(3,6)};
    p.initialize(d);
    kh = kmerhash("GGT",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 1);

    d = {Interval(4,7)};
    p.initialize(d);
    kh = kmerhash("GTA",3);
    idx->add_record(min(kh.first,kh.second), 1, p, 0);

    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    kh = kmerhash("AGC",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    kh = kmerhash("AGT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = kmerhash("GCT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kh = kmerhash("GTT",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(12,13), Interval(16,16), Interval(23,25)};
    p.initialize(d);
    kh = kmerhash("TTA",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 1);

    d = {Interval(23,26)};
    p.initialize(d);
    kh = kmerhash("TAA",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    d = {Interval(24,27)};
    p.initialize(d);
    kh = kmerhash("AAG",3);
    idx->add_record(min(kh.first,kh.second), 3, p, 0);

    d = {Interval(8,11)};
    p.initialize(d);
    kh = kmerhash("CTA",3);
    idx->add_record(min(kh.first,kh.second), 4, p, 0);

    d = {Interval(9,12)};
    p.initialize(d);
    kh = kmerhash("TAG",3);
    idx->add_record(min(kh.first,kh.second), 4, p, 0);

    d = {Interval(0,3)};
    p.initialize(d);
    kh = kmerhash("CTA",3);
    idx->add_record(min(kh.first,kh.second), 2, p, 0);

    d = {Interval(1,4)};
    p.initialize(d);
    kh = kmerhash("TAC",3);
    idx->add_record(min(kh.first,kh.second), 2, p, 0);

    d = {Interval(2,5)};
    p.initialize(d);
    kh = kmerhash("ACT",3);
    idx->add_record(min(kh.first,kh.second), 2, p, 0);

    // initialize pangraph;
    PanGraph *pg;
    pg = new PanGraph();

    pangraph_from_read_file("../test/test_cases/read0.fa", pg, idx, prgs, 1, 3, 1, 1);

    // create a pangraph object representing the truth we expect (prg 3 4 2 1)
    // note that prgs 1, 3, 4 share no 3mer, but 2 shares a 3mer with each of 2 other prgs
    PanGraph pg_exp;
    pg_exp.add_node(1,0);
    pg_exp.add_node(2,0);
    pg_exp.add_node(3,0);
    pg_exp.add_node(4,0);
    pg_exp.add_edge(3,4);
    pg_exp.add_edge(4,2);
    pg_exp.add_edge(2,1);

    EXPECT_EQ(pg_exp, *pg);
    delete idx;
    delete pg;
    delete lp0;
    delete lp1;
    delete lp2;
    delete lp3;
    delete lp4;
}
