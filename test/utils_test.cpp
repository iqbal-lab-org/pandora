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
    EXPECT_EQ(idx->minhash.size(), j);
    j = 1;
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
    j = 3;
    EXPECT_EQ(idx->minhash.size(), j);
    j = 2;
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
    j = 4;
    EXPECT_EQ(idx->minhash.size(), j);
    j = 3;
    EXPECT_EQ(prgs.size(), j);
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
    j = 4;
    EXPECT_EQ(idx->minhash.size(), j);
    j = 3;
    EXPECT_EQ(prgs.size(), j);
    uint64_t kh = kmerhash("AGC", 3);
    EXPECT_EQ(idx->minhash[kh].size(), j);
    kh = kmerhash("GCT", 3);
    EXPECT_EQ(idx->minhash[kh].size(), j);
    kh = kmerhash("AGT", 3);
    EXPECT_EQ(idx->minhash[kh].size(), j);
    j = 1;
    kh = kmerhash("GTT", 3);
    EXPECT_EQ(idx->minhash[kh].size(), j);
    delete idx;
}

TEST_F(UtilsTest, addReadHits){
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    MinimizerHits expected1;
    MinimizerHits expected2;
    MinimizerHit *m1, *m2, *m3, *m4;

    // initialize index as we would expect with example prgs 1 and 3 from above
    Index *idx;
    idx = new Index();
    deque<Interval> d = {Interval(0,3)};
    Path p;
    p.initialize(d);
    idx->add_record(kmerhash("AGC",3), 1, p);
    m1 = new MinimizerHit(0, Interval(0,3), 1, p, 1);
    expected1.hits.insert(m1);
    d = {Interval(1,4)};
    p.initialize(d);
    idx->add_record(kmerhash("GCT",3), 1, p);
    m2 = new MinimizerHit(0, Interval(1,4), 1, p, 1);
    expected2.hits.insert(m2);
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    idx->add_record(kmerhash("AGC",3), 3, p);
    m3 = new MinimizerHit(0, Interval(0,3), 3, p, 1);
    expected1.hits.insert(m3);
    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    idx->add_record(kmerhash("AGT",3), 3, p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("AGT",3), 3, p);
    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("GCT",3), 3, p);
    m4 = new MinimizerHit(0, Interval(1,4), 3, p, 1);
    expected2.hits.insert(m4);
    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("GTT",3), 3, p);

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
    idx->add_record(kmerhash("TAC",3), 1, p);

    d = {Interval(1,4)};
    p.initialize(d);
    idx->add_record(kmerhash("ACG",3), 1, p);

    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    idx->add_record(kmerhash("AGC",3), 3, p);

    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    idx->add_record(kmerhash("AGT",3), 3, p);

    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("AGT",3), 3, p);

    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("GCT",3), 3, p);

    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("GTT",3), 3, p);

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
    idx->add_record(kmerhash("TAC",3), 1, p);

    d = {Interval(1,4)};
    p.initialize(d);
    idx->add_record(kmerhash("ACG",3), 1, p);

    d = {Interval(2,5)};
    p.initialize(d);
    idx->add_record(kmerhash("CGG",3), 1, p);

    d = {Interval(3,6)};
    p.initialize(d);
    idx->add_record(kmerhash("GGT",3), 1, p);

    d = {Interval(4,7)};
    p.initialize(d);
    idx->add_record(kmerhash("GTA",3), 1, p);

    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    idx->add_record(kmerhash("AGC",3), 3, p);

    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    idx->add_record(kmerhash("AGT",3), 3, p);

    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("AGT",3), 3, p);

    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("GCT",3), 3, p);

    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("GTT",3), 3, p);

    d = {Interval(12,13), Interval(16,16), Interval(23,25)};
    p.initialize(d);
    idx->add_record(kmerhash("TTA",3), 3, p);

    d = {Interval(23,26)};
    p.initialize(d);
    idx->add_record(kmerhash("TAA",3), 3, p);

    d = {Interval(24,27)};
    p.initialize(d);
    idx->add_record(kmerhash("AAG",3), 3, p);

    d = {Interval(8,11)};
    p.initialize(d);
    idx->add_record(kmerhash("CTA",3), 4, p);

    d = {Interval(9,12)};
    p.initialize(d);
    idx->add_record(kmerhash("TAG",3), 4, p);

    d = {Interval(0,3)};
    p.initialize(d);
    idx->add_record(kmerhash("CTA",3), 2, p);

    d = {Interval(1,4)};
    p.initialize(d);
    idx->add_record(kmerhash("TAC",3), 2, p);

    d = {Interval(2,5)};
    p.initialize(d);
    idx->add_record(kmerhash("ACT",3), 2, p);

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
    // should give exactly the same results, but read the read from a file

    // initialize index as we would expect with example prgs
    Index *idx;
    idx = new Index();

    deque<Interval> d = {Interval(0,3)};
    Path p;
    p.initialize(d);
    idx->add_record(kmerhash("TAC",3), 1, p);

    d = {Interval(1,4)};
    p.initialize(d);
    idx->add_record(kmerhash("ACG",3), 1, p);

    d = {Interval(2,5)};
    p.initialize(d);
    idx->add_record(kmerhash("CGG",3), 1, p);

    d = {Interval(3,6)};
    p.initialize(d);
    idx->add_record(kmerhash("GGT",3), 1, p);

    d = {Interval(4,7)};
    p.initialize(d);
    idx->add_record(kmerhash("GTA",3), 1, p);

    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    idx->add_record(kmerhash("AGC",3), 3, p);

    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    idx->add_record(kmerhash("AGT",3), 3, p);

    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("AGT",3), 3, p);

    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("GCT",3), 3, p);

    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    idx->add_record(kmerhash("GTT",3), 3, p);

    d = {Interval(12,13), Interval(16,16), Interval(23,25)};
    p.initialize(d);
    idx->add_record(kmerhash("TTA",3), 3, p);

    d = {Interval(23,26)};
    p.initialize(d);
    idx->add_record(kmerhash("TAA",3), 3, p);

    d = {Interval(24,27)};
    p.initialize(d);
    idx->add_record(kmerhash("AAG",3), 3, p);

    d = {Interval(8,11)};
    p.initialize(d);
    idx->add_record(kmerhash("CTA",3), 4, p);

    d = {Interval(9,12)};
    p.initialize(d);
    idx->add_record(kmerhash("TAG",3), 4, p);

    d = {Interval(0,3)};
    p.initialize(d);
    idx->add_record(kmerhash("CTA",3), 2, p);

    d = {Interval(1,4)};
    p.initialize(d);
    idx->add_record(kmerhash("TAC",3), 2, p);

    d = {Interval(2,5)};
    p.initialize(d);
    idx->add_record(kmerhash("ACT",3), 2, p);

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
}
