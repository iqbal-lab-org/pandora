#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "utils.h"
#include "localPRG.h"
#include "interval.h"
#include "path.h"
#include "minihit.h"
#include "minihits.h"
#include "index.h"
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
    LocalPRG l1 = LocalPRG(1,"prg1", "AGCT");
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
    LocalPRG l2 = LocalPRG(2,"prg2", "A 5 GC 6 G 5 T");
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
    EXPECT_EQ(idx->minhash["AGC"].size(), j);
    EXPECT_EQ(idx->minhash["GCT"].size(), j);
    EXPECT_EQ(idx->minhash["AGT"].size(), j);
    j = 1;
    EXPECT_EQ(idx->minhash["GTT"].size(), j);
    delete idx;
}

TEST_F(UtilsTest, addReadHits){
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    MinimizerHits expected1 = MinimizerHits();
    MinimizerHits expected2 = MinimizerHits();
    MinimizerHit *m1, *m2, *m3, *m4;

    // initialize index as we would expect with example prgs 1 and 3 from above
    Index *idx;
    idx = new Index();
    deque<Interval> d = {Interval(0,3)};
    Path p = Path();
    p.initialize(d);
    idx->add_record("AGC", 1, p);
    m1 = new MinimizerHit(0, Interval(0,3), 1, p, 1);
    expected1.hits.insert(m1);
    d = {Interval(1,4)};
    p.initialize(d);
    idx->add_record("GCT", 1, p);
    m2 = new MinimizerHit(0, Interval(1,4), 1, p, 1);
    expected2.hits.insert(m2);
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    p.initialize(d);
    idx->add_record("AGC", 3, p);
    m3 = new MinimizerHit(0, Interval(0,3), 3, p, 1);
    expected1.hits.insert(m3);
    d = {Interval(0,1), Interval(4,5), Interval(12,13)};
    p.initialize(d);
    idx->add_record("AGT", 3, p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    idx->add_record("AGT", 3, p);
    d = {Interval(4,5), Interval(8,9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    idx->add_record("GCT", 3, p);
    m4 = new MinimizerHit(0, Interval(1,4), 3, p, 1);
    expected2.hits.insert(m4);
    d = {Interval(4,5), Interval(12,13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    idx->add_record("GTT", 3, p);

    cout << "add read1 hits" << endl;
    add_read_hits(0, "read1", "AGC", mhs, idx, 1, 3);
    cout << "check if correct" << endl;
    set<MinimizerHit*, pComp>::const_iterator it2 = expected1.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }
    
    // same if take w=2 as sketch of read AGCT should just contain AGC
    cout << "add read2 hits" << endl;
    delete mhs;
    mhs = new MinimizerHits();
    cout << "just cleared previous hits" << endl;
    add_read_hits(0, "read2", "AGCT", mhs, idx, 2, 3);
    cout << "check if correct" << endl;
    it2 = expected1.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    // now back to w = 1, add expected2 to expected1 as will get hits against both AGC and GCT
    cout << "clear previous hits" << endl;
    delete mhs;
    mhs = new MinimizerHits();
    cout << "add read3 hits" << endl;
    add_read_hits(0, "read3", "AGCT", mhs, idx, 1, 3);
    cout << "extend set with other set" << endl;
    expected1.hits.insert(expected2.hits.begin(), expected2.hits.end());
    cout << "check if correct" << endl;
    it2 = expected1.hits.begin();
    for (set<MinimizerHit*, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it)
    {
        EXPECT_EQ(**it2, **it);
        it2++;
    }
    cout << "delete stuff" << endl;
    delete idx;
    cout << "index deleted" << endl;
    delete m1, m2, m3, m4;
    cout << "hit pointers deleted" << endl;
    delete mhs;
    cout << "minihits poiter deleted" << endl;
}

TEST_F(UtilsTest, inferLocalPRGOrderForRead){
}

TEST_F(UtilsTest, pangraphFromReadFile){
}
