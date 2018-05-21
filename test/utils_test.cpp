#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "utils.h"
#include "localPRG.h"
#include "pangenome/pangraph.h"
#include "interval.h"
#include "path.h"
#include "minihit.h"
#include "minihits.h"
#include "index.h"
#include "inthash.h"
#include "seq.h"
#include <stdint.h>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

TEST(UtilsTest, split) {
    vector<string> v = {"abc", "def", "ghi"};
    EXPECT_EQ(v, split("abc, def, ghi", ", "));
    EXPECT_EQ(v, split("abc, def, ghi, ", ", "));
    EXPECT_EQ(v, split(", abc, def, ghi", ", "));
}

TEST(UtilsTest, revComplement) {
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

TEST(UtilsTest, readPrgFile) {
    vector<LocalPRG *> prgs;

    // simple case first, single prg with empty string sequence
    // doesn't get added to prgs 
    read_prg_file(prgs, "../../test/test_cases/prg0.fa");
    uint32_t j = 0;
    EXPECT_EQ(prgs.size(), j);

    // single prg with simple sequence
    read_prg_file(prgs, "../../test/test_cases/prg1.fa");
    LocalPRG l1(1, "prg1", "AGCT");
    j = 1;
    EXPECT_EQ(prgs.size(), j);
    j = 0;
    EXPECT_EQ(prgs[0]->id, j);
    EXPECT_EQ(prgs[0]->name, "prg1");
    EXPECT_EQ(prgs[0]->seq, "AGCT");
    EXPECT_EQ(prgs[0]->prg, l1.prg);

    // single prg with a variant site
    read_prg_file(prgs, "../../test/test_cases/prg2.fa");
    LocalPRG l2(2, "prg2", "A 5 GC 6 G 5 T");
    j = 2;
    EXPECT_EQ(prgs.size(), j);
    j = 0;
    EXPECT_EQ(prgs[1]->id, j);
    EXPECT_EQ(prgs[1]->name, "prg2");
    EXPECT_EQ(prgs[1]->seq, "A 5 GC 6 G 5 T");
    EXPECT_EQ(prgs[1]->prg, l2.prg);

    // single prg with a nested variant site
    read_prg_file(prgs, "../../test/test_cases/prg3.fa");
    LocalPRG l3 = LocalPRG(3, "prg3", "A 5 G 7 C 8 T 7  6 G 5 T");
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
    read_prg_file(prgs, "../../test/test_cases/prg0123.fa");
    j = 3;
    EXPECT_EQ(prgs.size(), j);
}

TEST(UtilsTest, addReadHits) {
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    MinimizerHits expected1;
    MinimizerHits expected2;
    MinimizerHits expected3;
    MinimizerHits expected4;

    // initialize index as we would expect with example prgs 1 and 3 from above
    KmerHash hash;
    Index *idx;
    idx = new Index();
    deque<Interval> d = {Interval(0, 3)};
    Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("AGC", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    MinimizerHitPtr m1(make_shared<MinimizerHit>(0, Interval(0, 3), 1, p, 0, 1));
    MinimizerHitPtr m2(make_shared<MinimizerHit>(0, Interval(1, 4), 1, p, 0, 0));
    expected1.hits.insert(m1);
    expected2.hits.insert(m2);
    d = {Interval(1, 4)};
    p.initialize(d);
    kh = hash.kmerhash("GCT", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    MinimizerHitPtr m3(make_shared<MinimizerHit>(0, Interval(1, 4), 1, p, 0, 1));
    MinimizerHitPtr m4(make_shared<MinimizerHit>(0, Interval(0, 3), 1, p, 0, 0));
    expected2.hits.insert(m3);
    expected1.hits.insert(m4);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kh = hash.kmerhash("AGC", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    MinimizerHitPtr m5(make_shared<MinimizerHit>(0, Interval(0, 3), 3, p, 0, 1));
    MinimizerHitPtr m6(make_shared<MinimizerHit>(0, Interval(1, 4), 3, p, 0, 0));
    expected1.hits.insert(m5);
    expected2.hits.insert(m6);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kh = hash.kmerhash("AGT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    MinimizerHitPtr m9(make_shared<MinimizerHit>(0, Interval(0, 3), 3, p, 0, 1));
    expected3.hits.insert(m9);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    MinimizerHitPtr m10(make_shared<MinimizerHit>(0, Interval(0, 3), 3, p, 0, 1));
    expected3.hits.insert(m10);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("GCT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    MinimizerHitPtr m7(make_shared<MinimizerHit>(0, Interval(1, 4), 3, p, 0, 1));
    MinimizerHitPtr m8(make_shared<MinimizerHit>(0, Interval(0, 3), 3, p, 0, 0));
    expected2.hits.insert(m7);
    expected1.hits.insert(m8);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("GTT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    MinimizerHitPtr m11(make_shared<MinimizerHit>(0, Interval(1, 4), 3, p, 0, 1));
    expected4.hits.insert(m11);

    Seq *s;
    s = new Seq(0, "read1", "AGC", 1, 3);
    add_read_hits(s, mhs, idx);
    mhs->sort();
    EXPECT_EQ(expected1.hits.size(), mhs->hits.size());
    set<MinimizerHitPtr, pComp>::const_iterator it2 = expected1.hits.begin();
    for (set<MinimizerHitPtr, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it) {
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
    mhs->sort();
    EXPECT_EQ(expected4.hits.size(), mhs->hits.size());
    it2 = expected4.hits.begin();
    for (set<MinimizerHitPtr, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it) {
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
    mhs->sort();
    EXPECT_EQ(expected3.hits.size(), mhs->hits.size());
    it2 = expected3.hits.begin();
    for (set<MinimizerHitPtr, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it) {
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
    mhs->sort();
    expected1.hits.insert(expected2.hits.begin(), expected2.hits.end());
    EXPECT_EQ(expected1.hits.size(), mhs->hits.size());
    it2 = expected1.hits.begin();
    for (set<MinimizerHitPtr, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it) {
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
    mhs->sort();
    EXPECT_EQ(expected1.hits.size(), mhs->hits.size());
    it2 = expected1.hits.begin();
    for (set<MinimizerHitPtr, pComp>::const_iterator it = mhs->hits.begin(); it != mhs->hits.end(); ++it) {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    expected1.hits.clear();
    expected2.hits.clear();
    expected3.hits.clear();
    expected4.hits.clear();
    idx->clear();
    delete idx;
    delete mhs;
    delete s;
}

TEST(UtilsTest, filter_clusters2) {
    deque<Interval> d = {Interval(0, 10)};
    Path p;
    p.initialize(d);

    set<MinimizerHitPtr, pComp> s;
    set<set<MinimizerHitPtr, pComp>, clusterComp> ss, ss_exp;

    MinimizerHitPtr mh;
    for (uint i = 0; i != 6; ++i) {
        mh = make_shared<MinimizerHit>(1, Interval(i, i + 10), 0, p, 0, 0);
        s.insert(mh);
    }
    ss.insert(s);
    ss_exp.insert(s);
    s.clear();
    for (uint i = 5; i != 15; ++i) {
        mh = make_shared<MinimizerHit>(1, Interval(i, i + 10), 1, p, 0, 0);
        s.insert(mh);
    }
    ss.insert(s);
    ss_exp.insert(s);
    s.clear();
    for (uint i = 3; i != 7; ++i) {
        mh = make_shared<MinimizerHit>(1, Interval(i, i + 10), 2, p, 0, 0);
        s.insert(mh);
    }
    ss.insert(s);

    filter_clusters2(ss, 20);

    EXPECT_EQ(ss_exp.size(), ss.size());

}

TEST(UtilsTest, simpleInferLocalPRGOrderForRead) {
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    KmerHash hash;

    // initialize a prgs object
    vector<LocalPRG *> prgs;
    LocalPRG *lp1;
    LocalPRG *lp3;
    lp1 = new LocalPRG(1, "1", "");
    lp3 = new LocalPRG(0, "0", "");
    prgs.push_back(lp3);
    prgs.push_back(lp1);

    // initialize index as we would expect with example prgs (variant of) 1 and 3 from above
    Index *idx;
    idx = new Index();

    vector<KmerNodePtr> v;
    KmerNodePtr kn;

    deque<Interval> d = {Interval(0, 0)};
    Path p;
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);

    d = {Interval(0, 3)};
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("TAC", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[0], v[1]);

    d = {Interval(1, 4)};
    p.initialize(d);
    kh = hash.kmerhash("ACG", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[1], v[2]);

    d = {Interval(4, 4)};
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[2], v[3]);

    d = {Interval(0, 0)};
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);

    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kh = hash.kmerhash("AGC", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[4], v[5]);

    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kh = hash.kmerhash("AGT", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[4], v[6]);

    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("ATT", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[4], v[7]);

    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("GCT", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[5], v[8]);

    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("GTT", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[6], v[9]);

    d = {Interval(12, 13), Interval(16, 16), Interval(23, 25)};
    p.initialize(d);
    kh = hash.kmerhash("TTA", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[9], v[10]);

    d = {Interval(23, 26)};
    p.initialize(d);
    kh = hash.kmerhash("TAA", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[11]);
    lp3->kmer_prg.add_edge(v[8], v[11]);
    lp3->kmer_prg.add_edge(v[10], v[11]);

    d = {Interval(24, 27)};
    p.initialize(d);
    kh = hash.kmerhash("AAG", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[11], v[12]);

    d = {Interval(27, 27)};
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[12], v[13]);

    // add read hits to mhs
    Seq *s;
    s = new Seq(0, "read1", "AGTTAAGTACG", 1, 3);
    add_read_hits(s, mhs, idx);
    delete s;
    //add_read_hits(0, "read1", "AGTTAAGTACG", mhs, idx, 1, 3);

    // initialize pangraph;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    infer_localPRG_order_for_reads(prgs, mhs, pg, 1, 100, 0.1, 1);

    // create a pangraph object representing the truth we expect (prg 3 then 1)
    pangenome::Graph pg_exp;
    MinimizerHits mhs_dummy;
    pg_exp.add_node(1, "1", 0, mhs_dummy.hits);
    pg_exp.add_node(0, "0", 0, mhs_dummy.hits);
    //pg_exp.add_edge(0,1,3,0);

    EXPECT_EQ(pg_exp, *pg);
    idx->clear();
    delete idx;
    delete pg;
    delete lp1;
    delete lp3;
    delete mhs;
}

TEST(UtilsTest, biggerInferLocalPRGOrderForRead) {
    // initialize minihits container
    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    KmerHash hash;

    // initialize a prgs object
    vector<LocalPRG *> prgs;
    LocalPRG *lp1;
    LocalPRG *lp2;
    LocalPRG *lp3;
    LocalPRG *lp4;
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

    vector<KmerNodePtr> v;
    KmerNodePtr kn;

    deque<Interval> d = {Interval(0, 0)};
    Path p;
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);

    d = {Interval(0, 3)};
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("TAC", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[0], v[1]);

    d = {Interval(1, 4)};
    p.initialize(d);
    kh = hash.kmerhash("ACG", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[1], v[2]);

    d = {Interval(2, 5)};
    p.initialize(d);
    kh = hash.kmerhash("CGG", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[2], v[3]);

    d = {Interval(3, 6)};
    p.initialize(d);
    kh = hash.kmerhash("GGT", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[3], v[4]);

    d = {Interval(4, 7)};
    p.initialize(d);
    kh = hash.kmerhash("GTA", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[4], v[5]);

    d = {Interval(7, 7)};
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[5], v[6]);

    d = {Interval(0, 0)};
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);

    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kh = hash.kmerhash("ACC", 3); // inconsistent
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[8]);

    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kh = hash.kmerhash("AGT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[9]);

    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("ATT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[10]);

    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("GCT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[8], v[11]);

    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("GTT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[9], v[12]);

    d = {Interval(12, 13), Interval(16, 16), Interval(23, 25)};
    p.initialize(d);
    kh = hash.kmerhash("TTA", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[12], v[13]);

    d = {Interval(23, 26)};
    p.initialize(d);
    kh = hash.kmerhash("TAT", 3);//inconsistent but I don't care
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[10], v[14]);
    lp3->kmer_prg.add_edge(v[11], v[14]);
    lp3->kmer_prg.add_edge(v[13], v[14]);

    d = {Interval(24, 27)};
    p.initialize(d);
    kh = hash.kmerhash("ATG", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[14], v[15]);

    d = {Interval(27, 27)};
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[15], v[16]);

    d = {Interval(8, 8)};
    p.initialize(d);
    kn = lp4->kmer_prg.add_node(p);
    v.push_back(kn);

    d = {Interval(8, 11)};
    p.initialize(d);
    kh = hash.kmerhash("CTA", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp4->kmer_prg.add_node(p);
    v.push_back(kn);
    lp4->kmer_prg.add_edge(v[17], v[18]);

    d = {Interval(9, 12)};
    p.initialize(d);
    kh = hash.kmerhash("TAG", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp4->kmer_prg.add_node(p);
    v.push_back(kn);
    lp4->kmer_prg.add_edge(v[18], v[19]);

    d = {Interval(12, 12)};
    p.initialize(d);
    kn = lp4->kmer_prg.add_node(p);
    v.push_back(kn);
    lp4->kmer_prg.add_edge(v[19], v[20]);

    d = {Interval(0, 0)};
    p.initialize(d);
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);

    d = {Interval(0, 3)};
    p.initialize(d);
    kh = hash.kmerhash("CTA", 3);
    idx->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[21], v[22]);

    d = {Interval(1, 4)};
    p.initialize(d);
    kh = hash.kmerhash("TAC", 3);
    idx->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[22], v[23]);

    d = {Interval(2, 5)};
    p.initialize(d);
    kh = hash.kmerhash("ACT", 3);
    idx->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[23], v[24]);

    d = {Interval(5, 5)};
    p.initialize(d);
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[24], v[25]);

    // add read hits to mhs
    Seq *s;
    s = new Seq(0, "read2", "AGTTATGCTAGCTACTTACGGTA", 1, 3);
    add_read_hits(s, mhs, idx);
    delete s;
    //add_read_hits(0, "read2", "AGTTATGCTAGCTACTTACGGTA", mhs, idx, 1, 3);

    // initialize pangraph;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    infer_localPRG_order_for_reads(prgs, mhs, pg, 1, 100, 0.1, 1);

    // create a pangraph object representing the truth we expect (prg 3 4 2 1)
    // note that prgs 1, 3, 4 share no 3mer, but 2 shares a 3mer with each of 2 other prgs
    pangenome::Graph pg_exp;
    MinimizerHits mhs_dummy;
    pg_exp.add_node(1, "1", 0, mhs_dummy.hits);
    pg_exp.add_node(2, "2", 0, mhs_dummy.hits);
    pg_exp.add_node(3, "3", 0, mhs_dummy.hits);
    pg_exp.add_node(0, "0", 0, mhs_dummy.hits);
    //pg_exp.add_edge(3,0,3,0);
    //pg_exp.add_edge(0,2,3,0);
    //pg_exp.add_edge(2,1,3,0);

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

TEST(UtilsTest, pangraphFromReadFile) {
    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    KmerHash hash;

    // initialize a prgs object
    vector<LocalPRG *> prgs;
    LocalPRG *lp1;
    LocalPRG *lp2;
    LocalPRG *lp3;
    LocalPRG *lp4;
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

    vector<KmerNodePtr> v;
    KmerNodePtr kn;

    deque<Interval> d = {Interval(0, 0)};
    Path p;
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);

    d = {Interval(0, 3)};
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("TAC", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[0], v[1]);

    d = {Interval(1, 4)};
    p.initialize(d);
    kh = hash.kmerhash("ACG", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[1], v[2]);

    d = {Interval(2, 5)};
    p.initialize(d);
    kh = hash.kmerhash("CGG", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[2], v[3]);

    d = {Interval(3, 6)};
    p.initialize(d);
    kh = hash.kmerhash("GGT", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[3], v[4]);

    d = {Interval(4, 7)};
    p.initialize(d);
    kh = hash.kmerhash("GTA", 3);
    idx->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[4], v[5]);

    d = {Interval(7, 7)};
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[5], v[6]);

    d = {Interval(0, 0)};
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);

    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kh = hash.kmerhash("ACC", 3); // inconsistent
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[8]);

    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kh = hash.kmerhash("AGT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[9]);

    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("ATT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[10]);

    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("GCT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[8], v[11]);

    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kh = hash.kmerhash("GTT", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[9], v[12]);

    d = {Interval(12, 13), Interval(16, 16), Interval(23, 25)};
    p.initialize(d);
    kh = hash.kmerhash("TTA", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[12], v[13]);

    d = {Interval(23, 26)};
    p.initialize(d);
    kh = hash.kmerhash("TAT", 3);//inconsistent but I don't care
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[10], v[14]);
    lp3->kmer_prg.add_edge(v[11], v[14]);
    lp3->kmer_prg.add_edge(v[13], v[14]);

    d = {Interval(24, 27)};
    p.initialize(d);
    kh = hash.kmerhash("ATG", 3);
    idx->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[14], v[15]);

    d = {Interval(27, 27)};
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[15], v[16]);

    d = {Interval(8, 8)};
    p.initialize(d);
    kn = lp4->kmer_prg.add_node(p);
    v.push_back(kn);

    d = {Interval(8, 11)};
    p.initialize(d);
    kh = hash.kmerhash("CTA", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp4->kmer_prg.add_node(p);
    v.push_back(kn);
    lp4->kmer_prg.add_edge(v[17], v[18]);

    d = {Interval(9, 12)};
    p.initialize(d);
    kh = hash.kmerhash("TAG", 3);
    idx->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp4->kmer_prg.add_node(p);
    v.push_back(kn);
    lp4->kmer_prg.add_edge(v[18], v[19]);

    d = {Interval(12, 12)};
    p.initialize(d);
    kn = lp4->kmer_prg.add_node(p);
    v.push_back(kn);
    lp4->kmer_prg.add_edge(v[19], v[20]);

    d = {Interval(0, 0)};
    p.initialize(d);
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);

    d = {Interval(0, 3)};
    p.initialize(d);
    kh = hash.kmerhash("CTA", 3);
    idx->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[21], v[22]);

    d = {Interval(1, 4)};
    p.initialize(d);
    kh = hash.kmerhash("TAC", 3);
    idx->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[22], v[23]);

    d = {Interval(2, 5)};
    p.initialize(d);
    kh = hash.kmerhash("ACT", 3);
    idx->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[23], v[24]);

    d = {Interval(5, 5)};
    p.initialize(d);
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[24], v[25]);

    // initialize pangraph;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pangraph_from_read_file("../../test/test_cases/read2.fa", mhs, pg, idx, prgs, 1, 3, 1, 0.1, 1);

    // create a pangraph object representing the truth we expect (prg 3 4 2 1)
    // note that prgs 1, 3, 4 share no 3mer, but 2 shares a 3mer with each of 2 other prgs
    pangenome::Graph pg_exp;
    MinimizerHits mhs_dummy;
    pg_exp.add_node(1, "1", 0, mhs_dummy.hits);
    pg_exp.add_node(2, "2", 0, mhs_dummy.hits);
    pg_exp.add_node(3, "3", 0, mhs_dummy.hits);
    pg_exp.add_node(0, "0", 0, mhs_dummy.hits);

    EXPECT_EQ(pg_exp, *pg);
    delete pg;

    pg = new pangenome::Graph();
    mhs->clear();
    pangraph_from_read_file("../../test/test_cases/read2.fq", mhs, pg, idx, prgs, 1, 3, 1, 0.1, 1);
    pg_exp.add_node(1, "1", 0, mhs_dummy.hits);
    pg_exp.add_node(2, "2", 0, mhs_dummy.hits);
    pg_exp.add_node(3, "3", 0, mhs_dummy.hits);
    pg_exp.add_node(0, "0", 0, mhs_dummy.hits);

    delete mhs;
    delete lp1;
    delete lp2;
    delete lp3;
    delete lp4;
    idx->clear();
    delete idx;
}

/*
TEST(UtilsTest, lognChoosek2)
{
    EXPECT_EQ(0.0, lognchoosek2(0,0,0));
    EXPECT_EQ(0.0, lognchoosek2(1,0,0));
    EXPECT_EQ(0.0, lognchoosek2(1,1,0));
}

//update_covgs_from_hits
*/
