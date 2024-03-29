#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "utils.h"
#include "localPRG.h"
#include "pangenome/pangraph.h"
#include "interval.h"
#include "prg/path.h"
#include "minihit.h"
#include "minihits.h"
#include "index.h"
#include "inthash.h"
#include "seq.h"
#include <stdint.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include "fatal_error.h"
#include "test_helpers_containers.h"
#include "test_helpers.h"

using namespace std;

const std::string TEST_CASE_DIR = "../../test/test_cases/";

TEST(UtilsTest, split)
{
    vector<string> v = { "abc", "def", "ghi" };
    EXPECT_EQ(v, split("abc, def, ghi", ", "));
    EXPECT_EQ(v, split("abc, def, ghi, ", ", "));
    EXPECT_EQ(v, split(", abc, def, ghi", ", "));
}

TEST(UtilsTest, revComplement)
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

TEST(UtilsTest, readPrgFile)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;

    // simple case first, single prg with empty string sequence
    // doesn't get added to prgs
    read_prg_file(prgs, TEST_CASE_DIR + "prg0.fa");
    uint32_t j = 0;
    EXPECT_EQ(prgs.size(), j);

    // single prg with simple sequence
    read_prg_file(prgs, TEST_CASE_DIR + "prg1.fa");
    LocalPRG l1(1, "prg1", "AGCT");
    j = 1;
    EXPECT_EQ(prgs.size(), j);
    j = 0;
    EXPECT_EQ(prgs[0]->id, j);
    EXPECT_EQ(prgs[0]->name, "prg1");
    EXPECT_EQ(prgs[0]->seq, "AGCT");
    EXPECT_EQ(prgs[0]->prg, l1.prg);

    // single prg with a variant site
    read_prg_file(prgs, TEST_CASE_DIR + "prg2.fa");
    LocalPRG l2(2, "prg2", "A 5 GC 6 G 5 T");
    j = 2;
    EXPECT_EQ(prgs.size(), j);
    j = 0;
    EXPECT_EQ(prgs[1]->id, j);
    EXPECT_EQ(prgs[1]->name, "prg2");
    EXPECT_EQ(prgs[1]->seq, "A 5 GC 6 G 5 T");
    EXPECT_EQ(prgs[1]->prg, l2.prg);

    // single prg with a nested variant site
    read_prg_file(prgs, TEST_CASE_DIR + "prg3.fa");
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
    read_prg_file(prgs, TEST_CASE_DIR + "prg0123.fa");
    j = 3;
    EXPECT_EQ(prgs.size(), j);
}

TEST(UtilsTest, readPrgFile_with_offset)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;

    // simple case first, single prg with empty string sequence
    // doesn't get added to prgs
    read_prg_file(prgs, TEST_CASE_DIR + "prg0.fa", 1);
    EXPECT_EQ(prgs.size(), (uint)0);

    // single prg with simple sequence
    read_prg_file(prgs, TEST_CASE_DIR + "prg1.fa", 1);
    LocalPRG l1(1, "prg1", "AGCT");
    EXPECT_EQ(prgs.size(), (uint)1);
    EXPECT_EQ(prgs[0]->id, (uint)1);
    EXPECT_EQ(prgs[0]->name, "prg1");
    EXPECT_EQ(prgs[0]->seq, "AGCT");
    EXPECT_EQ(prgs[0]->prg, l1.prg);

    // single prg with a variant site
    read_prg_file(prgs, TEST_CASE_DIR + "prg2.fa", 3);
    LocalPRG l2(2, "prg2", "A 5 GC 6 G 5 T");
    EXPECT_EQ(prgs.size(), (uint)2);
    EXPECT_EQ(prgs[1]->id, (uint)3);
    EXPECT_EQ(prgs[1]->name, "prg2");
    EXPECT_EQ(prgs[1]->seq, "A 5 GC 6 G 5 T");
    EXPECT_EQ(prgs[1]->prg, l2.prg);

    // single prg with a nested variant site
    read_prg_file(prgs, TEST_CASE_DIR + "prg3.fa", 40);
    LocalPRG l3 = LocalPRG(3, "prg3", "A 5 G 7 C 8 T 7  6 G 5 T");
    EXPECT_EQ(prgs.size(), (uint)3);
    EXPECT_EQ(prgs[2]->id, (uint)40);
    EXPECT_EQ(prgs[2]->name, "prg3");
    EXPECT_EQ(prgs[2]->seq, "A 5 G 7 C 8 T 7  6 G 5 T");
    EXPECT_EQ(prgs[2]->prg, l3.prg);

    // now a prg input file with all 4 in
    prgs.clear();
    EXPECT_EQ(prgs.size(), (uint)0);
    read_prg_file(prgs, TEST_CASE_DIR + "prg0123.fa", 6);
    EXPECT_EQ(prgs.size(), (uint)3);
    EXPECT_EQ(prgs[0]->id, (uint)6);
    EXPECT_EQ(prgs[1]->id, (uint)7);
    EXPECT_EQ(prgs[2]->id, (uint)8);
}

TEST(UtilsTest, addReadHits)
{
    // initialize minihits container
    auto minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits());
    MinimizerHits expected1;
    MinimizerHits expected2;
    MinimizerHits expected3;
    MinimizerHits expected4;

    // initialize index as we would expect with example prgs 1 and 3 from above
    pandora::KmerHash hash;
    auto index = std::make_shared<Index>();
    deque<Interval> d = { Interval(0, 3) };
    prg::Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("AGC", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    Minimizer min1(0, 0, 3, 1); // kmer, start, end, strand
    MiniRecord mr1(1, p, 0, 1);
    MinimizerHitPtr m1(make_shared<MinimizerHit>(0, min1, mr1));

    Minimizer min2(0, 1, 4, 0); // kmer, start, end, strand
    MiniRecord mr2(1, p, 0, 1);
    MinimizerHitPtr m2(make_shared<MinimizerHit>(0, min2, mr2));
    expected1.insert(m1);
    expected2.insert(m2);
    d = { Interval(1, 4) };
    p.initialize(d);
    kh = hash.kmerhash("GCT", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));

    Minimizer min3(0, 1, 4, 1); // kmer, start, end, strand
    MiniRecord mr3(1, p, 0, 1);
    MinimizerHitPtr m3(make_shared<MinimizerHit>(0, min3, mr3));

    Minimizer min4(0, 0, 3, 0); // kmer, start, end, strand
    MiniRecord mr4(1, p, 0, 1);
    MinimizerHitPtr m4(make_shared<MinimizerHit>(0, min4, mr4));

    expected2.insert(m3);
    expected1.insert(m4);
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    p.initialize(d);
    kh = hash.kmerhash("AGC", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    Minimizer min5(0, 0, 3, 1); // kmer, start, end, strand
    MiniRecord mr5(3, p, 0, 1);
    MinimizerHitPtr m5(make_shared<MinimizerHit>(0, min5, mr5));

    Minimizer min6(0, 1, 4, 0); // kmer, start, end, strand
    MiniRecord mr6(3, p, 0, 1);
    MinimizerHitPtr m6(make_shared<MinimizerHit>(0, min6, mr6));
    expected1.insert(m5);
    expected2.insert(m6);
    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13) };
    p.initialize(d);
    kh = hash.kmerhash("AGT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));

    Minimizer min9(0, 0, 3, 1); // kmer, start, end, strand
    MiniRecord mr9(3, p, 0, 1);
    MinimizerHitPtr m9(make_shared<MinimizerHit>(0, min9, mr9));
    expected3.insert(m9);
    d = { Interval(0, 1), Interval(19, 20), Interval(23, 24) };
    p.initialize(d);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));

    Minimizer min10(0, 0, 3, 1); // kmer, start, end, strand
    MiniRecord mr10(3, p, 0, 1);
    MinimizerHitPtr m10(make_shared<MinimizerHit>(0, min10, mr10));
    expected3.insert(m10);
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("GCT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));

    Minimizer min7(0, 1, 4, 1); // kmer, start, end, strand
    MiniRecord mr7(3, p, 0, 1);
    MinimizerHitPtr m7(make_shared<MinimizerHit>(0, min7, mr7));

    Minimizer min8(0, 0, 3, 0); // kmer, start, end, strand
    MiniRecord mr8(3, p, 0, 1);
    MinimizerHitPtr m8(make_shared<MinimizerHit>(0, min8, mr8));
    expected2.insert(m7);
    expected1.insert(m8);
    d = { Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("GTT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));

    Minimizer min11(0, 1, 4, 1); // kmer, start, end, strand
    MiniRecord mr11(3, p, 0, 1);
    MinimizerHitPtr m11(make_shared<MinimizerHit>(0, min11, mr11));
    expected4.insert(m11);

    Seq s(0, "read1", "AGC", 1, 3);
    add_read_hits(s, minimizer_hits, *index);
    EXPECT_EQ(expected1.size(), minimizer_hits->size());
    auto it2 = expected1.begin();
    for (auto it = minimizer_hits->begin();
         it != minimizer_hits->end(); ++it) {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    // if take w=2 as sketch of read AGTT should miss AGT, which occurs twice in PRG and
    // contain GTT

    minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits());
    uint32_t j = 0;
    EXPECT_EQ(j, minimizer_hits->size());
    s = Seq(0, "read2", "AGTT", 2, 3);
    add_read_hits(s, minimizer_hits, *index);
    EXPECT_EQ(expected4.size(), minimizer_hits->size());
    it2 = expected4.begin();
    for (auto it = minimizer_hits->begin();
         it != minimizer_hits->end(); ++it) {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    // but for w=1, only add one more hit, for GTT
    expected3.insert(m11);

    minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits());
    EXPECT_EQ(j, minimizer_hits->size());
    s = Seq(0, "read2", "AGTT", 1, 3);
    add_read_hits(s, minimizer_hits, *index);
    EXPECT_EQ(expected3.size(), minimizer_hits->size());
    it2 = expected3.begin();
    for (auto it = minimizer_hits->begin();
         it != minimizer_hits->end(); ++it) {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    // now back to w = 1, add expected2 to expected1 as will get hits against both AGC
    // and GCT
    minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits());
    j = 0;
    EXPECT_EQ(j, minimizer_hits->size());
    s = Seq(0, "read3", "AGCT", 1, 3);
    add_read_hits(s, minimizer_hits, *index);

    expected1.insert(expected2.begin(), expected2.end());
    EXPECT_EQ(expected1.size(), minimizer_hits->size());
    it2 = expected1.begin();
    for (auto it = minimizer_hits->begin(); it != minimizer_hits->end();
         ++it) {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    // same for w = 2, add expected2 to expected1 as will get hits against both because
    // AGC and GCT are joint minimums
    minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits());
    j = 0;
    EXPECT_EQ(j, minimizer_hits->size());
    s = Seq(0, "read3", "AGCT", 2, 3);
    add_read_hits(s, minimizer_hits, *index);
    EXPECT_EQ(expected1.size(), minimizer_hits->size());
    it2 = expected1.begin();
    for (auto it = minimizer_hits->begin(); it != minimizer_hits->end();
         ++it) {
        EXPECT_EQ(**it2, **it);
        it2++;
    }

    expected1.clear();
    expected2.clear();
    expected3.clear();
    expected4.clear();
    index->clear();
}

TEST(UtilsTest, simpleInferLocalPRGOrderForRead)
{
    // initialize minihits container
    auto minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits());
    pandora::KmerHash hash;

    // initialize a prgs object
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    auto lp1 = std::make_shared<LocalPRG>(LocalPRG(1, "1", ""));
    auto lp3 = std::make_shared<LocalPRG>(LocalPRG(0, "0", ""));
    prgs.push_back(lp3);
    prgs.push_back(lp1);

    // initialize index as we would expect with example prgs (variant of) 1 and 3 from
    // above
    auto index = std::make_shared<Index>();

    vector<KmerNodePtr> v;
    KmerNodePtr kn;

    deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);

    d = { Interval(0, 3) };
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("TAC", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[0], v[1]);

    d = { Interval(1, 4) };
    p.initialize(d);
    kh = hash.kmerhash("ACG", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[1], v[2]);

    d = { Interval(4, 4) };
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[2], v[3]);

    d = { Interval(0, 0) };
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);

    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    p.initialize(d);
    kh = hash.kmerhash("AGC", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[4], v[5]);

    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13) };
    p.initialize(d);
    kh = hash.kmerhash("AGT", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[4], v[6]);

    d = { Interval(0, 1), Interval(19, 20), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("ATT", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[4], v[7]);

    d = { Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("GCT", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[5], v[8]);

    d = { Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("GTT", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[6], v[9]);

    d = { Interval(12, 13), Interval(16, 16), Interval(23, 25) };
    p.initialize(d);
    kh = hash.kmerhash("TTA", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[9], v[10]);

    d = { Interval(23, 26) };
    p.initialize(d);
    kh = hash.kmerhash("TAA", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[11]);
    lp3->kmer_prg.add_edge(v[8], v[11]);
    lp3->kmer_prg.add_edge(v[10], v[11]);

    d = { Interval(24, 27) };
    p.initialize(d);
    kh = hash.kmerhash("AAG", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[11], v[12]);

    d = { Interval(27, 27) };
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[12], v[13]);

    // add read hits to mhs
    Seq s(0, "read1", "AGTTAAGTACG", 1, 3);
    add_read_hits(s, minimizer_hits, *index);
    // add_read_hits(0, "read1", "AGTTAAGTACG", mhs, index, 1, 3);

    // initialize pangraph;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    ClusterDefFile cluster_def_file("", true);
    ClusterFilterFile cluster_filter_file("", true);
    MinimizerHitClusters minimizer_hit_clusters = get_minimizer_hit_clusters(
        "sample", s, prgs, minimizer_hits, pangraph, 1, 100, 0.1,
        cluster_def_file, cluster_filter_file, 1);
    add_clusters_to_pangraph(minimizer_hit_clusters, pangraph, prgs);

    // create a pangraph object representing the truth we expect (prg 3 then 1)
    pangenome::Graph pg_exp;
    pg_exp.add_node(lp1);
    pg_exp.add_node(lp3);
    // pg_exp.add_edge(0,1,3,0);

    EXPECT_EQ(pg_exp, *pangraph);
    index->clear();
}

TEST(UtilsTest, biggerInferLocalPRGOrderForRead)
{
    // initialize minihits container
    auto minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits());
    pandora::KmerHash hash;

    // initialize a prgs object
    std::vector<std::shared_ptr<LocalPRG>> prgs;

    auto lp1 = std::make_shared<LocalPRG>(LocalPRG(1, "1", ""));
    auto lp3 = std::make_shared<LocalPRG>(LocalPRG(3, "3", ""));
    auto lp0 = std::make_shared<LocalPRG>(LocalPRG(0, "", ""));
    auto lp2 = std::make_shared<LocalPRG>(LocalPRG(2, "2", ""));
    prgs.push_back(lp0);
    prgs.push_back(lp1);
    prgs.push_back(lp2);
    prgs.push_back(lp3);

    // initialize index as we would expect with example prgs
    auto index = std::make_shared<Index>();

    vector<KmerNodePtr> v;
    KmerNodePtr kn;

    deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);

    d = { Interval(0, 3) };
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("TAC", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[0], v[1]);

    d = { Interval(1, 4) };
    p.initialize(d);
    kh = hash.kmerhash("ACG", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[1], v[2]);

    d = { Interval(2, 5) };
    p.initialize(d);
    kh = hash.kmerhash("CGG", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[2], v[3]);

    d = { Interval(3, 6) };
    p.initialize(d);
    kh = hash.kmerhash("GGT", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[3], v[4]);

    d = { Interval(4, 7) };
    p.initialize(d);
    kh = hash.kmerhash("GTA", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[4], v[5]);

    d = { Interval(7, 7) };
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[5], v[6]);

    d = { Interval(0, 0) };
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);

    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    p.initialize(d);
    kh = hash.kmerhash("ACC", 3); // inconsistent
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[8]);

    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13) };
    p.initialize(d);
    kh = hash.kmerhash("AGT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[9]);

    d = { Interval(0, 1), Interval(19, 20), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("ATT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[10]);

    d = { Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("GCT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[8], v[11]);

    d = { Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("GTT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[9], v[12]);

    d = { Interval(12, 13), Interval(16, 16), Interval(23, 25) };
    p.initialize(d);
    kh = hash.kmerhash("TTA", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[12], v[13]);

    d = { Interval(23, 26) };
    p.initialize(d);
    kh = hash.kmerhash("TAT", 3); // inconsistent but I don't care
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[10], v[14]);
    lp3->kmer_prg.add_edge(v[11], v[14]);
    lp3->kmer_prg.add_edge(v[13], v[14]);

    d = { Interval(24, 27) };
    p.initialize(d);
    kh = hash.kmerhash("ATG", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[14], v[15]);

    d = { Interval(27, 27) };
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[15], v[16]);

    d = { Interval(8, 8) };
    p.initialize(d);
    kn = lp0->kmer_prg.add_node(p);
    v.push_back(kn);

    d = { Interval(8, 11) };
    p.initialize(d);
    kh = hash.kmerhash("CTA", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp0->kmer_prg.add_node(p);
    v.push_back(kn);
    lp0->kmer_prg.add_edge(v[17], v[18]);

    d = { Interval(9, 12) };
    p.initialize(d);
    kh = hash.kmerhash("TAG", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp0->kmer_prg.add_node(p);
    v.push_back(kn);
    lp0->kmer_prg.add_edge(v[18], v[19]);

    d = { Interval(12, 12) };
    p.initialize(d);
    kn = lp0->kmer_prg.add_node(p);
    v.push_back(kn);
    lp0->kmer_prg.add_edge(v[19], v[20]);

    d = { Interval(0, 0) };
    p.initialize(d);
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);

    d = { Interval(0, 3) };
    p.initialize(d);
    kh = hash.kmerhash("CTA", 3);
    index->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[21], v[22]);

    d = { Interval(1, 4) };
    p.initialize(d);
    kh = hash.kmerhash("TAC", 3);
    index->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[22], v[23]);

    d = { Interval(2, 5) };
    p.initialize(d);
    kh = hash.kmerhash("ACT", 3);
    index->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[23], v[24]);

    d = { Interval(5, 5) };
    p.initialize(d);
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[24], v[25]);

    // add read hits to mhs
    Seq s(0, "read2", "AGTTATGCTAGCTACTTACGGTA", 1, 3);
    add_read_hits(s, minimizer_hits, *index);
    // add_read_hits(0, "read2", "AGTTATGCTAGCTACTTACGGTA", mhs, index, 1, 3);

    // initialize pangraph;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    ClusterDefFile cluster_def_file("", true);
    ClusterFilterFile cluster_filter_file("", true);
    MinimizerHitClusters minimizer_hit_clusters = get_minimizer_hit_clusters(
        "sample", s, prgs, minimizer_hits, pangraph, 1, 100, 0.1,
        cluster_def_file, cluster_filter_file, 1);
    add_clusters_to_pangraph(minimizer_hit_clusters, pangraph, prgs);

    // create a pangraph object representing the truth we expect (prg 3 4 2 1)
    // note that prgs 1, 3, 4 share no 3mer, but 2 shares a 3mer with each of 2 other
    // prgs
    pangenome::Graph pg_exp;
    pg_exp.add_node(lp1);
    pg_exp.add_node(lp2);
    pg_exp.add_node(lp3);
    pg_exp.add_node(lp0);
    // pg_exp.add_edge(3,0,3,0);
    // pg_exp.add_edge(0,2,3,0);
    // pg_exp.add_edge(2,1,3,0);

    EXPECT_EQ(pg_exp, *pangraph);
    index->clear();
}

/*TEST(UtilsTest, read_all_readfile) {

    MinimizerHits *mhs;
    mhs = new MinimizerHits();
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    Index *index;
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    pangraph_from_read_file(TEST_CASE_DIR + "reads.fq.gz", mhs, pg, index, prgs,
1, 3, 1, 0.1); delete mhs; delete pg;
}*/

void setup_index(
    std::vector<std::shared_ptr<LocalPRG>>& prgs, std::shared_ptr<Index>& index)
{
    auto lp1 = std::make_shared<LocalPRG>(LocalPRG(1, "1", ""));
    auto lp3 = std::make_shared<LocalPRG>(LocalPRG(3, "3", ""));
    auto lp0 = std::make_shared<LocalPRG>(LocalPRG(0, "0", ""));
    auto lp2 = std::make_shared<LocalPRG>(LocalPRG(2, "2", ""));
    prgs.push_back(lp0);
    prgs.push_back(lp1);
    prgs.push_back(lp2);
    prgs.push_back(lp3);

    pandora::KmerHash hash;
    vector<KmerNodePtr> v;
    KmerNodePtr kn;

    deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);

    d = { Interval(0, 3) };
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("TAC", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[0], v[1]);

    d = { Interval(1, 4) };
    p.initialize(d);
    kh = hash.kmerhash("ACG", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[1], v[2]);

    d = { Interval(2, 5) };
    p.initialize(d);
    kh = hash.kmerhash("CGG", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[2], v[3]);

    d = { Interval(3, 6) };
    p.initialize(d);
    kh = hash.kmerhash("GGT", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[3], v[4]);

    d = { Interval(4, 7) };
    p.initialize(d);
    kh = hash.kmerhash("GTA", 3);
    index->add_record(min(kh.first, kh.second), 1, p, 0, (kh.first < kh.second));
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[4], v[5]);

    d = { Interval(7, 7) };
    p.initialize(d);
    kn = lp1->kmer_prg.add_node(p);
    v.push_back(kn);
    lp1->kmer_prg.add_edge(v[5], v[6]);

    d = { Interval(0, 0) };
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);

    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    p.initialize(d);
    kh = hash.kmerhash("ACC", 3); // inconsistent
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[8]);

    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13) };
    p.initialize(d);
    kh = hash.kmerhash("AGT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[9]);

    d = { Interval(0, 1), Interval(19, 20), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("ATT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[7], v[10]);

    d = { Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("GCT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[8], v[11]);

    d = { Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kh = hash.kmerhash("GTT", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[9], v[12]);

    d = { Interval(12, 13), Interval(16, 16), Interval(23, 25) };
    p.initialize(d);
    kh = hash.kmerhash("TTA", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[12], v[13]);

    d = { Interval(23, 26) };
    p.initialize(d);
    kh = hash.kmerhash("TAT", 3); // inconsistent but I don't care
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[10], v[14]);
    lp3->kmer_prg.add_edge(v[11], v[14]);
    lp3->kmer_prg.add_edge(v[13], v[14]);

    d = { Interval(24, 27) };
    p.initialize(d);
    kh = hash.kmerhash("ATG", 3);
    index->add_record(min(kh.first, kh.second), 3, p, 0, (kh.first < kh.second));
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[14], v[15]);

    d = { Interval(27, 27) };
    p.initialize(d);
    kn = lp3->kmer_prg.add_node(p);
    v.push_back(kn);
    lp3->kmer_prg.add_edge(v[15], v[16]);

    d = { Interval(8, 8) };
    p.initialize(d);
    kn = lp0->kmer_prg.add_node(p);
    v.push_back(kn);

    d = { Interval(8, 11) };
    p.initialize(d);
    kh = hash.kmerhash("CTA", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp0->kmer_prg.add_node(p);
    v.push_back(kn);
    lp0->kmer_prg.add_edge(v[17], v[18]);

    d = { Interval(9, 12) };
    p.initialize(d);
    kh = hash.kmerhash("TAG", 3);
    index->add_record(min(kh.first, kh.second), 0, p, 0, (kh.first < kh.second));
    kn = lp0->kmer_prg.add_node(p);
    v.push_back(kn);
    lp0->kmer_prg.add_edge(v[18], v[19]);

    d = { Interval(12, 12) };
    p.initialize(d);
    kn = lp0->kmer_prg.add_node(p);
    v.push_back(kn);
    lp0->kmer_prg.add_edge(v[19], v[20]);

    d = { Interval(0, 0) };
    p.initialize(d);
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);

    d = { Interval(0, 3) };
    p.initialize(d);
    kh = hash.kmerhash("CTA", 3);
    index->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[21], v[22]);

    d = { Interval(1, 4) };
    p.initialize(d);
    kh = hash.kmerhash("TAC", 3);
    index->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[22], v[23]);

    d = { Interval(2, 5) };
    p.initialize(d);
    kh = hash.kmerhash("ACT", 3);
    index->add_record(min(kh.first, kh.second), 2, p, 0, (kh.first < kh.second));
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[23], v[24]);

    d = { Interval(5, 5) };
    p.initialize(d);
    kn = lp2->kmer_prg.add_node(p);
    v.push_back(kn);
    lp2->kmer_prg.add_edge(v[24], v[25]);
}

TEST(UtilsTest, pangraphFromReadFile_Fa)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;

    auto index = std::make_shared<Index>();
    setup_index(prgs, index);

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    SampleData sample_data{"reads_2", TEST_CASE_DIR + "read2.fa"};

    fs::path outdir("reads_2_outdir");
    fs::create_directories(outdir);
    pangraph_from_read_file(sample_data, pangraph, index, prgs, 1, 3, 1,
        0.1, outdir, 1);
    fs::remove_all(outdir);

    // create a pangraph object representing the truth we expect (prg 3 4 2 1)
    // note that prgs 1, 3, 4 share no 3mer, but 2 shares a 3mer with each of 2 other
    // prgs
    pangenome::Graph pg_exp;
    pg_exp.add_node(prgs[1]);
    pg_exp.add_node(prgs[2]);
    pg_exp.add_node(prgs[3]);
    pg_exp.add_node(prgs[0]);

    EXPECT_EQ(pg_exp, *pangraph);

    index->clear();
}

TEST(UtilsTest, pangraphFromReadFile_Fq)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;

    auto index = std::make_shared<Index>();
    setup_index(prgs, index);

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    fs::path outdir("reads_2_outdir");
    fs::create_directories(outdir);
    SampleData sample_data{"reads_2", TEST_CASE_DIR + "read2.fa"};
    pangraph_from_read_file(sample_data, pangraph, index, prgs, 1, 3, 1,
        0.1, outdir, 1);
    fs::remove_all(outdir);

    // create a pangraph object representing the truth we expect (prg 3 4 2 1)
    // note that prgs 1, 3, 4 share no 3mer, but 2 shares a 3mer with each of 2 other
    // prgs
    pangenome::Graph pg_exp;
    pg_exp.add_node(prgs[1]);
    pg_exp.add_node(prgs[2]);
    pg_exp.add_node(prgs[3]);
    pg_exp.add_node(prgs[0]);

    EXPECT_EQ(pg_exp, *pangraph);

    index->clear();
}

TEST(StrToGsTest, HandlesEmptyStr)
{
    const char* str { "" };

    const auto actual { strtogs(str) };
    const auto expected { 0 };

    EXPECT_EQ(actual, expected);
}

TEST(StrToGsTest, DecimalConvertsToInt)
{
    const char* str { "3.3" };

    const auto actual { strtogs(str) };
    const auto expected { 3 };

    EXPECT_EQ(actual, expected);
}

TEST(StrToGsTest, DecimalRoundsUpToInt)
{
    const char* str { "3.7" };

    const auto actual { strtogs(str) };
    const auto expected { 4 };

    EXPECT_EQ(actual, expected);
}

TEST(StrToGsTest, IntReturnsSame)
{
    const char* str { "310" };

    const auto actual { strtogs(str) };
    const auto expected { 310 };

    EXPECT_EQ(actual, expected);
}

TEST(StrToGsTest, NegativeFails)
{
    const char* str { "-4" };

    EXPECT_THROW(strtogs(str), std::logic_error);
}

TEST(StrToGsTest, Kilo)
{
    const auto* str { "3.1k" };

    auto actual { strtogs(str) };
    const auto expected { 3100 };

    EXPECT_EQ(actual, expected);

    str = "3.1K";

    actual = strtogs(str);

    EXPECT_EQ(actual, expected);
}

TEST(StrToGsTest, Mega)
{
    const auto* str { "3m" };

    auto actual { strtogs(str) };
    const auto expected { 3e6 };

    EXPECT_EQ(actual, expected);

    str = "3M";

    actual = strtogs(str);

    EXPECT_EQ(actual, expected);
}

TEST(StrToGsTest, Giga)
{
    const auto* str { "3g" };

    auto actual { strtogs(str) };
    const auto expected { 3e9 };

    EXPECT_EQ(actual, expected);

    str = "3G";

    actual = strtogs(str);

    EXPECT_EQ(actual, expected);
}

TEST(StrToGsTest, TooLargeFails)
{
    const auto* str { "20g" };

    EXPECT_THROW(strtogs(str), std::runtime_error);
}

TEST(TransformCliGsizeTest, GenomeSizeStringProducesIntToStr)
{
    const std::string str { "2k" };

    const auto actual { transform_cli_gsize(str) };
    const auto* const expected { "2000" };

    EXPECT_EQ(actual, expected);
}

TEST(MakeAbsoluteTest, CurrentDirAsPeriodReturnsAbsolutePath)
{
    const std::string str { "." };

    const auto actual { make_absolute(str) };
    const auto expected { fs::current_path() / "." };

    EXPECT_EQ(actual, expected.string());
}

TEST(MakeAbsoluteTest, RelativeDirReturnsAbsolutePath)
{
    const std::string str { "../abs" };

    const auto actual { make_absolute(str) };
    const auto expected { fs::current_path() / "../abs" };

    EXPECT_EQ(actual, expected.string());
}

TEST(MakeAbsoluteTest, FileReturnsAbsolutePathToFile)
{
    const std::string str { "foo.bar" };

    const auto actual { make_absolute(str) };
    const auto expected { fs::current_path() / str };

    EXPECT_EQ(actual, expected.string());
}

TEST(load_read_index, read_index_does_not_exist___expects_FatalRuntimeError)
{
    ASSERT_EXCEPTION(load_read_index(fs::path("nonexistent_read_index.tsv")),
        FatalRuntimeError, "Unable to open read index file");
}

TEST(load_read_index, read_index_has_three_samples)
{
    std::vector<SampleData> actual
        = load_read_index(fs::path("../../test/test_cases/sample_read_index.tsv"));
    std::vector<SampleData> expected { {
        std::make_pair("sample_1", "reads_1.fastq"),
        std::make_pair("sample_2", "reads_2.fastq"),
        std::make_pair("sample_3", "reads_3.fastq"),
    } };

    EXPECT_EQ(actual, expected);
}

TEST(load_read_index, read_index_has_three_samples_and_no_empty_line_at_end)
{
    std::vector<std::pair<SampleIdText, SampleFpath>> actual = load_read_index(
        fs::path("../../test/test_cases/sample_read_index_no_empty_line_at_end.tsv"));
    std::vector<std::pair<SampleIdText, SampleFpath>> expected { {
        std::make_pair("sample_1", "reads_1.fastq"),
        std::make_pair("sample_2", "reads_2.fastq"),
        std::make_pair("sample_3", "reads_3.fastq"),
    } };

    EXPECT_EQ(actual, expected);
}

TEST(load_read_index, read_index_has_three_samples_and_two_are_repeated)
{
    std::vector<SampleData> actual = load_read_index(
        fs::path("../../test/test_cases/sample_read_index_with_repeated_samples.tsv"));
    std::vector<SampleData> expected { {
        std::make_pair("sample_1", "first_reads_1.fastq"),
        std::make_pair("sample_2", "second_reads_2.fastq"),
        std::make_pair("sample_3", "fourth_reads_3.fastq"),
    } };

    EXPECT_EQ(actual, expected);
}

TEST(load_read_index, read_index_has_missing_column)
{
    ASSERT_EXCEPTION(
        load_read_index(fs::path("../../test/test_cases/malformatted_read_index.tsv")),
        FatalRuntimeError, "Malformatted read index file entry for sample_3");
}

TEST(splitAmbiguous, noAmbiguous)
{
    const std::string s("ACGT");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs { s };
    const std::vector<std::size_t> expected_offsets { 0 };
    const auto expected = make_pair(expected_substrs, expected_offsets);

    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, emptySequence)
{
    const std::string s("");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs;
    const std::vector<std::size_t> expected_offsets;
    const auto expected = make_pair(expected_substrs, expected_offsets);


    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, allAmbiguous)
{
    const std::string s("NXDW");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs;
    const std::vector<std::size_t> expected_offsets;
    const auto expected = make_pair(expected_substrs, expected_offsets);


    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, firstLetterIsAmbiguous)
{
    const std::string s("NACGT");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs { "ACGT" };
    const std::vector<std::size_t> expected_offsets { 1 };
    const auto expected = make_pair(expected_substrs, expected_offsets);

    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, firstTwoLettersAreAmbiguous)
{
    const std::string s("NWACGT");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs { "ACGT" };
    const std::vector<std::size_t> expected_offsets { 2 };
    const auto expected = make_pair(expected_substrs, expected_offsets);

    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, lastLetterIsAmbiguous)
{
    const std::string s("ACGTN");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs { "ACGT" };
    const std::vector<std::size_t> expected_offsets { 0 };
    const auto expected = make_pair(expected_substrs, expected_offsets);

    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, lastTwoLettersAreAmbiguous)
{
    const std::string s("ACGTNW");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs { "ACGT" };
    const std::vector<std::size_t> expected_offsets { 0 };
    const auto expected = make_pair(expected_substrs, expected_offsets);

    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, ambiguousBaseInMiddle)
{
    const std::string s("ACNGT");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs { "AC", "GT" };
    const std::vector<std::size_t> expected_offsets { 0, 3 };
    const auto expected = make_pair(expected_substrs, expected_offsets);

    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, ambiguousBaseOffCentre)
{
    const std::string s("AWCGT");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs { "A", "CGT" };
    const std::vector<std::size_t> expected_offsets { 0, 2 };
    const auto expected = make_pair(expected_substrs, expected_offsets);

    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, twoAmbiguousInMiddle)
{
    const std::string s("AWXCGT");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs { "A", "CGT" };
    const std::vector<std::size_t> expected_offsets { 0, 3 };
    const auto expected = make_pair(expected_substrs, expected_offsets);

    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, twoAmbiguousSpacedOut)
{
    const std::string s("AWCNGT");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs { "A", "C", "GT" };
    const std::vector<std::size_t> expected_offsets { 0, 2, 4 };
    const auto expected = make_pair(expected_substrs, expected_offsets);

    EXPECT_EQ(actual, expected);
}

TEST(splitAmbiguous, twoAmbiguousSpacedOutRuns)
{
    const std::string s("AWXCNXGT");

    const auto actual(split_ambiguous(s));

    const std::vector<std::string> expected_substrs { "A", "C", "GT" };
    const std::vector<std::size_t> expected_offsets { 0, 3, 6 };
    const auto expected = make_pair(expected_substrs, expected_offsets);


    EXPECT_EQ(actual, expected);
}
