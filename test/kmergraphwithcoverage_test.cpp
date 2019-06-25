#include "gtest/gtest.h"
#include "test_macro.cpp"

#include "interval.h"
#include "prg/path.h"
#include "kmergraphwithcoverage.h"
#include "kmernode.h"
#include "localPRG.h"
#include <stdint.h>
#include <iostream>
#include <cmath>

using namespace prg;

/*TEST(KmerGraphWithCoverageTest, equals) {
    KmerGraph kg1, kg2;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kg1.add_node(p1);
    auto m1 = kg2.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg1.add_node(p2);
    auto m2 = kg2.add_node(p2);
    kg1.add_edge(n1, n2);
    kg2.add_edge(m1, m2);

    d = {Interval(2, 5)};
    p3.initialize(d);
    auto m3 = kg2.add_node(p3);

    // same as themselves, different if different numbers of nodes
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1 == kg2), false);
    EXPECT_EQ((kg2 == kg1), false);

    auto n3 = kg1.add_node(p3);
    kg2.add_edge(m1, m3);

    // same as themselves, different if different numbers of edges
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1 == kg2), false);
    EXPECT_EQ((kg2 == kg1), false);

    kg1.add_edge(n2, n3);

    // same as themselves, different if edges in different places
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1 == kg2), false);
    EXPECT_EQ((kg2 == kg1), false);
}

TEST(KmerGraphWithCoverageTest, copy) {
    KmerGraphWithCoverage kg1;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kg1.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg1.add_node(p2);
    kg1.add_edge(n1, n2);

    KmerGraphWithCoverage kg2(kg1);

    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
}

TEST(KmerGraphWithCoverageTest, assign) {
    KmerGraphWithCoverage kg1;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    auto n = kg1.add_node(p);
    d = {Interval(0, 3)};
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kg1.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg1.add_node(p2);
    kg1.add_edge(n1, n2);
    d = {Interval(11, 14)};
    p3.initialize(d);
    auto n3 = kg1.add_node(p3);
    kg1.add_edge(n1, n3);
    d = {Interval(15, 18)};
    p1.initialize(d);
    n1 = kg1.add_node(p1);
    kg1.add_edge(n2, n1);
    d = {Interval(20, 20)};
    p.initialize(d);
    n = kg1.add_node(p);
    kg1.add_edge(n1, n);
    kg1.add_edge(n3, n);

    KmerGraphWithCoverage kg2 = kg1;

    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
}*/

TEST(KmerGraphWithCoverageTest, set_exp_depth_covg){
    KmerGraph kg;
    KmerGraphWithCoverage kgc(&kg);
    kgc.set_exp_depth_covg(10);
    EXPECT_EQ(kgc.exp_depth_covg, (uint)10);
}

TEST(KmerGraphWithCoverageTest, set_p) {
    KmerGraph kg;
    KmerGraphWithCoverage kgc(&kg);
    EXPECT_DEATH(kgc.set_p(0.4), "");
    kgc.kmer_prg->k = 3;
    EXPECT_DEATH(kgc.set_p(0), "");
    EXPECT_DEATH(kgc.set_p(1), "");
    kgc.set_p(0.5);
    EXPECT_EQ(1 / exp(1.5) - 0.00001 <= kgc.p and 1 / exp(1.5) + 0.00001 >= kgc.p, true);
}

TEST(KmerGraphWithCoverageTest, set_nb) {
    KmerGraph kg;
    KmerGraphWithCoverage kgc(&kg);
    kgc.set_nb(0,0);
    EXPECT_FLOAT_EQ(kgc.nb_p, 0.015); //unchanged
}

TEST(KmerGraphWithCoverageTest, prob_failNoNodes) {
    uint32_t sample_id = 0;
    KmerGraph kg;
    KmerGraphWithCoverage kgc(&kg);
    EXPECT_DEATH(kgc.prob(0, sample_id), "");
}

TEST(KmerGraphWithCoverageTest, prob_failNoP) {
    uint32_t sample_id = 0;
    KmerGraph kg;

    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    KmerGraphWithCoverage kgc(&kg);

    EXPECT_DEATH(kgc.prob(0, sample_id), "");
}

TEST(KmerGraphWithCoverageTest, prob_failNoNumReads) {
    uint32_t sample_id = 0;
    KmerGraph kg;

    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);

    KmerGraphWithCoverage kgc(&kg);
    kgc.kmer_prg->k = 3;
    kgc.set_p(0.5);

    EXPECT_DEATH(kgc.prob(0, sample_id), "");
}

TEST(KmerGraphWithCoverageTest, prob_simple) {
    uint32_t sample_id = 0;
    KmerGraph kg;

    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);

    KmerGraphWithCoverage kgc(&kg);
    kgc.kmer_prg->k = 3;
    kgc.set_p(0.5);
    kgc.num_reads = 1;

    EXPECT_EQ(kgc.kmer_prg->nodes.size(), (uint) 1);
    EXPECT_EQ(0, kgc.prob(0, sample_id));
}

TEST(KmerGraphWithCoverageTest, prob_realNodeCovgs) {
    uint32_t sample_id = 0;
    KmerGraph kg;

    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);

    KmerGraphWithCoverage kgc(&kg);
    kgc.kmer_prg->k = 3;
    kgc.set_p(0.5);
    kgc.num_reads = 1;

    EXPECT_EQ(kgc.kmer_prg->nodes.size(), (uint) 3);

    EXPECT_EQ(kgc.prob(1, sample_id), kgc.prob(1, sample_id));
    EXPECT_EQ(kgc.prob(2, sample_id), kgc.prob(2, sample_id));

}


TEST(KmerGraphWithCoverageTest, findMaxPath_InvalidProbModel) {
    KmerGraph kg;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 7;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[0], kg.nodes[3]);
    kg.add_edge(kg.nodes[3], kg.nodes[4]);
    kg.add_edge(kg.nodes[0], kg.nodes[5]);
    kg.add_edge(kg.nodes[2], kg.nodes[6]);
    kg.add_edge(kg.nodes[4], kg.nodes[6]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);

    KmerGraphWithCoverage kgc(&kg);
    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kgc.set_covg(1,4, 0, sample_id);
    kgc.set_covg(2,3, 0, sample_id);

    kgc.num_reads = 5;
    kgc.kmer_prg->k = 3;

    vector<KmerNodePtr> mp;
    kgc.set_p(0.01);
    EXPECT_DEATH(kgc.find_max_path(mp, "exp", max_num_kmers_to_average, sample_id),"");
}

TEST(KmerGraphWithCoverageTest, findMaxPathSimple) {
    KmerGraph kg;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 7;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[0], kg.nodes[3]);
    kg.add_edge(kg.nodes[3], kg.nodes[4]);
    kg.add_edge(kg.nodes[0], kg.nodes[5]);
    kg.add_edge(kg.nodes[2], kg.nodes[6]);
    kg.add_edge(kg.nodes[4], kg.nodes[6]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);

    KmerGraphWithCoverage kgc(&kg);
    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kgc.set_covg(1,4, 0, sample_id);
    kgc.set_covg(2,3, 0, sample_id);

    kgc.num_reads = 5;
    kgc.kmer_prg->k = 3;

    vector<KmerNodePtr> mp;
    kgc.set_p(0.01);
    kgc.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kg.nodes[1], kg.nodes[2]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kgc.set_covg(1, 0, 0, sample_id);
    kgc.set_covg(2, 0, 0, sample_id);
    kgc.set_covg(5, 5, 1, sample_id);
    kgc.set_p(0.01);
    kgc.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    exp_order = {kg.nodes[5]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);
}

TEST(KmerGraphWithCoverageTest, findMaxPathSimple_WithMaxKmersInAvg) {
    KmerGraph kg;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 7;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[0], kg.nodes[3]);
    kg.add_edge(kg.nodes[3], kg.nodes[4]);
    kg.add_edge(kg.nodes[0], kg.nodes[5]);
    kg.add_edge(kg.nodes[2], kg.nodes[6]);
    kg.add_edge(kg.nodes[4], kg.nodes[6]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);

    KmerGraphWithCoverage kgc(&kg);
    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 1;

    kgc.set_covg(1,4, 0, sample_id);
    kgc.set_covg(2,3, 0, sample_id);

    kgc.num_reads = 5;
    kgc.kmer_prg->k = 3;

    vector<KmerNodePtr> mp;
    kgc.set_p(0.01);
    kgc.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kg.nodes[1], kg.nodes[2]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kgc.set_covg(1, 0, 0, sample_id);
    kgc.set_covg(2, 0, 0, sample_id);
    kgc.set_covg(5, 5, 1, sample_id);
    kgc.set_p(0.01);
    kgc.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    exp_order = {kg.nodes[5]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);
}
/*
TEST(KmerGraphWithCoverageTest, findMaxPath2Level_bin) {
    KmerGraph kg;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(8, 9), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(12, 13), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(16, 18), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 10;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[2], kg.nodes[3]);
    kg.add_edge(kg.nodes[0], kg.nodes[4]);
    kg.add_edge(kg.nodes[4], kg.nodes[5]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);
    kg.add_edge(kg.nodes[3], kg.nodes[7]);
    kg.add_edge(kg.nodes[6], kg.nodes[7]);
    kg.add_edge(kg.nodes[0], kg.nodes[8]);
    kg.add_edge(kg.nodes[7], kg.nodes[9]);
    kg.add_edge(kg.nodes[8], kg.nodes[9]);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;
    kg.setup_coverages(1);

    kg.nodes[4]->set_covg(4, 0, sample_id);
    kg.nodes[5]->set_covg(3, 0, sample_id);
    kg.nodes[6]->set_covg(5, 0, sample_id);
    kg.nodes[7]->set_covg(4, 0, sample_id);

    kg.num_reads = 5;
    kg.k = 3;

    std::vector<KmerNodePtr> mp;
    kg.set_p(0.01);
    kg.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kg.nodes[4]->set_covg(0, 0, sample_id);
    kg.nodes[5]->set_covg(0, 0, sample_id);
    kg.nodes[6]->set_covg(0, 0, sample_id);
    kg.nodes[7]->set_covg(0, 0, sample_id);
    kg.nodes[8]->set_covg(5, 1, sample_id);
    kg.set_p(0.01);
    kg.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    exp_order = {kg.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);
}

TEST(KmerGraphWithCoverageTest, findMaxPath2Level_nbin) {
    KmerGraph kg;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(8, 9), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(12, 13), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(16, 18), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 10;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[2], kg.nodes[3]);
    kg.add_edge(kg.nodes[0], kg.nodes[4]);
    kg.add_edge(kg.nodes[4], kg.nodes[5]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);
    kg.add_edge(kg.nodes[3], kg.nodes[7]);
    kg.add_edge(kg.nodes[6], kg.nodes[7]);
    kg.add_edge(kg.nodes[0], kg.nodes[8]);
    kg.add_edge(kg.nodes[7], kg.nodes[9]);
    kg.add_edge(kg.nodes[8], kg.nodes[9]);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;
    kg.setup_coverages(1);

    kg.nodes[4]->set_covg(4, 0, sample_id);
    kg.nodes[5]->set_covg(3, 0, sample_id);
    kg.nodes[6]->set_covg(5, 0, sample_id);
    kg.nodes[7]->set_covg(4, 0, sample_id);

    kg.num_reads = 5;
    kg.k = 3;

    std::vector<KmerNodePtr> mp;
    kg.set_p(0.01);
    kg.find_max_path(mp, "nbin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kg.nodes[4]->set_covg(0, 0, sample_id);
    kg.nodes[5]->set_covg(0, 0, sample_id);
    kg.nodes[6]->set_covg(0, 0, sample_id);
    kg.nodes[7]->set_covg(0, 0, sample_id);
    kg.nodes[8]->set_covg(5, 1, sample_id);
    kg.set_p(0.01);
    kg.find_max_path(mp, "nbin", max_num_kmers_to_average, sample_id);
    exp_order = {kg.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);
}

TEST(KmerGraphWithCoverageTest, findMaxPath2Level_lin) {
    KmerGraph kg;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(8, 9), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(12, 13), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(16, 18), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 10;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[2], kg.nodes[3]);
    kg.add_edge(kg.nodes[0], kg.nodes[4]);
    kg.add_edge(kg.nodes[4], kg.nodes[5]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);
    kg.add_edge(kg.nodes[3], kg.nodes[7]);
    kg.add_edge(kg.nodes[6], kg.nodes[7]);
    kg.add_edge(kg.nodes[0], kg.nodes[8]);
    kg.add_edge(kg.nodes[7], kg.nodes[9]);
    kg.add_edge(kg.nodes[8], kg.nodes[9]);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;
    kg.setup_coverages(1);

    kg.nodes[4]->set_covg(4, 0, sample_id);
    kg.nodes[5]->set_covg(3, 0, sample_id);
    kg.nodes[6]->set_covg(5, 0, sample_id);
    kg.nodes[7]->set_covg(4, 0, sample_id);

    kg.num_reads = 5;
    kg.k = 3;

    std::vector<KmerNodePtr> mp;
    kg.set_p(0.01);
    kg.find_max_path(mp, "lin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kg.nodes[4]->set_covg(0, 0, sample_id);
    kg.nodes[5]->set_covg(0, 0, sample_id);
    kg.nodes[6]->set_covg(0, 0, sample_id);
    kg.nodes[7]->set_covg(0, 0, sample_id);
    kg.nodes[8]->set_covg(5, 1, sample_id);
    kg.set_p(0.01);
    kg.find_max_path(mp, "lin", max_num_kmers_to_average, sample_id);
    exp_order = {kg.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);
}
*/
/*
TEST(KmerGraphWithCoverageTest, find_max_paths_2Level) {
    KmerGraph kg;
    deque<Interval> d = {Interval(0, 0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(8, 9), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(12, 13), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(16, 18), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 10;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[2], kg.nodes[3]);
    kg.add_edge(kg.nodes[0], kg.nodes[4]);
    kg.add_edge(kg.nodes[4], kg.nodes[5]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);
    kg.add_edge(kg.nodes[3], kg.nodes[7]);
    kg.add_edge(kg.nodes[6], kg.nodes[7]);
    kg.add_edge(kg.nodes[0], kg.nodes[8]);
    kg.add_edge(kg.nodes[7], kg.nodes[9]);
    kg.add_edge(kg.nodes[8], kg.nodes[9]);

    kg.nodes[4]->covg[0] += 4;
    kg.nodes[5]->covg[0] += 3;
    kg.nodes[6]->covg[0] += 5;
    kg.nodes[7]->covg[0] += 4;
    kg.nodes[8]->covg[1] += 5;

    kg.num_reads = 10;
    kg.k = 3;

    kg.set_p(0.01);
    uint32_t sample_id = 0;
    vector<vector<KmerNodePtr>> mps = kg.find_max_paths(2, sample_id);
    EXPECT_EQ((uint) 2, mps.size());
    vector<KmerNodePtr> exp_order = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mps[1]);

    exp_order = {kg.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mps[0]);
}
 */
/*
TEST(KmerGraphWithCoverageTest, random_paths) {
    KmerGraph kg;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(8, 9), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(12, 13), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(16, 18), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 10;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[2], kg.nodes[3]);
    kg.add_edge(kg.nodes[0], kg.nodes[4]);
    kg.add_edge(kg.nodes[4], kg.nodes[5]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);
    kg.add_edge(kg.nodes[3], kg.nodes[7]);
    kg.add_edge(kg.nodes[6], kg.nodes[7]);
    kg.add_edge(kg.nodes[0], kg.nodes[8]);
    kg.add_edge(kg.nodes[7], kg.nodes[9]);
    kg.add_edge(kg.nodes[8], kg.nodes[9]);

    vector<vector<KmerNodePtr>> rps;
    vector<KmerNodePtr> exp_order1 = {kg.nodes[1], kg.nodes[2], kg.nodes[3], kg.nodes[7]};
    vector<KmerNodePtr> exp_order2 = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7]};
    vector<KmerNodePtr> exp_order3 = {kg.nodes[8]};

    rps = kg.get_random_paths(10);
    for (uint i = 0; i != rps.size(); ++i) {
        for (uint j = 0; j != rps[i].size(); ++j) {
            if (rps[i][j]->id == 1) {
                EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order1, rps[i]);
            } else if (rps[i][j]->id == 4) {
                EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order2, rps[i]);
            } else if (rps[i][j]->id == 8) {
                EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order3, rps[i]);
            }
        }
    }
}

TEST(KmerGraphWithCoverageTest, path_prob) {
    KmerGraph kg;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(8, 9), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(12, 13), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(16, 18), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 10;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[2], kg.nodes[3]);
    kg.add_edge(kg.nodes[0], kg.nodes[4]);
    kg.add_edge(kg.nodes[4], kg.nodes[5]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);
    kg.add_edge(kg.nodes[3], kg.nodes[7]);
    kg.add_edge(kg.nodes[6], kg.nodes[7]);
    kg.add_edge(kg.nodes[0], kg.nodes[8]);
    kg.add_edge(kg.nodes[7], kg.nodes[9]);
    kg.add_edge(kg.nodes[8], kg.nodes[9]);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;
    kg.setup_coverages(1);

    kg.nodes[4]->set_covg(4, 0, sample_id);
    kg.nodes[5]->set_covg(3, 0, sample_id);
    kg.nodes[6]->set_covg(5, 0, sample_id);
    kg.nodes[7]->set_covg(4, 0, sample_id);

    kg.num_reads = 5;
    kg.k = 3;

    std::vector<KmerNodePtr> mp;
    kg.set_p(0.01);
    float mp_p = kg.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7], kg.nodes[9]};
    float exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kg.prob(exp_order[i]->id, sample_id);
    }
    exp_p /= 4;
    EXPECT_EQ(mp_p, exp_p);

    mp.clear();
    kg.nodes[4]->set_covg(0, 0, sample_id);
    kg.nodes[5]->set_covg(0, 0, sample_id);
    kg.nodes[6]->set_covg(0, 0, sample_id);
    kg.nodes[7]->set_covg(0, 0, sample_id);
    kg.nodes[8]->set_covg(5, 1, sample_id);
    kg.set_p(0.01);
    mp_p = kg.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    exp_order = {kg.nodes[8], kg.nodes[9]};
    exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kg.prob(exp_order[i]->id, sample_id);
    }
    EXPECT_EQ(mp_p, exp_p);
}
*/
/*
TEST(KmerGraphWithCoverageTest, path_probs) {
    KmerGraph kg;
    deque<Interval> d = {Interval(0, 0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(8, 9), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 17)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(12, 13), Interval(16, 18)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(16, 18), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    uint j = 10;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[2], kg.nodes[3]);
    kg.add_edge(kg.nodes[0], kg.nodes[4]);
    kg.add_edge(kg.nodes[4], kg.nodes[5]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);
    kg.add_edge(kg.nodes[3], kg.nodes[7]);
    kg.add_edge(kg.nodes[6], kg.nodes[7]);
    kg.add_edge(kg.nodes[0], kg.nodes[8]);
    kg.add_edge(kg.nodes[7], kg.nodes[9]);
    kg.add_edge(kg.nodes[8], kg.nodes[9]);

    kg.nodes[4]->covg[0] += 4;
    kg.nodes[5]->covg[0] += 3;
    kg.nodes[6]->covg[0] += 5;
    kg.nodes[7]->covg[0] += 4;
    kg.nodes[8]->covg[1] += 5;

    kg.num_reads = 10;
    kg.k = 3;
    kg.set_p(0.01);
    vector<vector<KmerNodePtr>> mps = kg.find_max_paths(2, <#initializer#>);
    EXPECT_EQ((uint) 2, mps.size());

    // check get right answer
    vector<KmerNodePtr> exp_nodes = {kg.nodes[4], kg.nodes[5], kg.nodes[6], kg.nodes[7], kg.nodes[8]};
    float exp_p = 0;
    for (uint i = 0; i != exp_nodes.size(); ++i) {
        exp_p += kg.prob(exp_nodes[i]->id, 5);
    }
    exp_p /= 5;
    EXPECT_EQ(kg.prob_paths(mps), exp_p);
}
 */
/*
TEST(KmerGraphWithCoverageTest, save_covg_dist) {
    KmerGraph kg;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p, p1, p2;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 3)};
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);
    d = {Interval(4, 4)};
    p.initialize(d);
    kg.add_node(p);
    kg.setup_coverages(1);
    kg.nodes[1]->set_covg(5, 1, 0);
    kg.nodes[2]->set_covg(4, 1, 0);
    kg.nodes[1]->num_AT = 4;
    kg.nodes[2]->num_AT = 6;

    kg.save_covg_dist("test_cases/kmergraph_test.covg.txt");
}

TEST(KmerGraphWithCoverageTest, save) {
    auto l = std::make_shared<LocalPRG>(LocalPRG(1, "test localPRG", "ACGT"));

    KmerGraphWithCoverage kg;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);
    EXPECT_EQ((uint) 0, kg.nodes[0]->num_AT);

    //kg.save("../test/test_cases/kmergraph_test.gfa");
    kg.save("kmergraph_test.gfa", l);
}

TEST(KmerGraphWithCoverageTest, save_no_prg) {
    KmerGraphWithCoverage kg;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);
    EXPECT_EQ((uint) 0, kg.nodes[0]->num_AT);

    //kg.save("../test/test_cases/kmergraph_test.gfa");
    kg.save("kmergraph_test2.gfa");
}

TEST(KmerGraphWithCoverageTest, load) {
    KmerGraphWithCoverage kg, read_kg;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);

    //read_kg.load("../test/test_cases/kmergraph_test.gfa");
    read_kg.load("kmergraph_test2.gfa");
    EXPECT_EQ(kg, read_kg);
}

TEST(KmerGraphWithCoverageTest, load_prg) {
    KmerGraphWithCoverage read_kg;
    EXPECT_DEATH(read_kg.load("kmergraph_test.gfa"), "");
}
*/