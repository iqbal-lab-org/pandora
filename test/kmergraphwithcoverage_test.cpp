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
#include <cstdint>

using namespace prg;

/*TEST(KmerGraphWithCoverageTest, equals) {
    KmerGraph kmergraph1, kmergraph2;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kmergraph1.add_node(p1);
    auto m1 = kmergraph2.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kmergraph1.add_node(p2);
    auto m2 = kmergraph2.add_node(p2);
    kmergraph1.add_edge(n1, n2);
    kmergraph2.add_edge(m1, m2);

    d = {Interval(2, 5)};
    p3.initialize(d);
    auto m3 = kmergraph2.add_node(p3);

    // same as themselves, different if different numbers of nodes
    EXPECT_EQ(kmergraph1, kmergraph1);
    EXPECT_EQ(kmergraph2, kmergraph2);
    EXPECT_EQ((kmergraph1 == kmergraph2), false);
    EXPECT_EQ((kmergraph2 == kmergraph1), false);

    auto n3 = kmergraph1.add_node(p3);
    kmergraph2.add_edge(m1, m3);

    // same as themselves, different if different numbers of edges
    EXPECT_EQ(kmergraph1, kmergraph1);
    EXPECT_EQ(kmergraph2, kmergraph2);
    EXPECT_EQ((kmergraph1 == kmergraph2), false);
    EXPECT_EQ((kmergraph2 == kmergraph1), false);

    kmergraph1.add_edge(n2, n3);

    // same as themselves, different if edges in different places
    EXPECT_EQ(kmergraph1, kmergraph1);
    EXPECT_EQ(kmergraph2, kmergraph2);
    EXPECT_EQ((kmergraph1 == kmergraph2), false);
    EXPECT_EQ((kmergraph2 == kmergraph1), false);
}

TEST(KmerGraphWithCoverageTest, copy) {
    KmerGraphWithCoverage kmergraph1;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kmergraph1.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kmergraph1.add_node(p2);
    kmergraph1.add_edge(n1, n2);

    KmerGraphWithCoverage kmergraph2(kmergraph1);

    EXPECT_EQ(kmergraph1, kmergraph1);
    EXPECT_EQ(kmergraph2, kmergraph2);
}

TEST(KmerGraphWithCoverageTest, assign) {
    KmerGraphWithCoverage kmergraph1;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    auto n = kmergraph1.add_node(p);
    d = {Interval(0, 3)};
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kmergraph1.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kmergraph1.add_node(p2);
    kmergraph1.add_edge(n1, n2);
    d = {Interval(11, 14)};
    p3.initialize(d);
    auto n3 = kmergraph1.add_node(p3);
    kmergraph1.add_edge(n1, n3);
    d = {Interval(15, 18)};
    p1.initialize(d);
    n1 = kmergraph1.add_node(p1);
    kmergraph1.add_edge(n2, n1);
    d = {Interval(20, 20)};
    p.initialize(d);
    n = kmergraph1.add_node(p);
    kmergraph1.add_edge(n1, n);
    kmergraph1.add_edge(n3, n);

    KmerGraphWithCoverage kmergraph2 = kmergraph1;

    EXPECT_EQ(kmergraph1, kmergraph1);
    EXPECT_EQ(kmergraph2, kmergraph2);
}*/

TEST(KmerGraphWithCoverageTest, get_and_set_increment_covg){
    //create a kmergraph with 5 nodes
    const uint32_t nb_of_nodes = 5;
    KmerGraph kmergraph;

    auto create_kmergraph_helper = [&](uint32_t node_index) {
        std::deque<Interval> interval = {Interval(node_index, node_index)};
        prg::Path path;
        path.initialize(interval);
        kmergraph.add_node(path);
    };

    for (uint32_t node_index=0; node_index < nb_of_nodes; ++node_index)
        create_kmergraph_helper(node_index);

    //create a KmerGraphWithCoverage on 10 samples
    const uint32_t nb_of_samples = 10;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph, nb_of_samples);

    //function that will help to check the coverages in the tests
    using Nodeindex_Strand_Sampleindex_Tuple = std::tuple<uint32_t, bool, uint32_t>;
    std::map<Nodeindex_Strand_Sampleindex_Tuple, uint16_t> expected_coverage;
    auto check_coverages = [&]() {
        //verifies all <node, strand, sample> coverage matches expected_coverage
        //if the tuple does not exist, then the expected_coverage is assumed to be 0 (as per std::map::operator[] and uint16_t default constructor)
        for (uint32_t node_index=0; node_index < nb_of_nodes; ++node_index) {
            for (uint32_t sample_index=0; sample_index<nb_of_samples; ++sample_index) {
                EXPECT_EQ(kmergraph_with_coverage.get_covg(node_index, 0, sample_index), expected_coverage[make_tuple(node_index, 0, sample_index)]);
                EXPECT_EQ(kmergraph_with_coverage.get_covg(node_index, 1, sample_index), expected_coverage[make_tuple(node_index, 1, sample_index)]);
            }
        }
    };

    //1. expect all coverages to be 0
    check_coverages();

    //2. lets set some random coverages
    auto set_covg_helper = [&] (uint32_t node_id, uint16_t value, bool strand, uint32_t sample_id) {
        kmergraph_with_coverage.set_covg(node_id, value, strand, sample_id);
        expected_coverage[make_tuple(node_id, strand, sample_id)] = value;
    };
    set_covg_helper(0,100,0,5);
    set_covg_helper(0,150,1,5);
    set_covg_helper(1,789,1,9);
    set_covg_helper(2,120,0,3);
    set_covg_helper(3,130,0,2);
    set_covg_helper(4,780,1,7);
    check_coverages();

    //3. setup the maximum coverage
    set_covg_helper(4,UINT16_MAX,0,1);
    set_covg_helper(4,UINT16_MAX,1,1);
    check_coverages();

    //4. setup the minimum coverage
    set_covg_helper(3,0,0,1);
    set_covg_helper(3,0,1,1);
    check_coverages();

    //5. test increment coverage
    auto increment_covg_helper = [&] (uint32_t node_id, bool strand, uint32_t sample_id) {
        auto old_covg = kmergraph_with_coverage.get_covg(node_id, strand, sample_id);
        kmergraph_with_coverage.increment_covg(node_id, strand, sample_id);
        expected_coverage[make_tuple(node_id, strand, sample_id)] = old_covg+1;
    };
    increment_covg_helper(0,0,5);
    increment_covg_helper(0,1,5);
    increment_covg_helper(1,1,9);
    increment_covg_helper(2,0,3);
    increment_covg_helper(3,0,2);
    increment_covg_helper(4,1,7);
    increment_covg_helper(3,0,1);
    increment_covg_helper(3,1,1);
    check_coverages();


    //6. test incrementing coverage above the maximum value
    set_covg_helper(2,UINT16_MAX,0,9); //first set to max value
    check_coverages();
    kmergraph_with_coverage.increment_covg(2, 0, 9); //now try to increment
    EXPECT_EQ(kmergraph_with_coverage.get_covg(2, 0, 9), UINT16_MAX); //should still be UINT16_MAX
    check_coverages(); //to be sure nothing changed
}


TEST(KmerGraphWithCoverageTest, set_exp_depth_covg){
    KmerGraph kmergraph;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.set_exp_depth_covg(10);
    EXPECT_EQ(kmergraph_with_coverage.exp_depth_covg, (uint)10);
}

TEST(KmerGraphWithCoverageTest, set_p) {
    KmerGraph kmergraph;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    EXPECT_DEATH(kmergraph_with_coverage.set_binomial_parameter_p(0.4), "");
    kmergraph_with_coverage.kmer_prg->k = 3;
    EXPECT_DEATH(kmergraph_with_coverage.set_binomial_parameter_p(0), "");
    EXPECT_DEATH(kmergraph_with_coverage.set_binomial_parameter_p(1), "");
    kmergraph_with_coverage.set_binomial_parameter_p(0.5);
    EXPECT_EQ(1 / exp(1.5) - 0.00001 <= kmergraph_with_coverage.binomial_parameter_p and 1 / exp(1.5) + 0.00001 >= kmergraph_with_coverage.binomial_parameter_p, true);
}

TEST(KmerGraphWithCoverageTest, set_nb) {
    KmerGraph kmergraph;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.set_negative_binomial_parameters(0, 0);
    EXPECT_FLOAT_EQ(kmergraph_with_coverage.negative_binomial_parameter_p, 0.015); //unchanged
}

TEST(KmerGraphWithCoverageTest, prob_failNoNodes) {
    uint32_t sample_id = 0;
    KmerGraph kmergraph;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    EXPECT_DEATH(kmergraph_with_coverage.bin_prob(0, sample_id), "");
}

TEST(KmerGraphWithCoverageTest, prob_failNoP) {
    uint32_t sample_id = 0;
    KmerGraph kmergraph;

    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    EXPECT_DEATH(kmergraph_with_coverage.bin_prob(0, sample_id), "");
}

TEST(KmerGraphWithCoverageTest, prob_failNoNumReads) {
    uint32_t sample_id = 0;
    KmerGraph kmergraph;

    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.kmer_prg->k = 3;
    kmergraph_with_coverage.set_binomial_parameter_p(0.5);

    EXPECT_DEATH(kmergraph_with_coverage.bin_prob(0, sample_id), "");
}

TEST(KmerGraphWithCoverageTest, prob_simple) {
    uint32_t sample_id = 0;
    KmerGraph kmergraph;

    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.kmer_prg->k = 3;
    kmergraph_with_coverage.set_binomial_parameter_p(0.5);
    kmergraph_with_coverage.num_reads = 1;

    EXPECT_EQ(kmergraph_with_coverage.kmer_prg->nodes.size(), (uint) 1);
    EXPECT_EQ(0, kmergraph_with_coverage.bin_prob(0, sample_id));
}

TEST(KmerGraphWithCoverageTest, prob_realNodeCovgs) {
    uint32_t sample_id = 0;
    KmerGraph kmergraph;

    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kmergraph.add_node(p);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.kmer_prg->k = 3;
    kmergraph_with_coverage.set_binomial_parameter_p(0.5);
    kmergraph_with_coverage.num_reads = 1;

    EXPECT_EQ(kmergraph_with_coverage.kmer_prg->nodes.size(), (uint) 3);

    EXPECT_EQ(kmergraph_with_coverage.bin_prob(1, sample_id), kmergraph_with_coverage.bin_prob(1, sample_id));
    EXPECT_EQ(kmergraph_with_coverage.bin_prob(2, sample_id), kmergraph_with_coverage.bin_prob(2, sample_id));

}

KmerGraph setup_simple_kmergraph(){
    KmerGraph kmergraph;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kmergraph.add_node(p);

    assert(kmergraph.nodes.size() == 7);

    kmergraph.add_edge(kmergraph.nodes[0], kmergraph.nodes[1]);
    kmergraph.add_edge(kmergraph.nodes[1], kmergraph.nodes[2]);
    kmergraph.add_edge(kmergraph.nodes[0], kmergraph.nodes[3]);
    kmergraph.add_edge(kmergraph.nodes[3], kmergraph.nodes[4]);
    kmergraph.add_edge(kmergraph.nodes[0], kmergraph.nodes[5]);
    kmergraph.add_edge(kmergraph.nodes[2], kmergraph.nodes[6]);
    kmergraph.add_edge(kmergraph.nodes[4], kmergraph.nodes[6]);
    kmergraph.add_edge(kmergraph.nodes[5], kmergraph.nodes[6]);
    return kmergraph;
}

TEST(KmerGraphWithCoverageTest, findMaxPath_InvalidProbModel) {
    KmerGraph kmergraph = setup_simple_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_covg(1,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(2,3, 0, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    vector<KmerNodePtr> mp;
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    EXPECT_DEATH(kmergraph_with_coverage.find_max_path(mp, "exp", max_num_kmers_to_average, sample_id),"");
}

TEST(KmerGraphWithCoverageTest, findMaxPathSimple) {
    KmerGraph kmergraph = setup_simple_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_covg(1,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(2,3, 0, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    vector<KmerNodePtr> mp;
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    kmergraph_with_coverage.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kmergraph.nodes[1], kmergraph.nodes[2]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kmergraph_with_coverage.set_covg(1, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(2, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(5, 5, 1, sample_id);
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    kmergraph_with_coverage.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    exp_order = {kmergraph.nodes[5]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);
}

TEST(KmerGraphWithCoverageTest, findMaxPathSimple_WithMaxKmersInAvg) {
    KmerGraph kmergraph = setup_simple_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 1;

    kmergraph_with_coverage.set_covg(1,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(2,3, 0, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    vector<KmerNodePtr> mp;
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    kmergraph_with_coverage.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kmergraph.nodes[1], kmergraph.nodes[2]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kmergraph_with_coverage.set_covg(1, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(2, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(5, 5, 1, sample_id);
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    kmergraph_with_coverage.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    exp_order = {kmergraph.nodes[5]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);
}

KmerGraph setup_2level_kmergraph() {
    KmerGraph kmergraph;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(8, 9), Interval(16, 18)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 17)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(12, 13), Interval(16, 18)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(16, 18), Interval(23, 24)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kmergraph.add_node(p);

    assert(kmergraph.nodes.size() == 10);

    kmergraph.add_edge(kmergraph.nodes[0], kmergraph.nodes[1]);
    kmergraph.add_edge(kmergraph.nodes[1], kmergraph.nodes[2]);
    kmergraph.add_edge(kmergraph.nodes[2], kmergraph.nodes[3]);
    kmergraph.add_edge(kmergraph.nodes[0], kmergraph.nodes[4]);
    kmergraph.add_edge(kmergraph.nodes[4], kmergraph.nodes[5]);
    kmergraph.add_edge(kmergraph.nodes[5], kmergraph.nodes[6]);
    kmergraph.add_edge(kmergraph.nodes[3], kmergraph.nodes[7]);
    kmergraph.add_edge(kmergraph.nodes[6], kmergraph.nodes[7]);
    kmergraph.add_edge(kmergraph.nodes[0], kmergraph.nodes[8]);
    kmergraph.add_edge(kmergraph.nodes[7], kmergraph.nodes[9]);
    kmergraph.add_edge(kmergraph.nodes[8], kmergraph.nodes[9]);
    return kmergraph;
}

TEST(KmerGraphWithCoverageTest, findMaxPath2Level_bin) {
    KmerGraph kmergraph = setup_2level_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_covg(4,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(5,3, 0, sample_id);
    kmergraph_with_coverage.set_covg(6,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(7,3, 0, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    std::vector<KmerNodePtr> mp;
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);

    auto mp_p = kmergraph_with_coverage.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kmergraph.nodes[4], kmergraph.nodes[5], kmergraph.nodes[6], kmergraph.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    float exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.bin_prob(exp_order[i]->id, sample_id);
    }
    exp_p /= 4;
    EXPECT_EQ(mp_p, exp_p);

    mp.clear();
    kmergraph_with_coverage.set_covg(4, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(5, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(6, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(7, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(8, 5, 1, sample_id);
    mp_p = kmergraph_with_coverage.find_max_path(mp, "bin", max_num_kmers_to_average, sample_id);
    exp_order = {kmergraph.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.bin_prob(exp_order[i]->id, sample_id);
    }
    EXPECT_EQ(mp_p, exp_p);
}

TEST(KmerGraphWithCoverageTest, findMaxPath2Level_nbin) {
    KmerGraph kmergraph = setup_2level_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_covg(4,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(5,3, 0, sample_id);
    kmergraph_with_coverage.set_covg(6,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(7,3, 0, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    std::vector<KmerNodePtr> mp;
    auto mp_p = kmergraph_with_coverage.find_max_path(mp, "nbin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kmergraph.nodes[4], kmergraph.nodes[5], kmergraph.nodes[6], kmergraph.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    float exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.nbin_prob(exp_order[i]->id, sample_id);
    }
    exp_p /= 4;
    EXPECT_EQ(mp_p, exp_p);

    mp.clear();
    kmergraph_with_coverage.set_covg(4, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(5, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(6, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(7, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(8, 5, 1, sample_id);
    mp_p = kmergraph_with_coverage.find_max_path(mp, "nbin", max_num_kmers_to_average, sample_id);
    exp_order = {kmergraph.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.nbin_prob(exp_order[i]->id, sample_id);
    }
    EXPECT_EQ(mp_p, exp_p);
}

TEST(KmerGraphWithCoverageTest, findMaxPath2Level_lin) {
    KmerGraph kmergraph = setup_2level_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_covg(4,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(5,3, 0, sample_id);
    kmergraph_with_coverage.set_covg(6,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(7,3, 0, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    std::vector<KmerNodePtr> mp;
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    auto mp_p = kmergraph_with_coverage.find_max_path(mp, "lin", max_num_kmers_to_average, sample_id);
    vector<KmerNodePtr> exp_order = {kmergraph.nodes[4], kmergraph.nodes[5], kmergraph.nodes[6], kmergraph.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    float exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.lin_prob(exp_order[i]->id, sample_id);
    }
    exp_p /= 4;
    EXPECT_EQ(mp_p, exp_p);

    mp.clear();
    kmergraph_with_coverage.set_covg(4, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(5, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(6, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(7, 0, 0, sample_id);
    kmergraph_with_coverage.set_covg(8, 5, 1, sample_id);
    mp_p = kmergraph_with_coverage.find_max_path(mp, "lin", max_num_kmers_to_average, sample_id);
    exp_order = {kmergraph.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mp);

    exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.lin_prob(exp_order[i]->id, sample_id);
    }
    EXPECT_EQ(mp_p, exp_p);
}

/*
TEST(KmerGraphWithCoverageTest, find_max_paths_2Level) {
    KmerGraph kmergraph = setup_2level_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_covg(4,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(5,3, 0, sample_id);
    kmergraph_with_coverage.set_covg(6,5, 0, sample_id);
    kmergraph_with_coverage.set_covg(7,4, 0, sample_id);
    kmergraph_with_coverage.set_covg(8,5, 1, sample_id);

    kmergraph_with_coverage.num_reads = 10;
    kmergraph_with_coverage.kmer_prg->k = 3;
    kmergraph_with_coverage.set_p(0.01);

    vector<vector<KmerNodePtr>> mps = kmergraph_with_coverage.find_max_paths(2, sample_id);
    EXPECT_EQ((uint) 2, mps.size());
    vector<KmerNodePtr> exp_order = {kmergraph.nodes[4], kmergraph.nodes[5], kmergraph.nodes[6], kmergraph.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mps[1]);

    exp_order = {kmergraph.nodes[8]};
    EXPECT_ITERABLE_EQ(vector<KmerNodePtr>, exp_order, mps[0]);
}
 */

TEST(KmerGraphWithCoverageTest, random_paths) {
    KmerGraph kmergraph = setup_2level_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    vector<vector<KmerNodePtr>> rps;
    vector<KmerNodePtr> exp_order1 = {kmergraph.nodes[1], kmergraph.nodes[2], kmergraph.nodes[3], kmergraph.nodes[7]};
    vector<KmerNodePtr> exp_order2 = {kmergraph.nodes[4], kmergraph.nodes[5], kmergraph.nodes[6], kmergraph.nodes[7]};
    vector<KmerNodePtr> exp_order3 = {kmergraph.nodes[8]};

    rps = kmergraph_with_coverage.get_random_paths(10);
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

/*
TEST(KmerGraphWithCoverageTest, path_probs) {
    KmerGraph kmergraph;
    deque<Interval> d = {Interval(0, 0)};
    Path p;
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(8, 9), Interval(16, 18)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 17)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(12, 13), Interval(16, 18)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(16, 18), Interval(23, 24)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kmergraph.add_node(p);
    uint j = 10;
    EXPECT_EQ(j, kmergraph.nodes.size());

    kmergraph.add_edge(kmergraph.nodes[0], kmergraph.nodes[1]);
    kmergraph.add_edge(kmergraph.nodes[1], kmergraph.nodes[2]);
    kmergraph.add_edge(kmergraph.nodes[2], kmergraph.nodes[3]);
    kmergraph.add_edge(kmergraph.nodes[0], kmergraph.nodes[4]);
    kmergraph.add_edge(kmergraph.nodes[4], kmergraph.nodes[5]);
    kmergraph.add_edge(kmergraph.nodes[5], kmergraph.nodes[6]);
    kmergraph.add_edge(kmergraph.nodes[3], kmergraph.nodes[7]);
    kmergraph.add_edge(kmergraph.nodes[6], kmergraph.nodes[7]);
    kmergraph.add_edge(kmergraph.nodes[0], kmergraph.nodes[8]);
    kmergraph.add_edge(kmergraph.nodes[7], kmergraph.nodes[9]);
    kmergraph.add_edge(kmergraph.nodes[8], kmergraph.nodes[9]);

    kmergraph.nodes[4]->covg[0] += 4;
    kmergraph.nodes[5]->covg[0] += 3;
    kmergraph.nodes[6]->covg[0] += 5;
    kmergraph.nodes[7]->covg[0] += 4;
    kmergraph.nodes[8]->covg[1] += 5;

    kmergraph.num_reads = 10;
    kmergraph.k = 3;
    kmergraph.set_p(0.01);
    vector<vector<KmerNodePtr>> mps = kmergraph.find_max_paths(2, <#initializer#>);
    EXPECT_EQ((uint) 2, mps.size());

    // check get right answer
    vector<KmerNodePtr> exp_nodes = {kmergraph.nodes[4], kmergraph.nodes[5], kmergraph.nodes[6], kmergraph.nodes[7], kmergraph.nodes[8]};
    float exp_p = 0;
    for (uint i = 0; i != exp_nodes.size(); ++i) {
        exp_p += kmergraph.prob(exp_nodes[i]->id, 5);
    }
    exp_p /= 5;
    EXPECT_EQ(kmergraph.prob_paths(mps), exp_p);
}
 */

TEST(KmerGraphWithCoverageTest, save_covg_dist) {
    KmerGraph kmergraph;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p, p1, p2;
    p.initialize(d);
    kmergraph.add_node(p);
    d = {Interval(0, 3)};
    p1.initialize(d);
    auto n1 = kmergraph.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kmergraph.add_node(p2);
    kmergraph.add_edge(n1, n2);
    d = {Interval(4, 4)};
    p.initialize(d);
    kmergraph.add_node(p);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.set_covg(1, 5, 1, 0);
    kmergraph_with_coverage.set_covg(2, 4, 1, 0);

    kmergraph_with_coverage.save_covg_dist("test_cases/kmergraph_test.covg.txt");
}

TEST(KmerGraphWithCoverageTest, save_no_prg) {
    KmerGraph kmergraph;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kmergraph.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kmergraph.add_node(p2);
    kmergraph.add_edge(n1, n2);
    EXPECT_EQ((uint) 0, kmergraph.nodes[0]->num_AT);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.set_covg(0, 4, 1, 0);
    kmergraph_with_coverage.set_covg(1, 5, 0, 0);

    kmergraph_with_coverage.save("kmergraphwithcoverage_test.gfa");
}

TEST(KmerGraphWithCoverageTest, load) {
    KmerGraph kmergraph;
    deque<Interval> d = {Interval(0, 3)};
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kmergraph.add_node(p1);
    d = {Interval(1, 4)};
    p2.initialize(d);
    auto n2 = kmergraph.add_node(p2);
    kmergraph.add_edge(n1, n2);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.set_covg(1, 5, 1, 0);
    kmergraph_with_coverage.set_covg(0, 4, 1, 0);

    KmerGraphWithCoverage read_kmergraph_with_coverage(&kmergraph);
    read_kmergraph_with_coverage.load("kmergraphwithcoverage_test.gfa");
    //EXPECT_EQ(kmergraph_with_coverage, read_kmergraph_with_coverage);
}

TEST(KmerGraphWithCoverageTest, load_prg) {
    KmerGraph kmergraph;
    KmerGraphWithCoverage read_kmergraph_with_coverage(&kmergraph);
    EXPECT_DEATH(read_kmergraph_with_coverage.load("kmergraph_test.gfa"), "");
}
