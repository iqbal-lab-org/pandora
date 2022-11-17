#include "gtest/gtest.h"
#include "test_macro.cpp"

#include "interval.h"
#include "prg/path.h"
#include "kmergraphwithcoverage.h"
#include "localPRG.h"
#include <stdint.h>
#include <iostream>
#include <cmath>
#include "test_helpers_containers.h"
#include "test_helpers.h"

using namespace prg;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// helper functions - used inside the tests to simplify them, and help with their
// readability and understanding e.g. inside a test, it is better for a reader to read
// one line: KmerGraph kmergraph = create_kmergraph(nb_of_nodes); than 5+ lines which
// achieves this
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// create a kmergraph with the given nb of nodes
KmerGraph create_kmergraph(uint32_t nb_of_nodes)
{
    KmerGraph kmergraph;

    for (uint32_t node_index = 0; node_index < nb_of_nodes; ++node_index) {
        std::deque<Interval> interval = { Interval(node_index, node_index) };
        prg::Path path;
        path.initialize(interval);
        kmergraph.add_node(path);
    }

    return kmergraph;
}

using Nodeindex_Strand_Sampleindex_Tuple
    = std::tuple<uint32_t, pandora::Strand, uint32_t>;

// function that will help to check the coverages in the tests
void check_coverages(const KmerGraphWithCoverage& kmergraph_with_coverage,
    std::map<Nodeindex_Strand_Sampleindex_Tuple, uint16_t>& expected_coverage,
    uint32_t nb_of_nodes, uint32_t nb_of_samples)
{
    // verifies all <node, strand, sample> coverage matches expected_coverage
    // if the tuple does not exist, then the expected_coverage is assumed to be 0 (as
    // per std::map::operator[] and uint16_t default constructor)
    for (uint32_t node_index = 0; node_index < nb_of_nodes; ++node_index) {
        for (uint32_t sample_index = 0; sample_index < nb_of_samples; ++sample_index) {
            EXPECT_EQ(
                kmergraph_with_coverage.get_reverse_covg(node_index, sample_index),
                expected_coverage[std::make_tuple(
                    node_index, pandora::Strand::Reverse, sample_index)]);
            EXPECT_EQ(
                kmergraph_with_coverage.get_forward_covg(node_index, sample_index),
                expected_coverage[std::make_tuple(
                    node_index, pandora::Strand::Forward, sample_index)]);
        }
    }
}

// set the coverage in both kmergraph_with_coverage and expected_coverage
void set_covg_helper(KmerGraphWithCoverage& kmergraph_with_coverage,
    std::map<Nodeindex_Strand_Sampleindex_Tuple, uint16_t>& expected_coverage,
    uint32_t node_id, uint16_t value, pandora::Strand strand, uint32_t sample_id)
{
    if (strand == pandora::Strand::Forward) {
        kmergraph_with_coverage.set_forward_covg(node_id, value, sample_id);
    } else {
        kmergraph_with_coverage.set_reverse_covg(node_id, value, sample_id);
    }

    expected_coverage[std::make_tuple(node_id, strand, sample_id)] = value;
}

// increment the coverage in both kmergraph_with_coverage and expected_coverage
void increment_covg_helper(KmerGraphWithCoverage& kmergraph_with_coverage,
    std::map<Nodeindex_Strand_Sampleindex_Tuple, uint16_t>& expected_coverage,
    uint32_t node_id, bool is_forward, uint32_t sample_id)
{
    auto old_covg { is_forward
            ? kmergraph_with_coverage.get_forward_covg(node_id, sample_id)
            : kmergraph_with_coverage.get_reverse_covg(node_id, sample_id) };

    if (is_forward) {
        kmergraph_with_coverage.increment_forward_covg(node_id, sample_id);
        expected_coverage[std::make_tuple(node_id, pandora::Strand::Forward, sample_id)]
            = old_covg + 1;
    } else {
        kmergraph_with_coverage.increment_reverse_covg(node_id, sample_id);
        expected_coverage[std::make_tuple(node_id, pandora::Strand::Reverse, sample_id)]
            = old_covg + 1;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST(KmerGraphWithCoverageTest, noCoverageSetAllCoverageShouldBe0)
{
    const uint32_t nb_of_nodes = 5;
    KmerGraph kmergraph = create_kmergraph(nb_of_nodes);
    const uint32_t nb_of_samples = 10;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph, nb_of_samples);
    std::map<Nodeindex_Strand_Sampleindex_Tuple, uint16_t> expected_coverage;

    // expect all coverages to be 0
    check_coverages(
        kmergraph_with_coverage, expected_coverage, nb_of_nodes, nb_of_samples);
}

TEST(KmerGraphWithCoverageTest, severalRamdomCoveragesSet)
{
    const uint32_t nb_of_nodes = 5;
    KmerGraph kmergraph = create_kmergraph(nb_of_nodes);
    const uint32_t nb_of_samples = 10;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph, nb_of_samples);
    std::map<Nodeindex_Strand_Sampleindex_Tuple, uint16_t> expected_coverage;

    // set some random coverages and verify if they are all correctly set
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 0, 100,
        pandora::Strand::Forward, 5);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 0, 150,
        pandora::Strand::Reverse, 5);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 1, 789,
        pandora::Strand::Reverse, 9);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 2, 120,
        pandora::Strand::Forward, 3);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 3, 130,
        pandora::Strand::Forward, 2);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 4, 780,
        pandora::Strand::Reverse, 7);

    check_coverages(
        kmergraph_with_coverage, expected_coverage, nb_of_nodes, nb_of_samples);
}

TEST(KmerGraphWithCoverageTest, maximumCoverageSet)
{
    const uint32_t nb_of_nodes = 5;
    KmerGraph kmergraph = create_kmergraph(nb_of_nodes);
    const uint32_t nb_of_samples = 10;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph, nb_of_samples);
    std::map<Nodeindex_Strand_Sampleindex_Tuple, uint16_t> expected_coverage;

    // setup the maximum coverage
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 4, UINT16_MAX,
        pandora::Strand::Forward, 1);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 4, UINT16_MAX,
        pandora::Strand::Reverse, 1);

    check_coverages(
        kmergraph_with_coverage, expected_coverage, nb_of_nodes, nb_of_samples);
}

TEST(KmerGraphWithCoverageTest, minimumCoverageSet)
{
    const uint32_t nb_of_nodes = 5;
    KmerGraph kmergraph = create_kmergraph(nb_of_nodes);
    const uint32_t nb_of_samples = 10;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph, nb_of_samples);
    std::map<Nodeindex_Strand_Sampleindex_Tuple, uint16_t> expected_coverage;

    // setup the minimum coverage
    set_covg_helper(
        kmergraph_with_coverage, expected_coverage, 3, 0, pandora::Strand::Forward, 1);
    set_covg_helper(
        kmergraph_with_coverage, expected_coverage, 3, 0, pandora::Strand::Reverse, 1);

    check_coverages(
        kmergraph_with_coverage, expected_coverage, nb_of_nodes, nb_of_samples);
}

TEST(KmerGraphWithCoverageTest, incrementCoverage)
{
    const uint32_t nb_of_nodes = 5;
    KmerGraph kmergraph = create_kmergraph(nb_of_nodes);
    const uint32_t nb_of_samples = 10;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph, nb_of_samples);
    std::map<Nodeindex_Strand_Sampleindex_Tuple, uint16_t> expected_coverage;

    // set some random coverages
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 0, 100,
        pandora::Strand::Forward, 5);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 0, 150,
        pandora::Strand::Reverse, 5);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 1, 789,
        pandora::Strand::Reverse, 9);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 2, 120,
        pandora::Strand::Forward, 3);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 3, 130,
        pandora::Strand::Forward, 2);
    set_covg_helper(kmergraph_with_coverage, expected_coverage, 4, 780,
        pandora::Strand::Reverse, 7);
    set_covg_helper(
        kmergraph_with_coverage, expected_coverage, 3, 0, pandora::Strand::Forward, 1);
    set_covg_helper(
        kmergraph_with_coverage, expected_coverage, 3, 0, pandora::Strand::Reverse, 1);
    // test increment coverage
    increment_covg_helper(kmergraph_with_coverage, expected_coverage, 0, true, 5);
    increment_covg_helper(kmergraph_with_coverage, expected_coverage, 0, false, 5);
    increment_covg_helper(kmergraph_with_coverage, expected_coverage, 1, false, 9);
    increment_covg_helper(kmergraph_with_coverage, expected_coverage, 2, true, 3);
    increment_covg_helper(kmergraph_with_coverage, expected_coverage, 3, true, 2);
    increment_covg_helper(kmergraph_with_coverage, expected_coverage, 4, false, 7);
    increment_covg_helper(kmergraph_with_coverage, expected_coverage, 3, true, 1);
    increment_covg_helper(kmergraph_with_coverage, expected_coverage, 3, false, 1);

    check_coverages(
        kmergraph_with_coverage, expected_coverage, nb_of_nodes, nb_of_samples);
}

TEST(KmerGraphWithCoverageTest, incrementCoverageAboveMaximumValue)
{
    const uint32_t nb_of_nodes = 5;
    KmerGraph kmergraph = create_kmergraph(nb_of_nodes);
    const uint32_t nb_of_samples = 10;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph, nb_of_samples);
    std::map<Nodeindex_Strand_Sampleindex_Tuple, uint16_t> expected_coverage;

    set_covg_helper(kmergraph_with_coverage, expected_coverage, 2, UINT16_MAX,
        pandora::Strand::Reverse, 9);
    check_coverages(
        kmergraph_with_coverage, expected_coverage, nb_of_nodes, nb_of_samples);
    kmergraph_with_coverage.increment_reverse_covg(2, 9);

    EXPECT_EQ(kmergraph_with_coverage.get_reverse_covg(2, 9), UINT16_MAX);
    check_coverages(
        kmergraph_with_coverage, expected_coverage, nb_of_nodes, nb_of_samples);
}

TEST(KmerGraphWithCoverageTest, set_exp_depth_covg)
{
    KmerGraph kmergraph;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.set_exp_depth_covg(10);
    EXPECT_EQ(kmergraph_with_coverage.exp_depth_covg, (uint)10);
}

TEST(KmerGraphWithCoverageTest, set_p)
{
    KmerGraph kmergraph;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    ASSERT_EXCEPTION(kmergraph_with_coverage.set_binomial_parameter_p(0.4),
        FatalRuntimeError, "Error setting binomial parameter p, invalid parameters");
    kmergraph_with_coverage.kmer_prg->k = 3;
    ASSERT_EXCEPTION(kmergraph_with_coverage.set_binomial_parameter_p(0),
        FatalRuntimeError, "Error setting binomial parameter p, invalid parameters");
    ASSERT_EXCEPTION(kmergraph_with_coverage.set_binomial_parameter_p(1),
        FatalRuntimeError, "Error setting binomial parameter p, invalid parameters");
    kmergraph_with_coverage.set_binomial_parameter_p(0.5);
    EXPECT_EQ(1 / exp(1.5) - 0.00001 <= kmergraph_with_coverage.binomial_parameter_p
            and 1 / exp(1.5) + 0.00001 >= kmergraph_with_coverage.binomial_parameter_p,
        true);
}

TEST(KmerGraphWithCoverageTest, set_nb)
{
    KmerGraph kmergraph;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.set_negative_binomial_parameters(0, 0);
    EXPECT_FLOAT_EQ(
        kmergraph_with_coverage.negative_binomial_parameter_p, 0.015); // unchanged
}

TEST(KmerGraphWithCoverageTest, prob_failNoNodes)
{
    uint32_t sample_id = 0;
    KmerGraph kmergraph;
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    ASSERT_EXCEPTION(kmergraph_with_coverage.bin_prob(0, sample_id), FatalRuntimeError,
        "Impossible to compute bin_prob, no reads were mapped to this kmer graph");
}

TEST(KmerGraphWithCoverageTest, prob_failNoP)
{
    uint32_t sample_id = 0;
    KmerGraph kmergraph;

    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    ASSERT_EXCEPTION(kmergraph_with_coverage.bin_prob(0, sample_id), FatalRuntimeError,
        "Impossible to compute bin_prob, no reads were mapped to this kmer graph");
}

TEST(KmerGraphWithCoverageTest, prob_failNoNumReads)
{
    uint32_t sample_id = 0;
    KmerGraph kmergraph;

    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.kmer_prg->k = 3;
    kmergraph_with_coverage.set_binomial_parameter_p(0.5);

    ASSERT_EXCEPTION(kmergraph_with_coverage.bin_prob(0, sample_id), FatalRuntimeError,
        "Impossible to compute bin_prob, no reads were mapped to this kmer graph");
}

TEST(KmerGraphWithCoverageTest, prob_simple)
{
    uint32_t sample_id = 0;
    KmerGraph kmergraph;

    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.kmer_prg->k = 3;
    kmergraph_with_coverage.set_binomial_parameter_p(0.5);
    kmergraph_with_coverage.num_reads = 1;

    EXPECT_EQ(kmergraph_with_coverage.kmer_prg->nodes.size(), (uint)1);
    EXPECT_EQ(0, kmergraph_with_coverage.bin_prob(0, sample_id));
}

TEST(KmerGraphWithCoverageTest, prob_realNodeCovgs)
{
    uint32_t sample_id = 0;
    KmerGraph kmergraph;

    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13) };
    p.initialize(d);
    kmergraph.add_node(p);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.kmer_prg->k = 3;
    kmergraph_with_coverage.set_binomial_parameter_p(0.5);
    kmergraph_with_coverage.num_reads = 1;

    EXPECT_EQ(kmergraph_with_coverage.kmer_prg->nodes.size(), (uint)3);

    EXPECT_EQ(kmergraph_with_coverage.bin_prob(1, sample_id),
        kmergraph_with_coverage.bin_prob(1, sample_id));
    EXPECT_EQ(kmergraph_with_coverage.bin_prob(2, sample_id),
        kmergraph_with_coverage.bin_prob(2, sample_id));
}

KmerGraph setup_simple_kmergraph()
{
    KmerGraph kmergraph;
    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(0, 1), Interval(19, 20), Interval(23, 24) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(24, 24) };
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

TEST(KmerGraphWithCoverageTest, findMaxPath_InvalidProbModel)
{
    KmerGraph kmergraph = setup_simple_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_forward_covg(1, 4, sample_id);
    kmergraph_with_coverage.set_forward_covg(2, 3, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    std::vector<KmerNodePtr> mp;
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    ASSERT_EXCEPTION(kmergraph_with_coverage.find_max_path(
                         mp, "exp", max_num_kmers_to_average, sample_id),
        FatalRuntimeError, "Invalid probability model for kmer coverage distribution");
}

TEST(KmerGraphWithCoverageTest, findMaxPathSimple)
{
    KmerGraph kmergraph = setup_simple_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_forward_covg(1, 4, sample_id);
    kmergraph_with_coverage.set_forward_covg(2, 3, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    std::vector<KmerNodePtr> mp;
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    kmergraph_with_coverage.find_max_path(
        mp, "bin", max_num_kmers_to_average, sample_id);
    std::vector<KmerNodePtr> exp_order = { kmergraph.nodes[1], kmergraph.nodes[2] };
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kmergraph_with_coverage.set_forward_covg(1, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(2, 0, sample_id);
    kmergraph_with_coverage.set_reverse_covg(5, 5, sample_id);
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    kmergraph_with_coverage.find_max_path(
        mp, "bin", max_num_kmers_to_average, sample_id);
    exp_order = { kmergraph.nodes[5] };
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mp);
}

TEST(KmerGraphWithCoverageTest, findMaxPathSimple_WithMaxKmersInAvg)
{
    KmerGraph kmergraph = setup_simple_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 1;

    kmergraph_with_coverage.set_forward_covg(1, 4, sample_id);
    kmergraph_with_coverage.set_forward_covg(2, 3, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    std::vector<KmerNodePtr> mp;
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    kmergraph_with_coverage.find_max_path(
        mp, "bin", max_num_kmers_to_average, sample_id);
    std::vector<KmerNodePtr> exp_order = { kmergraph.nodes[1], kmergraph.nodes[2] };
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mp);

    mp.clear();
    kmergraph_with_coverage.set_forward_covg(1, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(2, 0, sample_id);
    kmergraph_with_coverage.set_reverse_covg(5, 5, sample_id);
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    kmergraph_with_coverage.find_max_path(
        mp, "bin", max_num_kmers_to_average, sample_id);
    exp_order = { kmergraph.nodes[5] };
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mp);
}

KmerGraph setup_2level_kmergraph()
{
    KmerGraph kmergraph;
    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(8, 9), Interval(16, 18) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(4, 5), Interval(12, 13), Interval(16, 17) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(12, 13), Interval(16, 18) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(16, 18), Interval(23, 24) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(0, 1), Interval(19, 20), Interval(23, 24) };
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(24, 24) };
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

TEST(KmerGraphWithCoverageTest, findMaxPath2Level_bin)
{
    KmerGraph kmergraph = setup_2level_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_forward_covg(4, 4, sample_id);
    kmergraph_with_coverage.set_forward_covg(5, 3, sample_id);
    kmergraph_with_coverage.set_forward_covg(6, 4, sample_id);
    kmergraph_with_coverage.set_forward_covg(7, 3, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    std::vector<KmerNodePtr> mp;
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);

    auto mp_p = kmergraph_with_coverage.find_max_path(
        mp, "bin", max_num_kmers_to_average, sample_id);
    std::vector<KmerNodePtr> exp_order = { kmergraph.nodes[4], kmergraph.nodes[5],
        kmergraph.nodes[6], kmergraph.nodes[7] };
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mp);

    float exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.bin_prob(exp_order[i]->id, sample_id);
    }
    exp_p /= 4;
    EXPECT_EQ(mp_p, exp_p);

    mp.clear();
    kmergraph_with_coverage.set_forward_covg(4, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(5, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(6, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(7, 0, sample_id);
    kmergraph_with_coverage.set_reverse_covg(8, 5, sample_id);
    mp_p = kmergraph_with_coverage.find_max_path(
        mp, "bin", max_num_kmers_to_average, sample_id);
    exp_order = { kmergraph.nodes[8] };
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mp);

    exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.bin_prob(exp_order[i]->id, sample_id);
    }
    EXPECT_EQ(mp_p, exp_p);
}

TEST(KmerGraphWithCoverageTest, findMaxPath2Level_nbin)
{
    KmerGraph kmergraph = setup_2level_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_forward_covg(4, 4, sample_id);
    kmergraph_with_coverage.set_forward_covg(5, 3, sample_id);
    kmergraph_with_coverage.set_forward_covg(6, 4, sample_id);
    kmergraph_with_coverage.set_forward_covg(7, 3, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    std::vector<KmerNodePtr> mp;
    auto mp_p = kmergraph_with_coverage.find_max_path(
        mp, "nbin", max_num_kmers_to_average, sample_id);
    std::vector<KmerNodePtr> exp_order = { kmergraph.nodes[4], kmergraph.nodes[5],
        kmergraph.nodes[6], kmergraph.nodes[7] };
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mp);

    float exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.nbin_prob(exp_order[i]->id, sample_id);
    }
    exp_p /= 4;
    EXPECT_EQ(mp_p, exp_p);

    mp.clear();
    kmergraph_with_coverage.set_forward_covg(4, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(5, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(6, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(7, 0, sample_id);
    kmergraph_with_coverage.set_reverse_covg(8, 5, sample_id);
    mp_p = kmergraph_with_coverage.find_max_path(
        mp, "nbin", max_num_kmers_to_average, sample_id);
    exp_order = { kmergraph.nodes[8] };
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mp);

    exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.nbin_prob(exp_order[i]->id, sample_id);
    }
    EXPECT_EQ(mp_p, exp_p);
}

TEST(KmerGraphWithCoverageTest, findMaxPath2Level_lin)
{
    KmerGraph kmergraph = setup_2level_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    uint32_t sample_id = 0;
    uint32_t max_num_kmers_to_average = 100;

    kmergraph_with_coverage.set_forward_covg(4, 4, sample_id);
    kmergraph_with_coverage.set_forward_covg(5, 3, sample_id);
    kmergraph_with_coverage.set_forward_covg(6, 4, sample_id);
    kmergraph_with_coverage.set_forward_covg(7, 3, sample_id);

    kmergraph_with_coverage.num_reads = 5;
    kmergraph_with_coverage.kmer_prg->k = 3;

    std::vector<KmerNodePtr> mp;
    kmergraph_with_coverage.set_binomial_parameter_p(0.01);
    auto mp_p = kmergraph_with_coverage.find_max_path(
        mp, "lin", max_num_kmers_to_average, sample_id);
    std::vector<KmerNodePtr> exp_order = { kmergraph.nodes[4], kmergraph.nodes[5],
        kmergraph.nodes[6], kmergraph.nodes[7] };
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mp);

    float exp_p = 0;
    for (uint i = 0; i != exp_order.size(); ++i) {
        exp_p += kmergraph_with_coverage.lin_prob(exp_order[i]->id, sample_id);
    }
    exp_p /= 4;
    EXPECT_EQ(mp_p, exp_p);

    mp.clear();
    kmergraph_with_coverage.set_forward_covg(4, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(5, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(6, 0, sample_id);
    kmergraph_with_coverage.set_forward_covg(7, 0, sample_id);
    kmergraph_with_coverage.set_reverse_covg(8, 5, sample_id);
    mp_p = kmergraph_with_coverage.find_max_path(
        mp, "lin", max_num_kmers_to_average, sample_id);
    exp_order = { kmergraph.nodes[8] };
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mp);

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

    std::vector<std::vector<KmerNodePtr>> mps = kmergraph_with_coverage.find_max_paths(2,
sample_id); EXPECT_EQ((uint) 2, mps.size()); std::vector<KmerNodePtr> exp_order =
{kmergraph.nodes[4], kmergraph.nodes[5], kmergraph.nodes[6], kmergraph.nodes[7]};
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mps[1]);

    exp_order = {kmergraph.nodes[8]};
    EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order, mps[0]);
}
 */

TEST(KmerGraphWithCoverageTest, random_paths)
{
    KmerGraph kmergraph = setup_2level_kmergraph();
    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);

    std::vector<std::vector<KmerNodePtr>> rps;
    std::vector<KmerNodePtr> exp_order1 = { kmergraph.nodes[1], kmergraph.nodes[2],
        kmergraph.nodes[3], kmergraph.nodes[7] };
    std::vector<KmerNodePtr> exp_order2 = { kmergraph.nodes[4], kmergraph.nodes[5],
        kmergraph.nodes[6], kmergraph.nodes[7] };
    std::vector<KmerNodePtr> exp_order3 = { kmergraph.nodes[8] };

    rps = kmergraph_with_coverage.get_random_paths(10);
    for (uint i = 0; i != rps.size(); ++i) {
        for (uint j = 0; j != rps[i].size(); ++j) {
            if (rps[i][j]->id == 1) {
                EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order1, rps[i]);
            } else if (rps[i][j]->id == 4) {
                EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order2, rps[i]);
            } else if (rps[i][j]->id == 8) {
                EXPECT_ITERABLE_EQ(std::vector<KmerNodePtr>, exp_order3, rps[i]);
            }
        }
    }
}

TEST(KmerGraphWithCoverageTest, save_covg_dist)
{
    KmerGraph kmergraph;
    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p, p1, p2;
    p.initialize(d);
    kmergraph.add_node(p);
    d = { Interval(0, 3) };
    p1.initialize(d);
    auto n1 = kmergraph.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kmergraph.add_node(p2);
    kmergraph.add_edge(n1, n2);
    d = { Interval(4, 4) };
    p.initialize(d);
    kmergraph.add_node(p);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.set_reverse_covg(1, 5, 0);
    kmergraph_with_coverage.set_reverse_covg(2, 4, 0);

    kmergraph_with_coverage.save_covg_dist("test_cases/kmergraph_test.covg.txt");
}

TEST(KmerGraphWithCoverageTest, save_no_prg)
{
    KmerGraph kmergraph;
    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kmergraph.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kmergraph.add_node(p2);
    kmergraph.add_edge(n1, n2);
    EXPECT_EQ((uint)0, kmergraph.nodes[0]->num_AT);

    KmerGraphWithCoverage kmergraph_with_coverage(&kmergraph);
    kmergraph_with_coverage.set_forward_covg(0, 4, 0);
    kmergraph_with_coverage.set_reverse_covg(1, 5, 0);

    const fs::path filepath { std::tmpnam(nullptr) };
    kmergraph_with_coverage.save(filepath);
    fs::ifstream ifs(filepath);
    const std::string actual_content(
        (std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));

    const std::string expected_content { "H\tVN:Z:1.0\tbn:Z:--linear --singlearr\n"
                                         "S\t0\t1{[0, 3)}\tFC:i:4\tRC:i:0\n"
                                         "L\t0\t+\t1\t+\t0M\n"
                                         "S\t1\t1{[1, 4)}\tFC:i:0\tRC:i:5\n" };

    EXPECT_EQ(actual_content, expected_content);
    ASSERT_TRUE(fs::remove(filepath));
}

TEST(KmerGraphWithCoverageTest, load)
{
    KmerGraph kmergraph;
    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kmergraph.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kmergraph.add_node(p2);
    kmergraph.add_edge(n1, n2);

    KmerGraphWithCoverage read_kmergraph_with_coverage(&kmergraph);
    const fs::path filepath { "deleteme.gfa" };
    {
        fs::ofstream outstream(filepath);
        outstream << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << std::endl;
        outstream << "S\t0\t1{[0, 3)}\tFC:i:4\tRC:i:0" << std::endl;
        outstream << "L\t0\t+\t1\t+\t0M" << std::endl;
        outstream << "S\t1\t1{[1, 4)}\tFC:i:0\tRC:i:5" << std::endl;
    }
    read_kmergraph_with_coverage.load(filepath.string());
    fs::remove(filepath);
    EXPECT_EQ(read_kmergraph_with_coverage.get_reverse_covg(1, 0), 5);
    EXPECT_EQ(read_kmergraph_with_coverage.get_forward_covg(1, 0), 0);
    EXPECT_EQ(read_kmergraph_with_coverage.get_forward_covg(0, 0), 4);
    EXPECT_EQ(read_kmergraph_with_coverage.get_reverse_covg(0, 0), 0);
}

TEST(KmerGraphWithCoverageTest, load_prg)
{
    KmerGraph kmergraph;
    KmerGraphWithCoverage read_kmergraph_with_coverage(&kmergraph);
    ASSERT_EXCEPTION(read_kmergraph_with_coverage.load("kmergraph_test.gfa"),
        FatalRuntimeError, "Error reading GFA");
}
