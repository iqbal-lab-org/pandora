#include "../test_helpers.h"
#include "../test_macro.cpp"
#include "denovo_discovery/candidate_region.h"
#include "fastaq.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <iostream>

TEST(CandidateRegionGetIntervalTest, noPaddingReturnsOriginalInterval)
{
    const CandidateRegion candidate { Interval(0, 2), "test" };

    const auto actual { candidate.get_interval() };
    const Interval expected { 0, 2 };

    EXPECT_EQ(actual, expected);
}

TEST(CandidateRegionGetIntervalTest,
    withPaddingLessThanIntervalStartReturnsOriginalIntervalWithPadding)
{
    const CandidateRegion candidate { Interval(3, 4), "test", 1 };

    const auto actual { candidate.get_interval() };
    const Interval expected { 2, 5 };

    EXPECT_EQ(actual, expected);
}

TEST(CandidateRegionGetIntervalTest,
    withPaddingGreaterThanIntervalStartReturnsOriginalIntervalWithPadding)
{
    const CandidateRegion candidate { Interval(3, 4), "test", 5 };

    const auto actual { candidate.get_interval() };
    const Interval expected { 0, 9 };

    EXPECT_EQ(actual, expected);
}

TEST(CandidateRegionGetNameTest, testCorrectValueRetrieved)
{
    const CandidateRegion candidate { Interval(0, 2), "test" };

    const auto& actual { candidate.get_name() };
    const std::string expected { "test" };

    EXPECT_EQ(actual, expected);
}

TEST(CandidateRegionGetIdTest, noPaddingTestCorrectIdRetrieved)
{
    const CandidateRegion candidate { Interval(0, 2), "test" };

    const auto actual { candidate.get_id() };
    const CandidateRegionIdentifier expected { candidate.get_interval(),
        candidate.get_name() };

    EXPECT_EQ(actual, expected);
}

TEST(CandidateRegionGetIdTest, withPaddingTestCorrectIdRetrieved)
{
    const CandidateRegion candidate { Interval(0, 2), "test", 6 };

    const auto actual { candidate.get_id() };
    const CandidateRegionIdentifier expected { candidate.get_interval(),
        candidate.get_name() };

    EXPECT_EQ(actual, expected);
}

TEST(GetMaxLikelihoodSequenceWithFlanksTest, noSequencesReturnEmptyString)
{
    const CandidateRegion candidate { Interval(0, 2), "test", 6 };

    const auto actual { candidate.get_max_likelihood_sequence_with_flanks() };
    const std::string expected;

    EXPECT_EQ(actual, expected);
}

TEST(GetMaxLikelihoodSequenceWithFlanksTest,
    noLeftFlankSequencesReturnMaxAndRightSequence)
{
    CandidateRegion candidate { Interval(0, 2), "test", 6 };
    candidate.max_likelihood_sequence = "max";
    candidate.right_flanking_sequence = "right";

    const auto actual { candidate.get_max_likelihood_sequence_with_flanks() };
    const std::string expected { "maxright" };

    EXPECT_EQ(actual, expected);
}

TEST(GetMaxLikelihoodSequenceWithFlanksTest,
    noRightFlankSequencesReturnLeftAndMaxSequence)
{
    CandidateRegion candidate { Interval(0, 2), "test", 6 };
    candidate.max_likelihood_sequence = "max";
    candidate.left_flanking_sequence = "left";

    const auto actual { candidate.get_max_likelihood_sequence_with_flanks() };
    const std::string expected { "leftmax" };

    EXPECT_EQ(actual, expected);
}

TEST(GetMaxLikelihoodSequenceWithFlanksTest, allSequencesPresentReturnAllJoined)
{
    CandidateRegion candidate { Interval(0, 2), "test", 6 };
    candidate.max_likelihood_sequence = "max";
    candidate.left_flanking_sequence = "left";
    candidate.right_flanking_sequence = "right";

    const auto actual { candidate.get_max_likelihood_sequence_with_flanks() };
    const std::string expected { "leftmaxright" };

    EXPECT_EQ(actual, expected);
}

TEST(CandidateRegionEqualityTest, identicalCandidateRegionsReturnsTrue)
{
    const CandidateRegion candidate1 { Interval(0, 2), "test" };
    const CandidateRegion candidate2 { Interval(0, 2), "test" };

    EXPECT_TRUE(candidate1 == candidate2);
}

TEST(CandidateRegionEqualityTest, differentCandidateRegionsReturnsFalse)
{
    const CandidateRegion candidate1 { Interval(0, 2), "test" };
    const CandidateRegion candidate2 { Interval(0, 1), "test" };

    EXPECT_FALSE(candidate1 == candidate2);
}

TEST(CandidateRegionInequalityTest, identicalCandidateRegionsReturnsFalse)
{
    const CandidateRegion candidate1 { Interval(0, 2), "test" };
    const CandidateRegion candidate2 { Interval(0, 2), "test" };

    EXPECT_FALSE(candidate1 != candidate2);
}

TEST(CandidateRegionInequalityTest, differentCandidateRegionsReturnsTrue)
{
    const CandidateRegion candidate1 { Interval(0, 2), "test" };
    const CandidateRegion candidate2 { Interval(0, 2), "test2" };

    EXPECT_TRUE(candidate1 != candidate2);
}

TEST(IdentifyLowCoverageIntervalsTest, emptyCovgsReturnEmpty)
{
    const auto min_len { 5 };
    const auto min_covg { 0 };
    const std::vector<uint32_t> covgs;

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest, singleCovgPositionAboveThresholdReturnEmpty)
{
    const auto min_len { 1 };
    const auto min_covg { 1 };
    const std::vector<uint32_t> covgs { 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    singleCovgPositionBelowThresholdReturnSingleInterval)
{
    const auto min_len { 1 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(0, 1) };

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest, allPositionsAboveThresholdReturnEmpty)
{
    const auto min_len { 1 };
    const auto min_covg { 1 };
    const std::vector<uint32_t> covgs { 2, 2, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    allPositionsBelowThresholdReturnIntervalForWholeVector)
{
    const auto min_len { 1 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2, 2, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(0, 4) };

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    allPositionsBelowThresholdButLessThanMinLengthReturnEmpty)
{
    const auto min_len { 10 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2, 2, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    twoPositionsAtStartBelowThresholdButLessThanMinLenReturnEmpty)
{
    const auto min_len { 3 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2, 2, 4, 4, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    twoPositionsInMiddleBelowThresholdButLessThanMinLenReturnEmpty)
{
    const auto min_len { 3 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 2, 2, 4, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    twoPositionsAtEndBelowThresholdButLessThanMinLenReturnEmpty)
{
    const auto min_len { 3 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 4, 4, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    twoPositionsAtStartBelowThresholdReturnSingleInterval)
{
    const auto min_len { 2 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2, 2, 4, 4, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(0, 2) };

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    twoPositionsInMiddleBelowThresholdReturnSingleInterval)
{
    const auto min_len { 1 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 2, 2, 4, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(1, 3) };

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    twoPositionsAtEndBelowThresholdReturnSingleInterval)
{
    const auto min_len { 2 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 4, 4, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(3, 5) };

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    twoRegionsAtStartAndEndBelowThresholdReturnTwoIntervals)
{
    const auto min_len { 2 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2, 2, 4, 4, 4, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(0, 2), Interval(5, 7) };

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest, twoRegionsBelowThresholdReturnTwoIntervals)
{
    const auto min_len { 2 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 2, 1, 1, 4, 1, 2, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(1, 4), Interval(5, 7) };

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    twoRegionsBelowThresholdOneLessThanMinLengthReturnOneInterval)
{
    const auto min_len { 3 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 2, 1, 1, 4, 1, 2, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(1, 4) };

    EXPECT_EQ(actual, expected);
}
TEST(IdentifyLowCoverageIntervalsTest,
    twoRegionsBelowThresholdOneSameAsMaxLengthReturnTwoIntervals)
{
    const auto min_len { 2 };
    const auto max_len { 4 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 2, 1, 1, 1, 4, 1, 2, 4 };

    const auto actual { identify_low_coverage_intervals(
        covgs, min_covg, min_len, max_len) };
    const std::vector<Interval> expected { Interval(1, 5), Interval(6, 8) };

    EXPECT_EQ(actual, expected);
}

TEST(IdentifyLowCoverageIntervalsTest,
    twoRegionsBelowThresholdOneGreaterThanMaxLengthReturnOneInterval)
{
    const auto min_len { 2 };
    const auto max_len { 4 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 2, 1, 1, 1, 2, 4, 1, 2, 4 };

    const auto actual { identify_low_coverage_intervals(
        covgs, min_covg, min_len, max_len) };
    const std::vector<Interval> expected { Interval(7, 9) };

    EXPECT_EQ(actual, expected);
}
TEST(FindCandidateRegionsForPanNodeTest, emptyPanNodeReturnsNoCandidates)
{
    const auto num_samples { 1 };
    const auto prg_id { 3 };
    auto local_prg_ptr { std::make_shared<LocalPRG>(prg_id, "test", "") };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path;
    const std::vector<KmerNodePtr> kmer_node_max_likelihood_path;

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(
        local_prg_ptr, local_prg_ptr->id, num_samples) };

    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr,
        kmer_node_max_likelihood_path, local_node_max_likelihood_path };

    const CandidateRegions expected;
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components) };

    EXPECT_EQ(actual, expected);
}

TEST(FindCandidateRegionsForPanNodeTest, noCoverageReturnWholePrgAsCandidate)
{
    const auto num_samples { 1 };
    const auto prg_id { 3 };
    auto local_prg_ptr { std::make_shared<LocalPRG>(
        prg_id, "test", "AAA 5 G 6 C 5 TTT") };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg_ptr->minimizer_sketch(index, w, k);
    const std::string expected_sequence { "AAAGTTT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path {
        local_prg_ptr->prg.nodes[0], local_prg_ptr->prg.nodes[1],
        local_prg_ptr->prg.nodes[3]
    };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path {
        local_prg_ptr->kmernode_path_from_localnode_path(local_node_max_likelihood_path)
    };

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(
        local_prg_ptr, local_prg_ptr->id, num_samples) };

    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr,
        kmer_node_max_likelihood_path, local_node_max_likelihood_path };

    CandidateRegion expected_candidate { Interval(0, 7), pangraph_node->get_name() };
    expected_candidate.max_likelihood_sequence = expected_sequence;
    const CandidateRegions expected { std::make_pair(
        expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components) };
    const auto actual_sequence {
        actual.at(expected_candidate.get_id()).max_likelihood_sequence
    };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
}

TEST(FindCandidateRegionsForPanNodeTest, highCoverageReturnEmpty)
{
    const auto num_samples { 1 };
    const auto prg_id { 3 };
    auto local_prg_ptr { std::make_shared<LocalPRG>(
        prg_id, "test", "AAA 5 G 6 C 5 TTT") };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg_ptr->minimizer_sketch(index, w, k);
    const std::string expected_sequence { "AAAGTTT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path {
        local_prg_ptr->prg.nodes[0], local_prg_ptr->prg.nodes[1],
        local_prg_ptr->prg.nodes[3]
    };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path {
        local_prg_ptr->kmernode_path_from_localnode_path(local_node_max_likelihood_path)
    };

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(
        local_prg_ptr, local_prg_ptr->id, num_samples) };

    for (const auto& kmer_node :
        pangraph_node->kmer_prg_with_coverage.kmer_prg->nodes) {
        pangraph_node->kmer_prg_with_coverage.set_covg(kmer_node->id, 100, false, 0);
    }

    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr,
        kmer_node_max_likelihood_path, local_node_max_likelihood_path };

    const CandidateRegions expected;
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components) };

    EXPECT_EQ(actual, expected);
}

TEST(
    FindCandidateRegionsForPanNodeTest, noCoverageOnFiveBasesReturnFiveBasesAsCandidate)
{
    const auto num_samples { 1 };
    const auto prg_id { 3 };
    auto local_prg_ptr { std::make_shared<LocalPRG>(
        prg_id, "test", "AAAA 5 GGG 6 CCC 5 TTTT") };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg_ptr->minimizer_sketch(index, w, k);
    const std::string expected_sequence { "AAAAGGGTTTT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path {
        local_prg_ptr->prg.nodes[0], local_prg_ptr->prg.nodes[1],
        local_prg_ptr->prg.nodes[3]
    };
    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 3, 5, 7, 9, 11, 13,
        14, 15 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());
    for (const auto& idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg_ptr->kmer_prg.nodes[idx]);
    }

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(
        local_prg_ptr, local_prg_ptr->id, num_samples) };
    const std::vector<int> kmer_node_idxs_for_high_covg { 0, 1, 14, 15 };
    for (const auto& idx : kmer_node_idxs_for_high_covg) {
        pangraph_node->kmer_prg_with_coverage.set_covg(idx, 100, false, 0);
    }

    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr,
        kmer_node_max_likelihood_path, local_node_max_likelihood_path };

    CandidateRegion expected_candidate { Interval(3, 8), pangraph_node->get_name() };
    const CandidateRegions expected { std::make_pair(
        expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components) };
    const auto actual_sequence {
        actual.at(expected_candidate.get_id()).get_max_likelihood_sequence_with_flanks()
    };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
}

TEST(FindCandidateRegionsForPanNodeTest,
    noCoverageOnFiveBasesReturnFiveBasesPlusPaddingAsCandidate)
{
    const auto prg_id { 3 };
    auto local_prg_ptr { std::make_shared<LocalPRG>(
        prg_id, "test", "AAAA 5 GGG 6 CCC 5 TTTT") };
    const std::string expected_sequence { "AAAAGGGTTTT" };
    const auto expected_max_likelihood_sequence { "AAGGGTT" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg_ptr->minimizer_sketch(index, w, k);

    const std::vector<int> local_node_idxs_for_max_path { 0, 1, 3 };
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    local_node_max_likelihood_path.reserve(local_node_idxs_for_max_path.size());

    for (const auto& idx : local_node_idxs_for_max_path) {
        local_node_max_likelihood_path.push_back(local_prg_ptr->prg.nodes[idx]);
    }

    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 3, 5, 7, 9, 11, 13,
        14, 15 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());

    for (const auto& idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg_ptr->kmer_prg.nodes[idx]);
    }

    const auto num_samples { 1 };
    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(
        local_prg_ptr, local_prg_ptr->id, num_samples) };
    const std::vector<int> kmer_node_idxs_for_high_covg { 0, 1, 14, 15 };

    for (const auto& idx : kmer_node_idxs_for_high_covg) {
        pangraph_node->kmer_prg_with_coverage.set_covg(idx, 100, false, 0);
    }

    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr,
        kmer_node_max_likelihood_path, local_node_max_likelihood_path };

    const auto interval_padding { 1 };
    CandidateRegion expected_candidate { Interval(3, 8), pangraph_node->get_name(),
        interval_padding };
    const CandidateRegions expected { std::make_pair(
        expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(
        pangraph_node_components, interval_padding) };
    const auto actual_sequence {
        actual.at(expected_candidate.get_id()).get_max_likelihood_sequence_with_flanks()
    };
    const auto actual_max_likelihood_sequence {
        actual.at(expected_candidate.get_id()).max_likelihood_sequence
    };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
    EXPECT_EQ(actual_max_likelihood_sequence, expected_max_likelihood_sequence);
}

TEST(FindCandidateRegionsForPanNodeTest,
    noCoverageOnStartFiveBasesReturnFiveBasesAsCandidate)
{
    const auto prg_id { 3 };
    auto local_prg_ptr { std::make_shared<LocalPRG>(
        prg_id, "test", "AAAA 5 GGG 6 CCC 5 TTTT") };
    const std::string expected_sequence { "AAAAGGGTTTT" };
    const auto expected_max_likelihood_sequence { "AAAAG" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg_ptr->minimizer_sketch(index, w, k);

    const std::vector<int> local_node_idxs_for_max_path { 0, 1, 3 };
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    local_node_max_likelihood_path.reserve(local_node_idxs_for_max_path.size());

    for (const auto& idx : local_node_idxs_for_max_path) {
        local_node_max_likelihood_path.push_back(local_prg_ptr->prg.nodes[idx]);
    }

    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 3, 5, 7, 9, 11, 13,
        14, 15 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());

    for (const auto& idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg_ptr->kmer_prg.nodes[idx]);
    }

    const auto num_samples { 1 };
    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(
        local_prg_ptr, local_prg_ptr->id, num_samples) };
    const std::vector<int> kmer_node_idxs_for_high_covg { 9, 11, 13, 14, 15 };

    for (const auto& idx : kmer_node_idxs_for_high_covg) {
        pangraph_node->kmer_prg_with_coverage.set_covg(idx, 100, false, 0);
    }

    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr,
        kmer_node_max_likelihood_path, local_node_max_likelihood_path };

    const auto interval_padding { 0 };
    CandidateRegion expected_candidate { Interval(0, 5), pangraph_node->get_name(),
        interval_padding };
    const CandidateRegions expected { std::make_pair(
        expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(
        pangraph_node_components, interval_padding) };
    const auto actual_sequence {
        actual.at(expected_candidate.get_id()).get_max_likelihood_sequence_with_flanks()
    };
    const auto actual_max_likelihood_sequence {
        actual.at(expected_candidate.get_id()).max_likelihood_sequence
    };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
    EXPECT_EQ(actual_max_likelihood_sequence, expected_max_likelihood_sequence);
}

TEST(FindCandidateRegionsForPanNodeTest,
    noCoverageOnEndFiveBasesReturnFiveBasesAsCandidate)
{
    const auto prg_id { 3 };
    auto local_prg_ptr { std::make_shared<LocalPRG>(
        prg_id, "test", "AAAA 5 GGG 6 CCC 5 TTTT") };
    const std::string expected_sequence { "AAAAGGGTTTT" };
    const auto expected_max_likelihood_sequence { "GTTTT" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg_ptr->minimizer_sketch(index, w, k);

    const std::vector<int> local_node_idxs_for_max_path { 0, 1, 3 };
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    local_node_max_likelihood_path.reserve(local_node_idxs_for_max_path.size());

    for (const auto& idx : local_node_idxs_for_max_path) {
        local_node_max_likelihood_path.push_back(local_prg_ptr->prg.nodes[idx]);
    }

    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 3, 5, 7, 9, 11, 13,
        14, 15 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());

    for (const auto& idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg_ptr->kmer_prg.nodes[idx]);
    }

    const auto num_samples { 1 };
    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(
        local_prg_ptr, local_prg_ptr->id, num_samples) };
    const std::vector<int> kmer_node_idxs_for_high_covg { 0, 1, 2, 3, 5 };

    for (const auto& idx : kmer_node_idxs_for_high_covg) {
        pangraph_node->kmer_prg_with_coverage.set_covg(idx, 100, false, 0);
    }

    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr,
        kmer_node_max_likelihood_path, local_node_max_likelihood_path };

    const auto interval_padding { 0 };
    CandidateRegion expected_candidate { Interval(6, 11), pangraph_node->get_name(),
        interval_padding };
    const CandidateRegions expected { std::make_pair(
        expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(
        pangraph_node_components, interval_padding) };
    const auto actual_sequence {
        actual.at(expected_candidate.get_id()).get_max_likelihood_sequence_with_flanks()
    };
    const auto actual_max_likelihood_sequence {
        actual.at(expected_candidate.get_id()).max_likelihood_sequence
    };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
    EXPECT_EQ(actual_max_likelihood_sequence, expected_max_likelihood_sequence);
}

TEST(FindCandidateRegionsForPanNodeTest,
    noCoverageOnFiveBasesWithinDoubleNestingReturnFiveBasesPlusPaddingAsCandidate)
{
    const auto prg_id { 3 };
    auto local_prg_ptr { std::make_shared<LocalPRG>(
        prg_id, "test", "AAAA 5 CCCC 6 GG 7 XXX 8 YYY 7 GG 5 TTTT") };
    const std::string expected_sequence { "AAAAGGYYYGGTTTT" };
    const auto expected_max_likelihood_sequence { "GGYYYGG" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg_ptr->minimizer_sketch(index, w, k);

    const std::vector<int> local_node_idxs_for_max_path { 0, 2, 4, 5, 6 };
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    local_node_max_likelihood_path.reserve(local_node_idxs_for_max_path.size());

    for (const auto& idx : local_node_idxs_for_max_path) {
        local_node_max_likelihood_path.push_back(local_prg_ptr->prg.nodes[idx]);
    }

    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 4, 6, 9, 12, 15, 18,
        21, 23, 24, 19, 22, 25 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());

    for (const auto& idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg_ptr->kmer_prg.nodes[idx]);
    }

    const auto num_samples { 1 };
    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(
        local_prg_ptr, local_prg_ptr->id, num_samples) };
    const std::vector<int> kmer_node_idxs_for_high_covg { 0, 1, 2, 4, 24, 19, 22, 25 };

    for (const auto& idx : kmer_node_idxs_for_high_covg) {
        pangraph_node->kmer_prg_with_coverage.set_covg(idx, 100, false, 0);
    }

    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr,
        kmer_node_max_likelihood_path, local_node_max_likelihood_path };

    const auto interval_padding { 1 };
    CandidateRegion expected_candidate { Interval(5, 10), pangraph_node->get_name(),
        interval_padding };
    const CandidateRegions expected { std::make_pair(
        expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(
        pangraph_node_components, interval_padding) };
    const auto actual_sequence {
        actual.at(expected_candidate.get_id()).get_max_likelihood_sequence_with_flanks()
    };
    const auto actual_max_likelihood_sequence {
        actual.at(expected_candidate.get_id()).max_likelihood_sequence
    };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
    EXPECT_EQ(actual_max_likelihood_sequence, expected_max_likelihood_sequence);
}

TEST(FindCandidateRegionsForPanNodeTest,
    noCoverageOnTwoFiveBaseRegionsWithinDoubleNestingReturnTwoRegionsAsCandidate)
{
    const auto prg_id { 3 };
    auto local_prg_ptr { std::make_shared<LocalPRG>(
        prg_id, "test", "AAAA 5 CCCC 6 GG 7 XXX 8 YYY 7 GG 5 TTTTT") };
    const std::string expected_sequence { "AAAAGGYYYGGTTTTT" };
    std::vector<std::string> expected_max_likelihood_sequences { "AAAAG", "YGGTT" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg_ptr->minimizer_sketch(index, w, k);

    const std::vector<int> local_node_idxs_for_max_path { 0, 2, 4, 5, 6 };
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    local_node_max_likelihood_path.reserve(local_node_idxs_for_max_path.size());

    for (const auto& idx : local_node_idxs_for_max_path) {
        local_node_max_likelihood_path.push_back(local_prg_ptr->prg.nodes[idx]);
    }

    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 4, 6, 9, 12, 15, 18,
        21, 23, 25, 19, 22, 24, 26 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());

    for (const auto& idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg_ptr->kmer_prg.nodes[idx]);
    }

    const auto num_samples { 1 };
    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(
        local_prg_ptr, local_prg_ptr->id, num_samples) };
    const std::vector<int> kmer_node_idxs_for_high_covg { 12, 24, 26 };

    for (const auto& idx : kmer_node_idxs_for_high_covg) {
        pangraph_node->kmer_prg_with_coverage.set_covg(idx, 100, false, 0);
    }

    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr,
        kmer_node_max_likelihood_path, local_node_max_likelihood_path };

    const auto interval_padding { 0 };
    CandidateRegion expected_candidate1 { Interval(0, 5), pangraph_node->get_name(),
        interval_padding };
    CandidateRegion expected_candidate2 { Interval(8, 13), pangraph_node->get_name(),
        interval_padding };
    const CandidateRegions expected { std::make_pair(expected_candidate1.get_id(),
                                          expected_candidate1),
        std::make_pair(expected_candidate2.get_id(), expected_candidate2) };
    const auto actual { find_candidate_regions_for_pan_node(
        pangraph_node_components, interval_padding) };

    EXPECT_EQ(actual, expected);
    std::vector<std::string> actual_max_likelihood_sequences;
    actual_max_likelihood_sequences.reserve(2);
    for (const auto& actual_candidate : actual) {
        EXPECT_EQ(actual_candidate.second.get_max_likelihood_sequence_with_flanks(),
            expected_sequence);
        actual_max_likelihood_sequences.push_back(
            actual_candidate.second.max_likelihood_sequence);
    }
    std::sort(
        actual_max_likelihood_sequences.begin(), actual_max_likelihood_sequences.end());
    std::sort(expected_max_likelihood_sequences.begin(),
        expected_max_likelihood_sequences.end());
    EXPECT_EQ(actual_max_likelihood_sequences, expected_max_likelihood_sequences);
}

TEST(AddPileupEntryForCandidateRegionTest, emptyReadWithReadCoordsPileupsEmpty)
{
    const auto read_sequence { "" };
    CandidateRegion candidate { Interval(1, 3), "test" };
    ReadCoordinate read_coord_1 { 0, 6, 10, true };

    candidate.read_coordinates.insert(read_coord_1);
    candidate.add_pileup_entry(read_sequence, read_coord_1);
    const auto& actual { candidate.pileup };
    const ReadPileup expected;

    EXPECT_EQ(actual, expected);
}

TEST(AddPileupEntryForCandidateRegionTest, oneReadOneReadCoordPileupHasOneEntry)
{
    const std::string read_sequence { "XXXFOOXXX" };
    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord { 0, 3, 6, true }; // FOO
    candidate.read_coordinates.insert(read_coord);

    candidate.add_pileup_entry(read_sequence, read_coord);

    const auto& actual { candidate.pileup };
    const ReadPileup expected { "FOO" };

    EXPECT_EQ(actual, expected);
}

TEST(AddPileupEntryForCandidateRegionTest,
    oneReadOneWholeReadCoordPileupHasOneEntrySameAsRead)
{
    const std::string read_sequence { "XXXFOOXXX" };
    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord { 0, 0, 10, true }; // XXXFOOXXX
    candidate.read_coordinates.insert(read_coord);

    candidate.add_pileup_entry(read_sequence, read_coord);

    const auto& actual { candidate.pileup };
    const ReadPileup expected { "XXXFOOXXX" };

    EXPECT_EQ(actual, expected);
}

TEST(AddPileupEntryForCandidateRegionTest,
    oneReadOneReadCoordThatRunsPastEndOfReadPileupHasOneEntryUpToEndOfRead)
{
    const std::string read_sequence { "XXXFOOXXX" };

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord { 0, 5, 20, true }; // OXXX
    candidate.read_coordinates.insert(read_coord);

    candidate.add_pileup_entry(read_sequence, read_coord);

    const auto& actual { candidate.pileup };
    const ReadPileup expected { "OXXX" };

    EXPECT_EQ(actual, expected);
}

TEST(AddPileupEntryForCandidateRegionTest,
    oneReadOneReverseReadCoordThatRunsPastEndOfReadPileupHasOneEntryUpToEndOfRead)
{
    const std::string read_sequence { "AATTCCGG" };

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord { 0, 5, 20, false }; // CCG
    candidate.read_coordinates.insert(read_coord);

    candidate.add_pileup_entry(read_sequence, read_coord);

    const auto& actual { candidate.pileup };
    const ReadPileup expected { "CCG" };

    EXPECT_EQ(actual, expected);
}

TEST(AddPileupEntryForCandidateRegionTest, oneReadOneReadCoordOutsideReadPileupEmpty)
{

    const std::string read_sequence { "XXXFOOXXX" };

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord { 0, 15, 20, true };
    candidate.read_coordinates.insert(read_coord);

    candidate.add_pileup_entry(read_sequence, read_coord);

    const auto& actual { candidate.pileup };
    const ReadPileup expected;

    EXPECT_EQ(actual, expected);
}

TEST(LoadAllCandidateRegionsPileupsFromFastq, emptyCandidateRegionReturnsEmptyPileup)
{
    Fastaq temp_fastq { false, true };
    auto read_name { "0" };
    std::string read_sequence { "AATTCCGG" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    read_name = "1";
    read_sequence = "GATTACAA";
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegions candidate_regions;
    auto pileup_construction_map = construct_pileup_construction_map(candidate_regions);
    load_all_candidate_regions_pileups_from_fastq(
        temp_reads_filepath, candidate_regions, pileup_construction_map);

    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    EXPECT_TRUE(candidate_regions.empty());
}

TEST(AddPileupEntryForCandidateRegionTest,
    oneCandidateZeroReadsInFilePileupHasZeroEntries)
{

    Fastaq temp_fastq { false, true };
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord1 { 0, 2, 4, true };
    ReadCoordinate read_coord2 { 1, 3, 6, false };
    candidate.read_coordinates.insert(read_coord1);
    candidate.read_coordinates.insert(read_coord2);

    CandidateRegions candidate_regions { std::make_pair(
        candidate.get_id(), candidate) };
    auto pileup_construction_map = construct_pileup_construction_map(candidate_regions);
    load_all_candidate_regions_pileups_from_fastq(
        temp_reads_filepath, candidate_regions, pileup_construction_map);

    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto& actual { candidate_regions.at(candidate.get_id()).pileup };
    const ReadPileup expected {};

    EXPECT_EQ(actual, expected);
}

TEST(AddPileupEntryForCandidateRegionTest,
    oneCandidateTwoReadsInFileOneReadInCandidatePileupHasOneEntry)
{

    Fastaq temp_fastq { false, true };
    auto read_name { "0" };
    std::string read_sequence { "AATTCCGG" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    read_name = "5";
    read_sequence = "GATTACAA";
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord1 { 0, 2, 4, true }; // TT
    ReadCoordinate read_coord2 { 4, 3, 6, false }; // no corresponding read
    candidate.read_coordinates.insert(read_coord1);
    candidate.read_coordinates.insert(read_coord2);

    CandidateRegions candidate_regions { std::make_pair(
        candidate.get_id(), candidate) };
    auto pileup_construction_map = construct_pileup_construction_map(candidate_regions);
    load_all_candidate_regions_pileups_from_fastq(
        temp_reads_filepath, candidate_regions, pileup_construction_map);

    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto& actual { candidate_regions.at(candidate.get_id()).pileup };
    const ReadPileup expected { "TT" };

    EXPECT_EQ(actual, expected);
}

TEST(AddPileupEntryForCandidateRegionTest,
    oneCandidateTwoReadsInFileTwoReadsInCandidatePileupHasTwoEntries)
{

    Fastaq temp_fastq { false, true };
    auto read_name { "0" };
    std::string read_sequence { "AATTCCGG" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    read_name = "5";
    read_sequence = "GATTACAA";
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord1 { 0, 2, 4, true }; // TT
    ReadCoordinate read_coord2 { 1, 3, 6, false }; // GTA
    candidate.read_coordinates.insert(read_coord1);
    candidate.read_coordinates.insert(read_coord2);

    CandidateRegions candidate_regions { std::make_pair(
        candidate.get_id(), candidate) };
    auto pileup_construction_map = construct_pileup_construction_map(candidate_regions);
    load_all_candidate_regions_pileups_from_fastq(
        temp_reads_filepath, candidate_regions, pileup_construction_map);

    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto& actual { candidate_regions.at(candidate.get_id()).pileup };
    const ReadPileup expected { "TT", "GTA" };

    EXPECT_EQ(actual, expected);
}

TEST(AddPileupEntryForCandidateRegionTest,
    twoCandidatesZeroReadsInFilePileupHasZeroEntries)
{

    Fastaq temp_fastq { false, true };
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate_1 { Interval(0, 3), "test" };
    ReadCoordinate read_coord1 { 0, 2, 4, true };
    ReadCoordinate read_coord2 { 1, 3, 6, false };
    candidate_1.read_coordinates.insert(read_coord1);
    candidate_1.read_coordinates.insert(read_coord2);

    CandidateRegion candidate_2 { Interval(0, 3), "test2" };
    ReadCoordinate read_coord3 { 0, 2, 4, true };
    ReadCoordinate read_coord4 { 1, 3, 6, false };
    candidate_2.read_coordinates.insert(read_coord3);
    candidate_2.read_coordinates.insert(read_coord4);

    CandidateRegions candidate_regions { std::make_pair(
                                             candidate_1.get_id(), candidate_1),
        std::make_pair(candidate_2.get_id(), candidate_2) };
    auto pileup_construction_map = construct_pileup_construction_map(candidate_regions);
    load_all_candidate_regions_pileups_from_fastq(
        temp_reads_filepath, candidate_regions, pileup_construction_map);

    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto& actual_1 { candidate_regions.at(candidate_1.get_id()).pileup };
    const auto& actual_2 { candidate_regions.at(candidate_2.get_id()).pileup };
    const ReadPileup expected {};

    EXPECT_EQ(actual_1, expected);
    EXPECT_EQ(actual_2, expected);
}

TEST(AddPileupEntryForCandidateRegionTest,
    twoCandidatesTwoReadsInFileOneReadInOneCandidatePileupHasOneEntryInOneCandidate)
{

    Fastaq temp_fastq { false, true };
    auto read_name { "0" };
    std::string read_sequence { "AATTCCGG" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    read_name = "5";
    read_sequence = "GATTACAA";
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate_1 { Interval(0, 3), "test" };
    ReadCoordinate read_coord1 { 0, 2, 4, true };
    candidate_1.read_coordinates.insert(read_coord1);

    CandidateRegion candidate_2 { Interval(0, 3), "test2" };
    ReadCoordinate read_coord2 { 4, 3, 6, false };
    candidate_2.read_coordinates.insert(read_coord2);

    CandidateRegions candidate_regions { std::make_pair(
                                             candidate_1.get_id(), candidate_1),
        std::make_pair(candidate_2.get_id(), candidate_2) };
    auto pileup_construction_map = construct_pileup_construction_map(candidate_regions);
    load_all_candidate_regions_pileups_from_fastq(
        temp_reads_filepath, candidate_regions, pileup_construction_map);

    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto& actual_1 { candidate_regions.at(candidate_1.get_id()).pileup };
    const auto& actual_2 { candidate_regions.at(candidate_2.get_id()).pileup };
    const ReadPileup expected_1 { "TT" };
    const ReadPileup expected_2 {};

    EXPECT_EQ(actual_1, expected_1);
    EXPECT_EQ(actual_2, expected_2);
}

TEST(AddPileupEntryForCandidateRegionTest,
    twoCandidatesTwoReadsInFileOneReadInEachCandidatePileupHasOneEntryInEachCandidate)
{

    Fastaq temp_fastq { false, true };
    auto read_name { "0" };
    std::string read_sequence { "AATTCCGG" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    read_name = "1";
    read_sequence = "GATTACAA";
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate_1 { Interval(0, 3), "test" };
    ReadCoordinate read_coord1 { 0, 2, 4, true }; // TT
    candidate_1.read_coordinates.insert(read_coord1);

    CandidateRegion candidate_2 { Interval(0, 3), "test2" };
    ReadCoordinate read_coord2 { 1, 3, 6, true }; // TAC
    candidate_2.read_coordinates.insert(read_coord2);

    CandidateRegions candidate_regions { std::make_pair(
                                             candidate_1.get_id(), candidate_1),
        std::make_pair(candidate_2.get_id(), candidate_2) };
    auto pileup_construction_map = construct_pileup_construction_map(candidate_regions);
    load_all_candidate_regions_pileups_from_fastq(
        temp_reads_filepath, candidate_regions, pileup_construction_map);

    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto& actual_1 { candidate_regions.at(candidate_1.get_id()).pileup };
    const auto& actual_2 { candidate_regions.at(candidate_2.get_id()).pileup };
    const ReadPileup expected_1 { "TT" };
    const ReadPileup expected_2 { "TAC" };

    EXPECT_EQ(actual_1, expected_1);
    EXPECT_EQ(actual_2, expected_2);
}

TEST(AddPileupEntryForCandidateRegionTest,
    twoCandidatesTwoReadsInFileOneReadInEachCandidatePileupHasOneEntryInEachCandidateUsingFourThreads)
{

    const uint32_t threads = 4;
    Fastaq temp_fastq { false, true };
    auto read_name { "0" };
    std::string read_sequence { "AATTCCGG" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    for (uint32_t i = 0; i < threads / 2 * 1000; ++i) {
        temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    }
    read_name = "1";
    read_sequence = "GATTACAA";
    for (uint32_t i = 0; i < threads / 2 * 1000; ++i) {
        temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    }
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate_1 { Interval(0, 3), "test" };
    ReadCoordinate read_coord1 { 0, 2, 4, true }; // TT
    candidate_1.read_coordinates.insert(read_coord1);

    CandidateRegion candidate_2 { Interval(0, 3), "test2" };
    ReadCoordinate read_coord2 { 3000, 3, 6, true }; // TAC
    candidate_2.read_coordinates.insert(read_coord2);

    CandidateRegions candidate_regions { std::make_pair(
                                             candidate_1.get_id(), candidate_1),
        std::make_pair(candidate_2.get_id(), candidate_2) };
    auto pileup_construction_map = construct_pileup_construction_map(candidate_regions);
    load_all_candidate_regions_pileups_from_fastq(
        temp_reads_filepath, candidate_regions, pileup_construction_map, threads);

    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto& actual_1 { candidate_regions.at(candidate_1.get_id()).pileup };
    const auto& actual_2 { candidate_regions.at(candidate_2.get_id()).pileup };
    const ReadPileup expected_1 { "TT" };
    const ReadPileup expected_2 { "TAC" };

    EXPECT_EQ(actual_1, expected_1);
    EXPECT_EQ(actual_2, expected_2);
}

std::string read_file_to_string(const fs::path& filepath)
{
    fs::ifstream f(filepath);
    std::stringstream buffer;
    buffer << f.rdbuf();
    return buffer.str();
}

TEST(WriteDenovoPathsToFileTest, noReadsDoesntWriteFile)
{
    const std::string name { "test" };
    CandidateRegion candidate { Interval(0, 1), name };
    const fs::path output_directory { "/tmp" };
    const fs::path filepath { output_directory / candidate.filename };

    candidate.write_denovo_paths_to_file(output_directory);

    EXPECT_FALSE(fs::exists(filepath));
}

TEST(WriteDenovoPathsToFileTest, twoReadsWritesTwoReadsToFile)
{
    const std::string name { "test" };
    CandidateRegion candidate { Interval(0, 1), name };
    const fs::path output_directory { "/tmp" };
    const fs::path filepath { output_directory / candidate.filename };
    const DenovoPaths paths { "shrubberies", "ni" };
    candidate.denovo_paths = paths;

    candidate.write_denovo_paths_to_file(output_directory);

    std::stringstream expected_ss;

    for (size_t i = 0; i < paths.size(); ++i) {
        expected_ss << ">test." << std::to_string(i) << std::endl
                    << paths.at(i) << std::endl;
    }

    const std::string expected { expected_ss.str() };
    const auto actual { read_file_to_string(filepath) };

    fs::remove(filepath);

    EXPECT_EQ(actual, expected);
}

TEST(ConstructPileupConstructionMapTest, emptyInEmptyOut)
{
    CandidateRegions empty_candidate_regions;
    const auto actual { construct_pileup_construction_map(empty_candidate_regions) };
    PileupConstructionMap expected;
    EXPECT_EQ(actual, expected);
}

void compare_maps(const PileupConstructionMap& map1, const PileupConstructionMap& map2)
{
    std::vector<ReadId> keys1, keys2;
    for (const auto& item : map1)
        keys1.push_back(item.first);
    for (const auto& item : map2)
        keys2.push_back(item.first);

    EXPECT_EQ(keys1, keys2);

    for (auto key : keys1) {
        auto vector1 = map1.at(key);
        auto vector2 = map2.at(key);

        // TODO: find the correct way to compare a vector of pair of pointers with
        // google test
        // TODO: or even better: compare the maps themselves
        EXPECT_TRUE(vector1.size() == vector2.size());
        EXPECT_TRUE(std::is_permutation(vector1.begin(), vector1.end(), vector2.begin(),
            [](const std::pair<CandidateRegion*, const ReadCoordinate*>& pair1,
                const std::pair<CandidateRegion*, const ReadCoordinate*>& pair2) {
                return *pair1.first == *pair2.first && *pair1.second == *pair2.second;
            }));
    }
}

TEST(ConstructPileupConstructionMapTest, oneCandidateRegionOneReadCoordinate)
{
    CandidateRegion candidate { Interval(0, 3), "test" };
    const uint32_t read_id = 0;
    ReadCoordinate read_coord { read_id, 5, 20, true };
    candidate.read_coordinates.insert(read_coord);
    CandidateRegions candidate_regions { std::make_pair(
        candidate.get_id(), candidate) };

    const auto actual { construct_pileup_construction_map(candidate_regions) };
    PileupConstructionMap expected;
    expected[read_id].emplace_back(&candidate, &read_coord);
    compare_maps(actual, expected);
}

TEST(ConstructPileupConstructionMapTest, oneCandidateRegionTwoReadCoordinatesSameReadId)
{
    CandidateRegion candidate { Interval(0, 3), "test" };
    const uint32_t read_id = 0;
    ReadCoordinate read_coord_1 { read_id, 5, 20, true };
    candidate.read_coordinates.insert(read_coord_1);
    ReadCoordinate read_coord_2 { read_id, 50, 100, true };
    candidate.read_coordinates.insert(read_coord_2);
    CandidateRegions candidate_regions { std::make_pair(
        candidate.get_id(), candidate) };

    const auto actual { construct_pileup_construction_map(candidate_regions) };
    PileupConstructionMap expected;
    expected[read_id].emplace_back(&candidate, &read_coord_1);
    expected[read_id].emplace_back(&candidate, &read_coord_2);

    compare_maps(actual, expected);
}

TEST(ConstructPileupConstructionMapTest,
    oneCandidateRegionTwoReadCoordinatesDifferentReadId)
{
    CandidateRegion candidate { Interval(0, 3), "test" };
    const uint32_t read_id_1 = 0;
    ReadCoordinate read_coord_1 { read_id_1, 5, 20, true };
    candidate.read_coordinates.insert(read_coord_1);
    const uint32_t read_id_2 = 2;
    ReadCoordinate read_coord_2 { read_id_2, 50, 100, true };
    candidate.read_coordinates.insert(read_coord_2);
    CandidateRegions candidate_regions { std::make_pair(
        candidate.get_id(), candidate) };

    const auto actual { construct_pileup_construction_map(candidate_regions) };
    PileupConstructionMap expected;
    expected[read_id_1].emplace_back(&candidate, &read_coord_1);
    expected[read_id_2].emplace_back(&candidate, &read_coord_2);

    compare_maps(actual, expected);
}

TEST(ConstructPileupConstructionMapTest,
    twoCandidateRegionsOneReadCoordinateEachDifferentReadId)
{
    CandidateRegion candidate_1 { Interval(0, 3), "test" };
    const uint32_t read_id_1 = 0;
    ReadCoordinate read_coord_1 { read_id_1, 5, 20, true };
    candidate_1.read_coordinates.insert(read_coord_1);
    CandidateRegion candidate_2 { Interval(10, 20), "test2" };
    const uint32_t read_id_2 = 2;
    ReadCoordinate read_coord_2 { read_id_2, 50, 100, true };
    candidate_2.read_coordinates.insert(read_coord_2);
    CandidateRegions candidate_regions { std::make_pair(
                                             candidate_1.get_id(), candidate_1),
        std::make_pair(candidate_2.get_id(), candidate_2) };

    const auto actual { construct_pileup_construction_map(candidate_regions) };
    PileupConstructionMap expected;
    expected[read_id_1].emplace_back(&candidate_1, &read_coord_1);
    expected[read_id_2].emplace_back(&candidate_2, &read_coord_2);

    compare_maps(actual, expected);
}

TEST(ConstructPileupConstructionMapTest,
    twoCandidateRegionsOneReadCoordinateEachSameReadId)
{
    CandidateRegion candidate_1 { Interval(0, 3), "test" };
    const uint32_t read_id = 0;
    ReadCoordinate read_coord_1 { read_id, 5, 20, true };
    candidate_1.read_coordinates.insert(read_coord_1);
    CandidateRegion candidate_2 { Interval(10, 20), "test2" };
    ReadCoordinate read_coord_2 { read_id, 50, 100, true };
    candidate_2.read_coordinates.insert(read_coord_2);
    CandidateRegions candidate_regions { std::make_pair(
                                             candidate_1.get_id(), candidate_1),
        std::make_pair(candidate_2.get_id(), candidate_2) };

    const auto actual { construct_pileup_construction_map(candidate_regions) };
    PileupConstructionMap expected;
    expected[read_id].emplace_back(&candidate_1, &read_coord_1);
    expected[read_id].emplace_back(&candidate_2, &read_coord_2);

    compare_maps(actual, expected);
}