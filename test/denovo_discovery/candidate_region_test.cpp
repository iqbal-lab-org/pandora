#include <iostream>
#include <cstdio>
#include "gtest/gtest.h"
#include "../test_macro.cpp"
#include "denovo_discovery/candidate_region.h"
#include "fastaq.h"


TEST(CandidateRegionGetIntervalTest, noPaddingReturnsOriginalInterval) {
    const CandidateRegion candidate { Interval(0, 2), "test" };

    const auto actual { candidate.get_interval() };
    const Interval expected { 0, 2 };

    EXPECT_EQ(actual, expected);
}


TEST(CandidateRegionGetIntervalTest, withPaddingLessThanIntervalStartReturnsOriginalIntervalWithPadding) {
    const CandidateRegion candidate { Interval(3, 4), "test", 1 };

    const auto actual { candidate.get_interval() };
    const Interval expected { 2, 5 };

    EXPECT_EQ(actual, expected);
}


TEST(CandidateRegionGetIntervalTest, withPaddingGreaterThanIntervalStartReturnsOriginalIntervalWithPadding) {
    const CandidateRegion candidate { Interval(3, 4), "test", 5 };

    const auto actual { candidate.get_interval() };
    const Interval expected { 0, 9 };

    EXPECT_EQ(actual, expected);
}


TEST(CandidateRegionGetNameTest, testCorrectValueRetrieved) {
    const CandidateRegion candidate { Interval(0, 2), "test" };

    const auto &actual { candidate.get_name() };
    const std::string expected { "test" };

    EXPECT_EQ(actual, expected);
}


TEST(CandidateRegionGetIdTest, noPaddingTestCorrectIdRetrieved) {
    const CandidateRegion candidate { Interval(0, 2), "test" };

    const auto actual { candidate.get_id() };
    const CandidateRegionIdentifier expected { candidate.get_interval(), candidate.get_name() };

    EXPECT_EQ(actual, expected);
}


TEST(CandidateRegionGetIdTest, withPaddingTestCorrectIdRetrieved) {
    const CandidateRegion candidate { Interval(0, 2), "test", 6 };

    const auto actual { candidate.get_id() };
    const CandidateRegionIdentifier expected { candidate.get_interval(), candidate.get_name() };

    EXPECT_EQ(actual, expected);
}


TEST(GetMaxLikelihoodSequenceWithFlanksTest, noSequencesReturnEmptyString) {
    const CandidateRegion candidate { Interval(0, 2), "test", 6 };

    const auto actual { candidate.get_max_likelihood_sequence_with_flanks() };
    const std::string expected;

    EXPECT_EQ(actual, expected);
}


TEST(GetMaxLikelihoodSequenceWithFlanksTest, noLeftFlankSequencesReturnMaxAndRightSequence) {
    CandidateRegion candidate { Interval(0, 2), "test", 6 };
    candidate.max_likelihood_sequence = "max";
    candidate.right_flanking_sequence = "right";

    const auto actual { candidate.get_max_likelihood_sequence_with_flanks() };
    const std::string expected { "maxright" };

    EXPECT_EQ(actual, expected);
}


TEST(GetMaxLikelihoodSequenceWithFlanksTest, noRightFlankSequencesReturnLeftAndMaxSequence) {
    CandidateRegion candidate { Interval(0, 2), "test", 6 };
    candidate.max_likelihood_sequence = "max";
    candidate.left_flanking_sequence = "left";

    const auto actual { candidate.get_max_likelihood_sequence_with_flanks() };
    const std::string expected { "leftmax" };

    EXPECT_EQ(actual, expected);
}


TEST(GetMaxLikelihoodSequenceWithFlanksTest, allSequencesPresentReturnAllJoined) {
    CandidateRegion candidate { Interval(0, 2), "test", 6 };
    candidate.max_likelihood_sequence = "max";
    candidate.left_flanking_sequence = "left";
    candidate.right_flanking_sequence = "right";

    const auto actual { candidate.get_max_likelihood_sequence_with_flanks() };
    const std::string expected { "leftmaxright" };

    EXPECT_EQ(actual, expected);
}


TEST(CandidateRegionEqualityTest, identicalCandidateRegionsReturnsTrue) {
    const CandidateRegion candidate1 { Interval(0, 2), "test" };
    const CandidateRegion candidate2 { Interval(0, 2), "test" };

    EXPECT_TRUE(candidate1 == candidate2);
}


TEST(CandidateRegionEqualityTest, differentCandidateRegionsReturnsFalse) {
    const CandidateRegion candidate1 { Interval(0, 2), "test" };
    const CandidateRegion candidate2 { Interval(0, 1), "test" };

    EXPECT_FALSE(candidate1 == candidate2);
}


TEST(CandidateRegionInequalityTest, identicalCandidateRegionsReturnsFalse) {
    const CandidateRegion candidate1 { Interval(0, 2), "test" };
    const CandidateRegion candidate2 { Interval(0, 2), "test" };

    EXPECT_FALSE(candidate1 != candidate2);
}


TEST(CandidateRegionInequalityTest, differentCandidateRegionsReturnsTrue) {
    const CandidateRegion candidate1 { Interval(0, 2), "test" };
    const CandidateRegion candidate2 { Interval(0, 2), "test2" };

    EXPECT_TRUE(candidate1 != candidate2);
}


TEST(IdentifyLowCoverageIntervalsTest, emptyCovgsReturnEmpty) {
    const auto min_len { 5 };
    const auto min_covg { 0 };
    const std::vector<uint32_t> covgs;

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, singleCovgPositionAboveThresholdReturnEmpty) {
    const auto min_len { 1 };
    const auto min_covg { 1 };
    const std::vector<uint32_t> covgs { 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, singleCovgPositionBelowThresholdReturnSingleInterval) {
    const auto min_len { 1 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(0, 1) };

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, allPositionsAboveThresholdReturnEmpty) {
    const auto min_len { 1 };
    const auto min_covg { 1 };
    const std::vector<uint32_t> covgs { 2, 2, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, allPositionsBelowThresholdReturnIntervalForWholeVector) {
    const auto min_len { 1 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2, 2, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(0, 4) };

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, allPositionsBelowThresholdButLessThanMinLengthReturnEmpty) {
    const auto min_len { 10 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2, 2, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, twoPositionsAtStartBelowThresholdButLessThanMinLenReturnEmpty) {
    const auto min_len { 3 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2, 2, 4, 4, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, twoPositionsInMiddleBelowThresholdButLessThanMinLenReturnEmpty) {
    const auto min_len { 3 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 2, 2, 4, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, twoPositionsAtEndBelowThresholdButLessThanMinLenReturnEmpty) {
    const auto min_len { 3 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 4, 4, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, twoPositionsAtStartBelowThresholdReturnSingleInterval) {
    const auto min_len { 2 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2, 2, 4, 4, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(0, 2) };

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, twoPositionsInMiddleBelowThresholdReturnSingleInterval) {
    const auto min_len { 1 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 2, 2, 4, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(1, 3) };

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, twoPositionsAtEndBelowThresholdReturnSingleInterval) {
    const auto min_len { 2 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 4, 4, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(3, 5) };

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, twoRegionsAtStartAndEndBelowThresholdReturnTwoIntervals) {
    const auto min_len { 2 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 2, 2, 4, 4, 4, 2, 2 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(0, 2), Interval(5, 7) };

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, twoRegionsBelowThresholdReturnTwoIntervals) {
    const auto min_len { 2 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 2, 1, 1, 4, 1, 2, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(1, 4), Interval(5, 7) };

    EXPECT_EQ(actual, expected);
}


TEST(IdentifyLowCoverageIntervalsTest, twoRegionsBelowThresholdOneLessThanMinLengthReturnOneInterval) {
    const auto min_len { 3 };
    const auto min_covg { 3 };
    const std::vector<uint32_t> covgs { 4, 2, 1, 1, 4, 1, 2, 4 };

    const auto actual { identify_low_coverage_intervals(covgs, min_covg, min_len) };
    const std::vector<Interval> expected { Interval(1, 4) };

    EXPECT_EQ(actual, expected);
}


TEST(FindCandidateRegionsForPanNodeTest, emptyPanNodeReturnsNoCandidates) {
    const auto num_samples { 1 };
    const auto prg_id { 3 };
    LocalPRG local_prg { prg_id, "test", "" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path;
    const std::vector<KmerNodePtr> kmer_node_max_likelihood_path;

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(0, prg_id, "test") };
    pangraph_node->kmer_prg = local_prg.kmer_prg;
    pangraph_node->kmer_prg.setup_coverages(num_samples);

    auto local_prg_ptr { std::make_shared<LocalPRG>(local_prg) };
    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr, kmer_node_max_likelihood_path,
                                                local_node_max_likelihood_path };

    const CandidateRegions expected;
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components) };

    EXPECT_EQ(actual, expected);

}


TEST(FindCandidateRegionsForPanNodeTest, noCoverageReturnWholePrgAsCandidate) {
    const auto num_samples { 1 };
    const auto prg_id { 3 };
    LocalPRG local_prg { prg_id, "test", "AAA 5 G 6 C 5 TTT" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg.minimizer_sketch(index, w, k);
    const std::string expected_sequence { "AAAGTTT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path {
            local_prg.kmernode_path_from_localnode_path(local_node_max_likelihood_path) };

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(0, prg_id, "test") };
    pangraph_node->kmer_prg = local_prg.kmer_prg;
    pangraph_node->kmer_prg.setup_coverages(num_samples);

    auto local_prg_ptr { std::make_shared<LocalPRG>(local_prg) };
    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr, kmer_node_max_likelihood_path,
                                                local_node_max_likelihood_path };

    CandidateRegion expected_candidate { Interval(0, 7), pangraph_node->get_name() };
    expected_candidate.max_likelihood_sequence = expected_sequence;
    const CandidateRegions expected { std::make_pair(expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components) };
    const auto actual_sequence { actual.at(expected_candidate.get_id()).max_likelihood_sequence };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
}


TEST(FindCandidateRegionsForPanNodeTest, highCoverageReturnEmpty) {
    const auto num_samples { 1 };
    const auto prg_id { 3 };
    LocalPRG local_prg { prg_id, "test", "AAA 5 G 6 C 5 TTT" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg.minimizer_sketch(index, w, k);
    const std::string expected_sequence { "AAAGTTT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path {
            local_prg.kmernode_path_from_localnode_path(local_node_max_likelihood_path) };

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(0, prg_id, "test") };
    pangraph_node->kmer_prg = local_prg.kmer_prg;
    pangraph_node->kmer_prg.setup_coverages(num_samples);
    for (const auto &kmer_node : pangraph_node->kmer_prg.nodes) {
        kmer_node->set_covg(100, false, 0);
    }

    auto local_prg_ptr { std::make_shared<LocalPRG>(local_prg) };
    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr, kmer_node_max_likelihood_path,
                                                local_node_max_likelihood_path };

    const CandidateRegions expected;
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components) };

    EXPECT_EQ(actual, expected);
}


TEST(FindCandidateRegionsForPanNodeTest, noCoverageOnFiveBasesReturnFiveBasesAsCandidate) {
    const auto num_samples { 1 };
    const auto prg_id { 3 };
    LocalPRG local_prg { prg_id, "test", "AAAA 5 GGG 6 CCC 5 TTTT" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg.minimizer_sketch(index, w, k);
    const std::string expected_sequence { "AAAAGGGTTTT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };
    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 3, 5, 7, 9, 11, 13, 14, 15 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());
    for (const auto &idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg.kmer_prg.nodes[idx]);
    }

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(0, prg_id, "test") };
    local_prg.kmer_prg.setup_coverages(num_samples);
    const std::vector<int> kmer_node_idxs_for_high_covg { 0, 1, 14, 15 };
    for (const auto &idx : kmer_node_idxs_for_high_covg) {
        local_prg.kmer_prg.nodes[idx]->set_covg(100, false, 0);
    }
    pangraph_node->kmer_prg = local_prg.kmer_prg;

    auto local_prg_ptr { std::make_shared<LocalPRG>(local_prg) };
    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr, kmer_node_max_likelihood_path,
                                                local_node_max_likelihood_path };

    CandidateRegion expected_candidate { Interval(3, 8), pangraph_node->get_name() };
    const CandidateRegions expected { std::make_pair(expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components) };
    const auto actual_sequence { actual.at(expected_candidate.get_id()).get_max_likelihood_sequence_with_flanks() };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
}


TEST(FindCandidateRegionsForPanNodeTest, noCoverageOnFiveBasesReturnFiveBasesPlusPaddingAsCandidate) {
    const auto prg_id { 3 };
    LocalPRG local_prg { prg_id, "test", "AAAA 5 GGG 6 CCC 5 TTTT" };
    const std::string expected_sequence { "AAAAGGGTTTT" };
    const auto expected_max_likelihood_sequence { "AAGGGTT" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg.minimizer_sketch(index, w, k);

    const std::vector<int> local_node_idxs_for_max_path { 0, 1, 3 };
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    local_node_max_likelihood_path.reserve(local_node_idxs_for_max_path.size());

    for (const auto &idx : local_node_idxs_for_max_path) {
        local_node_max_likelihood_path.push_back(local_prg.prg.nodes[idx]);
    }

    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 3, 5, 7, 9, 11, 13, 14, 15 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());

    for (const auto &idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg.kmer_prg.nodes[idx]);
    }

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(0, prg_id, "test") };
    const auto num_samples { 1 };
    local_prg.kmer_prg.setup_coverages(num_samples);
    const std::vector<int> kmer_node_idxs_for_high_covg { 0, 1, 14, 15 };

    for (const auto &idx : kmer_node_idxs_for_high_covg) {
        local_prg.kmer_prg.nodes[idx]->set_covg(100, false, 0);
    }

    pangraph_node->kmer_prg = local_prg.kmer_prg;

    auto local_prg_ptr { std::make_shared<LocalPRG>(local_prg) };
    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr, kmer_node_max_likelihood_path,
                                                local_node_max_likelihood_path };

    const auto interval_padding { 1 };
    CandidateRegion expected_candidate { Interval(3, 8), pangraph_node->get_name(), interval_padding };
    const CandidateRegions expected { std::make_pair(expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components, interval_padding) };
    const auto actual_sequence { actual.at(expected_candidate.get_id()).get_max_likelihood_sequence_with_flanks() };
    const auto actual_max_likelihood_sequence { actual.at(expected_candidate.get_id()).max_likelihood_sequence };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
    EXPECT_EQ(actual_max_likelihood_sequence, expected_max_likelihood_sequence);
}


TEST(FindCandidateRegionsForPanNodeTest, noCoverageOnStartFiveBasesReturnFiveBasesAsCandidate) {
    const auto prg_id { 3 };
    LocalPRG local_prg { prg_id, "test", "AAAA 5 GGG 6 CCC 5 TTTT" };
    const std::string expected_sequence { "AAAAGGGTTTT" };
    const auto expected_max_likelihood_sequence { "AAAAG" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg.minimizer_sketch(index, w, k);

    const std::vector<int> local_node_idxs_for_max_path { 0, 1, 3 };
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    local_node_max_likelihood_path.reserve(local_node_idxs_for_max_path.size());

    for (const auto &idx : local_node_idxs_for_max_path) {
        local_node_max_likelihood_path.push_back(local_prg.prg.nodes[idx]);
    }

    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 3, 5, 7, 9, 11, 13, 14, 15 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());

    for (const auto &idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg.kmer_prg.nodes[idx]);
    }

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(0, prg_id, "test") };
    const auto num_samples { 1 };
    local_prg.kmer_prg.setup_coverages(num_samples);
    const std::vector<int> kmer_node_idxs_for_high_covg { 9, 11, 13, 14, 15 };

    for (const auto &idx : kmer_node_idxs_for_high_covg) {
        local_prg.kmer_prg.nodes[idx]->set_covg(100, false, 0);
    }

    pangraph_node->kmer_prg = local_prg.kmer_prg;

    auto local_prg_ptr { std::make_shared<LocalPRG>(local_prg) };
    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr, kmer_node_max_likelihood_path,
                                                local_node_max_likelihood_path };

    const auto interval_padding { 0 };
    CandidateRegion expected_candidate { Interval(0, 5), pangraph_node->get_name(), interval_padding };
    const CandidateRegions expected { std::make_pair(expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components, interval_padding) };
    const auto actual_sequence { actual.at(expected_candidate.get_id()).get_max_likelihood_sequence_with_flanks() };
    const auto actual_max_likelihood_sequence { actual.at(expected_candidate.get_id()).max_likelihood_sequence };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
    EXPECT_EQ(actual_max_likelihood_sequence, expected_max_likelihood_sequence);
}


TEST(FindCandidateRegionsForPanNodeTest, noCoverageOnEndFiveBasesReturnFiveBasesAsCandidate) {
    const auto prg_id { 3 };
    LocalPRG local_prg { prg_id, "test", "AAAA 5 GGG 6 CCC 5 TTTT" };
    const std::string expected_sequence { "AAAAGGGTTTT" };
    const auto expected_max_likelihood_sequence { "GTTTT" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg.minimizer_sketch(index, w, k);

    const std::vector<int> local_node_idxs_for_max_path { 0, 1, 3 };
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    local_node_max_likelihood_path.reserve(local_node_idxs_for_max_path.size());

    for (const auto &idx : local_node_idxs_for_max_path) {
        local_node_max_likelihood_path.push_back(local_prg.prg.nodes[idx]);
    }

    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 3, 5, 7, 9, 11, 13, 14, 15 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());

    for (const auto &idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg.kmer_prg.nodes[idx]);
    }

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(0, prg_id, "test") };
    const auto num_samples { 1 };
    local_prg.kmer_prg.setup_coverages(num_samples);
    const std::vector<int> kmer_node_idxs_for_high_covg { 0, 1, 2, 3, 5 };

    for (const auto &idx : kmer_node_idxs_for_high_covg) {
        local_prg.kmer_prg.nodes[idx]->set_covg(100, false, 0);
    }

    pangraph_node->kmer_prg = local_prg.kmer_prg;

    auto local_prg_ptr { std::make_shared<LocalPRG>(local_prg) };
    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr, kmer_node_max_likelihood_path,
                                                local_node_max_likelihood_path };

    const auto interval_padding { 0 };
    CandidateRegion expected_candidate { Interval(6, 11), pangraph_node->get_name(), interval_padding };
    const CandidateRegions expected { std::make_pair(expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components, interval_padding) };
    const auto actual_sequence { actual.at(expected_candidate.get_id()).get_max_likelihood_sequence_with_flanks() };
    const auto actual_max_likelihood_sequence { actual.at(expected_candidate.get_id()).max_likelihood_sequence };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
    EXPECT_EQ(actual_max_likelihood_sequence, expected_max_likelihood_sequence);
}


TEST(FindCandidateRegionsForPanNodeTest,
     noCoverageOnFiveBasesWithinDoubleNestingReturnFiveBasesPlusPaddingAsCandidate) {
    const auto prg_id { 3 };
    LocalPRG local_prg { prg_id, "test", "AAAA 5 CCCC 6 GG 7 XXX 8 YYY 7 GG 5 TTTT" };
    const std::string expected_sequence { "AAAAGGYYYGGTTTT" };
    const auto expected_max_likelihood_sequence { "GGYYYGG" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg.minimizer_sketch(index, w, k);

    const std::vector<int> local_node_idxs_for_max_path { 0, 2, 4, 5, 6 };
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    local_node_max_likelihood_path.reserve(local_node_idxs_for_max_path.size());

    for (const auto &idx : local_node_idxs_for_max_path) {
        local_node_max_likelihood_path.push_back(local_prg.prg.nodes[idx]);
    }

    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 4, 6, 9, 12, 15, 18, 21, 23, 24, 19, 22, 25 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());

    for (const auto &idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg.kmer_prg.nodes[idx]);
    }

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(0, prg_id, "test") };
    const auto num_samples { 1 };
    local_prg.kmer_prg.setup_coverages(num_samples);
    const std::vector<int> kmer_node_idxs_for_high_covg { 0, 1, 2, 4, 24, 19, 22, 25 };

    for (const auto &idx : kmer_node_idxs_for_high_covg) {
        local_prg.kmer_prg.nodes[idx]->set_covg(100, false, 0);
    }

    pangraph_node->kmer_prg = local_prg.kmer_prg;

    auto local_prg_ptr { std::make_shared<LocalPRG>(local_prg) };
    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr, kmer_node_max_likelihood_path,
                                                local_node_max_likelihood_path };

    const auto interval_padding { 1 };
    CandidateRegion expected_candidate { Interval(5, 10), pangraph_node->get_name(), interval_padding };
    const CandidateRegions expected { std::make_pair(expected_candidate.get_id(), expected_candidate) };
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components, interval_padding) };
    const auto actual_sequence { actual.at(expected_candidate.get_id()).get_max_likelihood_sequence_with_flanks() };
    const auto actual_max_likelihood_sequence { actual.at(expected_candidate.get_id()).max_likelihood_sequence };

    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual_sequence, expected_sequence);
    EXPECT_EQ(actual_max_likelihood_sequence, expected_max_likelihood_sequence);
}


TEST(FindCandidateRegionsForPanNodeTest, noCoverageOnTwoFiveBaseRegionsWithinDoubleNestingReturnTwoRegionsAsCandidate) {
    const auto prg_id { 3 };
    LocalPRG local_prg { prg_id, "test", "AAAA 5 CCCC 6 GG 7 XXX 8 YYY 7 GG 5 TTTTT" };
    const std::string expected_sequence { "AAAAGGYYYGGTTTTT" };
    std::vector<std::string> expected_max_likelihood_sequences { "AAAAG", "YGGTT" };
    auto index { std::make_shared<Index>() };
    const auto w { 1 };
    const auto k { 3 };
    local_prg.minimizer_sketch(index, w, k);

    const std::vector<int> local_node_idxs_for_max_path { 0, 2, 4, 5, 6 };
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    local_node_max_likelihood_path.reserve(local_node_idxs_for_max_path.size());

    for (const auto &idx : local_node_idxs_for_max_path) {
        local_node_max_likelihood_path.push_back(local_prg.prg.nodes[idx]);
    }

    const std::vector<int> kmer_node_idxs_for_max_path { 0, 1, 2, 4, 6, 9, 12, 15, 18, 21, 23, 25, 19, 22, 24, 26 };
    std::vector<KmerNodePtr> kmer_node_max_likelihood_path;
    kmer_node_max_likelihood_path.reserve(kmer_node_idxs_for_max_path.size());

    for (const auto &idx : kmer_node_idxs_for_max_path) {
        kmer_node_max_likelihood_path.push_back(local_prg.kmer_prg.nodes[idx]);
    }

    PanNodePtr pangraph_node { std::make_shared<pangenome::Node>(0, prg_id, "test") };
    const auto num_samples { 1 };
    local_prg.kmer_prg.setup_coverages(num_samples);
    const std::vector<int> kmer_node_idxs_for_high_covg { 12, 24, 26 };

    for (const auto &idx : kmer_node_idxs_for_high_covg) {
        local_prg.kmer_prg.nodes[idx]->set_covg(100, false, 0);
    }

    pangraph_node->kmer_prg = local_prg.kmer_prg;

    auto local_prg_ptr { std::make_shared<LocalPRG>(local_prg) };
    const TmpPanNode pangraph_node_components { pangraph_node, local_prg_ptr, kmer_node_max_likelihood_path,
                                                local_node_max_likelihood_path };

    const auto interval_padding { 0 };
    CandidateRegion expected_candidate1 { Interval(0, 5), pangraph_node->get_name(), interval_padding };
    CandidateRegion expected_candidate2 { Interval(8, 13), pangraph_node->get_name(), interval_padding };
    const CandidateRegions expected { std::make_pair(expected_candidate1.get_id(), expected_candidate1),
                                      std::make_pair(expected_candidate2.get_id(), expected_candidate2) };
    const auto actual { find_candidate_regions_for_pan_node(pangraph_node_components, interval_padding) };

    EXPECT_EQ(actual, expected);
    std::vector<std::string> actual_max_likelihood_sequences;
    actual_max_likelihood_sequences.reserve(2);
    for (const auto &actual_candidate : actual) {
        EXPECT_EQ(actual_candidate.second.get_max_likelihood_sequence_with_flanks(), expected_sequence);
        actual_max_likelihood_sequences.push_back(actual_candidate.second.max_likelihood_sequence);
    }
    std::sort(actual_max_likelihood_sequences.begin(), actual_max_likelihood_sequences.end());
    std::sort(expected_max_likelihood_sequences.begin(), expected_max_likelihood_sequences.end());
    EXPECT_EQ(actual_max_likelihood_sequences, expected_max_likelihood_sequences);
}


TEST(GenerateReadPileupsForCandidateRegionTest, emptyReadsFileWithReadCoordsPileupsEmpty) {

    Fastaq temp_fastq { false, true };
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate { Interval(1, 3), "test" };

    ReadCoordinate read_coord_1 { 0, 6, 10, true };
    candidate.read_coordinates.insert(read_coord_1);

    candidate.generate_read_pileup(temp_reads_filepath);
    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto &actual { candidate.pileup };
    const ReadPileup expected;

    EXPECT_EQ(actual, expected);
}


TEST(GenerateReadPileupsForCandidateRegionTest, nonEmptyReadsFileWithNoReadCoordsPileupsEmpty) {

    Fastaq temp_fastq { false, true };
    const auto read_name { "0" };
    const auto read_sequence { "ABC" };
    const std::vector<uint32_t> read_covg { 1, 2, 3 };
    const auto global_covg { 2 };
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate { Interval(1, 3), "test" };

    candidate.generate_read_pileup(temp_reads_filepath);
    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto &actual { candidate.pileup };
    const ReadPileup expected;

    EXPECT_EQ(actual, expected);
}


TEST(GenerateReadPileupsForCandidateRegionTest, oneReadInFileOneReadCoordPileupHasOneEntry) {

    Fastaq temp_fastq { false, true };
    const auto read_name { "0" };
    const std::string read_sequence { "XXXFOOXXX" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord { 0, 3, 6, true }; // FOO
    candidate.read_coordinates.insert(read_coord);

    candidate.generate_read_pileup(temp_reads_filepath);
    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto &actual { candidate.pileup };
    const ReadPileup expected { "FOO" };

    EXPECT_EQ(actual, expected);
}


TEST(GenerateReadPileupsForCandidateRegionTest, oneReadInFileOneWholeReadCoordPileupHasOneEntrySameAsRead) {

    Fastaq temp_fastq { false, true };
    const auto read_name { "0" };
    const std::string read_sequence { "XXXFOOXXX" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord { 0, 0, 10, true }; // XXXFOOXXX
    candidate.read_coordinates.insert(read_coord);

    candidate.generate_read_pileup(temp_reads_filepath);
    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto &actual { candidate.pileup };
    const ReadPileup expected { "XXXFOOXXX" };

    EXPECT_EQ(actual, expected);
}


TEST(GenerateReadPileupsForCandidateRegionTest,
     oneReadInFileOneReadCoordThatRunsPastEndOfReadPileupHasOneEntryUpToEndOfRead) {

    Fastaq temp_fastq { false, true };
    const auto read_name { "0" };
    const std::string read_sequence { "XXXFOOXXX" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord { 0, 5, 20, true }; // OXXX
    candidate.read_coordinates.insert(read_coord);

    candidate.generate_read_pileup(temp_reads_filepath);
    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto &actual { candidate.pileup };
    const ReadPileup expected { "OXXX" };

    EXPECT_EQ(actual, expected);
}


TEST(GenerateReadPileupsForCandidateRegionTest,
     oneReadInFileOneReverseReadCoordThatRunsPastEndOfReadPileupHasOneEntryUpToEndOfRead) {

    Fastaq temp_fastq { false, true };
    const auto read_name { "0" };
    const std::string read_sequence { "AATTCCGG" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord { 0, 5, 20, false }; // CCG
    candidate.read_coordinates.insert(read_coord);

    candidate.generate_read_pileup(temp_reads_filepath);
    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto &actual { candidate.pileup };
    const ReadPileup expected { "CCG" };

    EXPECT_EQ(actual, expected);
}


TEST(GenerateReadPileupsForCandidateRegionTest, oneReadInFileOneReadCoordOutsideReadPileupEmpty) {

    Fastaq temp_fastq { false, true };
    const auto read_name { "0" };
    const std::string read_sequence { "XXXFOOXXX" };
    const auto global_covg { 2 };
    const std::vector<uint32_t> read_covg(read_sequence.length(), global_covg);
    temp_fastq.add_entry(read_name, read_sequence, read_covg, global_covg);
    const fs::path temp_reads_filepath { fs::unique_path() };
    temp_fastq.save(temp_reads_filepath.string());

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord { 0, 15, 20, true };
    candidate.read_coordinates.insert(read_coord);

    candidate.generate_read_pileup(temp_reads_filepath);
    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto &actual { candidate.pileup };
    const ReadPileup expected;

    EXPECT_EQ(actual, expected);
}


TEST(GenerateReadPileupsForCandidateRegionTest, twoReadsInFileOneReverseOneForwardReadCoordPileupHasTwoEntries) {

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

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord1 { 0, 2, 4, true }; // TT
    ReadCoordinate read_coord2 { 1, 3, 6, false }; // GTA
    candidate.read_coordinates.insert(read_coord1);
    candidate.read_coordinates.insert(read_coord2);

    candidate.generate_read_pileup(temp_reads_filepath);
    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    const auto &actual { candidate.pileup };
    const ReadPileup expected { "TT", "GTA" };

    EXPECT_EQ(actual, expected);
}


TEST(GenerateReadPileupsForCandidateRegionTest, twoReadsInFileThreeForwardReadCoordsPileupHasThreeEntries) {

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

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord1 { 0, 0, 2, true }; // AA
    ReadCoordinate read_coord2 { 1, 4, 6, true }; // AC
    ReadCoordinate read_coord3 { 0, 1, 6, true }; // ATTCC
    candidate.read_coordinates.insert(read_coord1);
    candidate.read_coordinates.insert(read_coord2);
    candidate.read_coordinates.insert(read_coord3);

    candidate.generate_read_pileup(temp_reads_filepath);
    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    auto actual { candidate.pileup };
    std::sort(actual.begin(), actual.end());
    const ReadPileup expected { "AA", "AC", "ATTCC" };

    EXPECT_EQ(actual, expected);
}


TEST(GenerateReadPileupsForCandidateRegionTest, twoReadsInFileThreeForwardTwoTheSameReadCoordsPileupHasTwoEntries) {

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

    CandidateRegion candidate { Interval(0, 3), "test" };
    ReadCoordinate read_coord1 { 0, 0, 2, true }; // AA
    ReadCoordinate read_coord2 { 1, 4, 6, true }; // AC
    ReadCoordinate read_coord3 { 0, 0, 2, true }; // AA duplicate
    candidate.read_coordinates.insert(read_coord1);
    candidate.read_coordinates.insert(read_coord2);
    candidate.read_coordinates.insert(read_coord3);

    candidate.generate_read_pileup(temp_reads_filepath);
    const auto temp_removed_successfully { fs::remove(temp_reads_filepath) };
    ASSERT_TRUE(temp_removed_successfully);

    auto actual { candidate.pileup };
    std::sort(actual.begin(), actual.end());
    const ReadPileup expected { "AA", "AC" };

    EXPECT_EQ(actual, expected);
}


std::string read_file_to_string(const fs::path &filepath) {
    fs::ifstream f(filepath);
    std::stringstream buffer;
    buffer << f.rdbuf();
    return buffer.str();
}


TEST(WriteDenovoPathsToFileTest, noReadsDoesntWriteFile) {
    const std::string name { "test" };
    CandidateRegion candidate { Interval(0, 1), name };
    const fs::path output_directory { "/tmp" };
    const fs::path filepath { output_directory / candidate.filename };

    candidate.write_denovo_paths_to_file(output_directory);

    EXPECT_FALSE(fs::exists(filepath));
}


TEST(WriteDenovoPathsToFileTest, twoReadsWritesTwoReadsToFile) {
    const std::string name { "test" };
    CandidateRegion candidate { Interval(0, 1), name };
    const fs::path output_directory { "/tmp" };
    const fs::path filepath { output_directory / candidate.filename };
    const DenovoPaths paths { "shrubberies", "ni" };
    candidate.denovo_paths = paths;

    candidate.write_denovo_paths_to_file(output_directory);

    std::stringstream expected_ss;

    for (size_t i = 0; i < paths.size(); ++i) {
        expected_ss << ">test." << std::to_string(i) << std::endl << paths.at(i) << std::endl;
    }

    const std::string expected { expected_ss.str() };
    const auto actual { read_file_to_string(filepath) };

    fs::remove(filepath);

    EXPECT_EQ(actual, expected);
}
