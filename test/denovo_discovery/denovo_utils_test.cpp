#include <iostream>
#include <set>
#include "gtest/gtest.h"
#include "../test_macro.cpp"
#include "denovo_discovery/denovo_utils.h"
#include "minihit.h"
#include "pangenome/panread.h"


using PanNodePtr = std::shared_ptr<pangenome::Node>;
using PanReadPtr = std::shared_ptr<pangenome::Read>;


TEST(PathComponentsEqivalenceOperatorTest, twoEqualPathComponentsReturnsTrue) {
    prg::Path flank_left;
    flank_left.initialize(Interval(0, 2));
    prg::Path slice;
    slice.initialize(Interval(3, 6));
    prg::Path flank_right;
    flank_right.initialize(Interval(7, 8));
    const PathComponents x { flank_left, slice, flank_right };
    const PathComponents y { flank_left, slice, flank_right };

    EXPECT_TRUE(x == y);
}


TEST(PathComponentsEqivalenceOperatorTest, twoDifferentPathComponentsReturnsFalse) {
    prg::Path flank_left_x;
    flank_left_x.initialize(Interval(0, 1));
    prg::Path flank_left_y;
    flank_left_x.initialize(Interval(0, 2));
    prg::Path slice;
    slice.initialize(Interval(3, 6));
    prg::Path flank_right;
    flank_right.initialize(Interval(7, 8));
    const PathComponents x { flank_left_x, slice, flank_right };
    const PathComponents y { flank_left_y, slice, flank_right };

    EXPECT_FALSE(x == y);
}


TEST(PathComponentsNonEqivalenceOperatorTest, twoEqualPathComponentsReturnsFalse) {
    prg::Path flank_left;
    flank_left.initialize(Interval(0, 2));
    prg::Path slice;
    slice.initialize(Interval(3, 6));
    prg::Path flank_right;
    flank_right.initialize(Interval(7, 8));
    const PathComponents x { flank_left, slice, flank_right };
    const PathComponents y { flank_left, slice, flank_right };

    EXPECT_FALSE(x != y);
}


TEST(PathComponentsNonEqivalenceOperatorTest, twoDifferentPathComponentsReturnsTrue) {
    prg::Path flank_left_x;
    flank_left_x.initialize(Interval(0, 1));
    prg::Path flank_left_y;
    flank_left_x.initialize(Interval(0, 2));
    prg::Path slice;
    slice.initialize(Interval(3, 6));
    prg::Path flank_right;
    flank_right.initialize(Interval(7, 8));
    const PathComponents x { flank_left_x, slice, flank_right };
    const PathComponents y { flank_left_y, slice, flank_right };

    EXPECT_TRUE(x != y);
}


TEST(FindIntervalInLocalPathTest, emptyIntervalReturnsEmptyResult) {
    const Interval interval;
    LocalPRG local_prg { 0, "test", "A" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0] };

    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };
    const PathComponents expected;

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, emptyPrgReturnsEmptyResult) {
    const Interval interval { 0, 5 };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path;

    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };
    const PathComponents expected;

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, emptyInputsReturnsEmptyResult) {
    const Interval interval;
    const std::vector<LocalNodePtr> local_node_max_likelihood_path;

    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };
    const PathComponents expected;

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, singleBaseIntervalSingleBasePrgReturnsSingleBase) {
    const Interval interval { 0, 1 };
    LocalPRG local_prg { 0, "test", "A" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0] };

    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };
    prg::Path expected_slice;
    expected_slice.initialize(interval);
    const PathComponents expected { prg::Path(), expected_slice, prg::Path() };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, singleBaseIntervalMultiBasePrgReturnsSingleBaseAndRightFlank) {
    const Interval interval { 0, 1 };
    LocalPRG local_prg { 0, "test", "AT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0] };

    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };
    prg::Path expected_slice;
    expected_slice.initialize(interval);
    prg::Path expected_right_flank;
    expected_right_flank.initialize(Interval(1, 2));
    const PathComponents expected { prg::Path(), expected_slice, expected_right_flank };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, singleBaseIntervalMultiBasePrgReturnsSingleBaseAndLeftFlank) {
    const Interval interval { 1, 2 };
    LocalPRG local_prg { 0, "test", "AT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0] };

    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };
    prg::Path expected_slice;
    expected_slice.initialize(interval);
    prg::Path expected_left_flank;
    expected_left_flank.initialize(Interval(0, 1));
    const PathComponents expected { expected_left_flank, expected_slice, prg::Path() };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, singleBaseIntervalMultiBasePrgReturnsSingleBaseAndBothFlanks) {
    const Interval interval { 1, 2 };
    LocalPRG local_prg { 0, "test", "TAT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(interval);
    expected_left_flank.initialize(Interval(0, 1));
    expected_right_flank.initialize(Interval(2, 3));
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, singleBaseIntervalMultiNodeSingleBasePrgReturnsSingleBaseAndBothFlanks) {
    const Interval interval { 1, 2 };
    LocalPRG local_prg { 0, "test", "T 5 A 6 C 5 T" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(4, 5));  // the A in the PRG
    expected_left_flank.initialize(Interval(0, 1));
    expected_right_flank.initialize(Interval(12, 13));
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, singleBaseIntervalMultiNodeSingleBasePrgReturnsSingleBaseAndLeftFlank) {
    const Interval interval { 2, 3 };
    LocalPRG local_prg { 0, "test", "T 5 A 6 C 5 T" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(12, 13));  // the T at the end of the PRG
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 1), Interval(4, 5) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, singleBaseIntervalMultiNodeSingleBasePrgReturnsSingleBaseAndRightFlank) {
    const Interval interval { 0, 1 };
    LocalPRG local_prg { 0, "test", "T 5 A 6 C 5 T" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(0, 1));  // the T at the start of the PRG
    expected_right_flank.initialize(std::vector<Interval> { Interval(4, 5), Interval(12, 13) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, singleBaseIntervalMultiNodeMultiBasePrgReturnsSingleBaseAndRightFlank) {
    const Interval interval { 0, 1 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 TT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(0, 1));  // the T at the start of the PRG
    expected_right_flank.initialize(std::vector<Interval> { Interval(1, 2), Interval(5, 7), Interval(15, 17) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest,
     singleBaseIntervalMultiNodeMultiBasePrgReturnsSingleBaseSingleLeftFlankAndMultiRightFlank) {
    const Interval interval { 1, 2 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 TT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(1, 2));  // the second T at the start of the PRG
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 1) });
    expected_right_flank.initialize(std::vector<Interval> { Interval(5, 7), Interval(15, 17) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, singleBaseIntervalMultiNodeMultiBasePrgReturnsSingleBaseAndLeftFlank) {
    const Interval interval { 5, 6 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 TT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(16, 17));  // the T at the end of the PRG
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 2), Interval(5, 7), Interval(15, 16) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest,
     singleBaseIntervalMultiNodeMultiBasePrgReturnsSingleBaseSingleRightFlankAndMultiLeftFlank) {
    const Interval interval { 4, 5 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 TT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(15, 16));  // the second last T at the end of the PRG
    expected_right_flank.initialize(std::vector<Interval> { Interval(16, 17) });
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 2), Interval(5, 7) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, singleBaseIntervalMultiNodeMultiBasePrgReturnsSingleBaseAndFlanks) {
    const Interval interval { 2, 3 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 TT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[2],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(10, 11));  // the first C in the PRG
    expected_right_flank.initialize(std::vector<Interval> { Interval(11, 12), Interval(15, 17) });
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 2) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, multiBaseIntervalMultiNodeMultiBasePrgReturnsSingleNodeAndRightFlanks) {
    const Interval interval { 0, 2 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 TT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(0, 2));  // the TT at the start of the PRG
    expected_right_flank.initialize(std::vector<Interval> { Interval(5, 7), Interval(15, 17) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, multiBaseIntervalMultiNodeMultiBasePrgReturnsSingleNodeAndLeftFlanks) {
    const Interval interval { 4, 6 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 TT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(15, 17));  // the TT at the end of the PRG
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 2), Interval(5, 7) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, multiBaseIntervalMultiNodeMultiBasePrgReturnsSingleNodeAndBothFlanks) {
    const Interval interval { 2, 4 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 TT" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[2],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(Interval(10, 12));  // the CC in the PRG
    expected_right_flank.initialize(std::vector<Interval> { Interval(15, 17) });
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 2) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest,
     nodeSpanningIntervalMultiNodeMultiBasePrgReturnsMultiNodeSingleLeftAndMultiRightFlanks) {
    const Interval interval { 1, 3 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 GG" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[2],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(std::vector<Interval> { Interval(1, 2), Interval(10, 11) });  // TC
    expected_right_flank.initialize(std::vector<Interval> { Interval(11, 12), Interval(15, 17) });
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 1) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest,
     nodeSpanningIntervalMultiNodeMultiBasePrgReturnsMultiNodeMultiLeftAndSingleRightFlanks) {
    const Interval interval { 3, 5 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 GG" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[2],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(std::vector<Interval> { Interval(11, 12), Interval(15, 16) });  // CG
    expected_right_flank.initialize(std::vector<Interval> { Interval(16, 17) });
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 2), Interval(10, 11) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest,
     nodeSpanningIntervalMultiNodeMultiBasePrgReturnsMultiNodeSingleLeftAndSingleRightFlanks) {
    const Interval interval { 1, 5 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 GG" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[2],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(std::vector<Interval> { Interval(1, 2), Interval(10, 12), Interval(15, 16) });  // TCCG
    expected_right_flank.initialize(std::vector<Interval> { Interval(16, 17) });
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 1) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, intervalSpanningWholePrgReturnsNoFlanks) {
    const Interval interval { 0, 6 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 GG" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(std::vector<Interval> { Interval(0, 2), Interval(5, 7), Interval(15, 17) });  // TTAAGG
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindIntervalInLocalPathTest, intervalSpanningPastEndOfPrgReturnsUpToEndOfPrg) {
    const Interval interval { 2, 8 };
    LocalPRG local_prg { 0, "test", "TT 5 AA 6 CC 5 GG" };
    const std::vector<LocalNodePtr> local_node_max_likelihood_path { local_prg.prg.nodes[0], local_prg.prg.nodes[1],
                                                                     local_prg.prg.nodes[3] };

    prg::Path expected_slice;
    prg::Path expected_left_flank;
    prg::Path expected_right_flank;
    expected_slice.initialize(std::vector<Interval> { Interval(5, 7), Interval(15, 17) });  // AAGG
    expected_left_flank.initialize(std::vector<Interval> { Interval(0, 2) });
    const PathComponents expected { expected_left_flank, expected_slice, expected_right_flank };
    const auto actual { find_interval_and_flanks_in_localpath(interval, local_node_max_likelihood_path) };

    EXPECT_EQ(actual, expected);
}


TEST(FindHitsInsidePathTest, emptyPathReturnsEmpty) {
    std::set<MinimizerHitPtr, pComp_path> hits;
    prg::Path local_path;

    const auto actual { find_hits_inside_path(hits, local_path) };
    const std::set<MinimizerHitPtr, pComp_path> expected;

    EXPECT_EQ(actual, expected);
}


TEST(FindHitsInsidePathTest, hitNotOnPathReturnEmpty) {
    LocalPRG local_prg(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    const std::vector<LocalNodePtr> local_max_likelihood_path { local_prg.prg.nodes[1], local_prg.prg.nodes[2],
                                                                local_prg.prg.nodes[4], local_prg.prg.nodes[6],
                                                                local_prg.prg.nodes[7] };  // A G C T CGG  TAT
    const uint32_t read_id { 0 };
    const uint32_t prg_id { 3 };
    const uint32_t knode_id { 0 };
    const Interval read_interval { 1, 4 };
    const bool is_forward { true };
    std::deque<Interval> intervals { Interval(7, 8), Interval(10, 12) };
    prg::Path prg_path;
    prg_path.initialize(intervals);

    std::set<MinimizerHitPtr, pComp_path> read_hits;
    MinimizerHitPtr minimizer_hit {
            std::make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, is_forward) };
    read_hits.insert(minimizer_hit);
    prg::Path local_path;
    for (const auto &node : local_max_likelihood_path) {
        local_path.add_end_interval(node->pos);
    }
    std::set<MinimizerHitPtr, pComp_path> actual { find_hits_inside_path(read_hits, local_path) };
    std::set<MinimizerHitPtr, pComp_path> expected;

    EXPECT_EQ(actual, expected);
}


TEST(FindHitsInsidePathTest, hitsBranchingFromPathReturnEmpty) {
    LocalPRG local_prg(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    const std::vector<LocalNodePtr> local_max_likelihood_path { local_prg.prg.nodes[1], local_prg.prg.nodes[2],
                                                                local_prg.prg.nodes[4], local_prg.prg.nodes[6],
                                                                local_prg.prg.nodes[7] };  // A G C T CGG  TAT
    const uint32_t read_id { 0 };
    const uint32_t prg_id { 3 };
    const uint32_t knode_id { 0 };
    const Interval read_interval { 1, 4 };
    const bool is_forward { true };
    std::deque<Interval> intervals { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
    prg::Path prg_path;
    prg_path.initialize(intervals);

    std::set<MinimizerHitPtr, pComp_path> read_hits;
    MinimizerHitPtr minimizer_hit {
            std::make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, is_forward) };
    read_hits.insert(minimizer_hit);

    intervals = { Interval(29, 30), Interval(31, 33) };
    prg_path.initialize(intervals);
    minimizer_hit = { std::make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, is_forward) };
    read_hits.insert(minimizer_hit);

    prg::Path local_path;
    for (const auto &node : local_max_likelihood_path) {
        local_path.add_end_interval(node->pos);
    }
    std::set<MinimizerHitPtr, pComp_path> actual { find_hits_inside_path(read_hits, local_path) };
    std::set<MinimizerHitPtr, pComp_path> expected;

    EXPECT_EQ(actual, expected);
}


TEST(FindHitsInsidePathTest, hitsOverlappingEdgesOfPathReturnEmpty) {
    LocalPRG local_prg(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    const std::vector<LocalNodePtr> local_max_likelihood_path { local_prg.prg.nodes[1], local_prg.prg.nodes[2],
                                                                local_prg.prg.nodes[4], local_prg.prg.nodes[6],
                                                                local_prg.prg.nodes[7] };  // A G C T CGG  TAT
    const uint32_t read_id { 0 };
    const uint32_t prg_id { 3 };
    const uint32_t knode_id { 0 };
    const Interval read_interval { 1, 4 };
    const bool is_forward { true };
    std::deque<Interval> intervals { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg::Path prg_path;
    prg_path.initialize(intervals);

    std::set<MinimizerHitPtr, pComp_path> read_hits;
    MinimizerHitPtr minimizer_hit {
            std::make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, is_forward) };
    read_hits.insert(minimizer_hit);

    intervals = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(intervals);
    minimizer_hit = { std::make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, is_forward) };
    read_hits.insert(minimizer_hit);

    intervals = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
    prg_path.initialize(intervals);
    minimizer_hit = { std::make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, is_forward) };
    read_hits.insert(minimizer_hit);

    prg::Path local_path;
    for (const auto &node : local_max_likelihood_path) {
        local_path.add_end_interval(node->pos);
    }
    std::set<MinimizerHitPtr, pComp_path> actual { find_hits_inside_path(read_hits, local_path) };
    std::set<MinimizerHitPtr, pComp_path> expected;

    EXPECT_EQ(actual, expected);
}


TEST(FindHitsInsidePathTest, hitsOnPathReturnCorrectHits) {
    LocalPRG local_prg(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    const std::vector<LocalNodePtr> local_max_likelihood_path { local_prg.prg.nodes[1], local_prg.prg.nodes[2],
                                                                local_prg.prg.nodes[4], local_prg.prg.nodes[6],
                                                                local_prg.prg.nodes[7] };  // A G C T CGG  TAT
    const uint32_t read_id { 0 };
    const uint32_t prg_id { 3 };
    const uint32_t knode_id { 0 };
    const Interval read_interval { 1, 4 };
    const bool is_forward { true };
    std::set<MinimizerHitPtr, pComp_path> expected;

    std::deque<Interval> intervals { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
    prg::Path prg_path;
    prg_path.initialize(intervals);

    std::set<MinimizerHitPtr, pComp_path> read_hits;
    MinimizerHitPtr minimizer_hit {
            std::make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, is_forward) };
    read_hits.insert(minimizer_hit);
    expected.insert(minimizer_hit);

    intervals = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(intervals);
    minimizer_hit = { std::make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, is_forward) };
    read_hits.insert(minimizer_hit);
    expected.insert(minimizer_hit);

    intervals = { Interval(16, 17), Interval(27, 29) };
    prg_path.initialize(intervals);
    minimizer_hit = { std::make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, is_forward) };
    read_hits.insert(minimizer_hit);
    expected.insert(minimizer_hit);

    intervals = { Interval(27, 30) };
    prg_path.initialize(intervals);
    minimizer_hit = { std::make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, is_forward) };
    read_hits.insert(minimizer_hit);
    expected.insert(minimizer_hit);

    prg::Path local_path;
    for (const auto &node : local_max_likelihood_path) {
        local_path.add_end_interval(node->pos);
    }
    std::set<MinimizerHitPtr, pComp_path> actual { find_hits_inside_path(read_hits, local_path) };

    EXPECT_EQ(actual, expected);
}


TEST(ExtractReadsTest, read_coordinate_ordering) {
    ReadCoordinate read_coord1 { 0, 6, 10, true };
    ReadCoordinate read_coord2 { 1, 16, 22, true };
    ReadCoordinate read_coord3 { 2, 3, 6, true };
    ReadCoordinate read_coord4 { 2, 8, 11, true };
    ReadCoordinate read_coord5 { 3, 2, 6, true };
    ReadCoordinate read_coord6 { 3, 6, 12, true };
    ReadCoordinate read_coord7 { 3, 2, 6, false };
    ReadCoordinate read_coord8 { 3, 6, 12, false };

    set<ReadCoordinate> read_coords = { read_coord8, read_coord7, read_coord6, read_coord5, read_coord4, read_coord3,
                                        read_coord2, read_coord1 };
    vector<ReadCoordinate> exp_read_coords = { read_coord1, read_coord2, read_coord3, read_coord4, read_coord5,
                                               read_coord7, read_coord6, read_coord8 };

    EXPECT_EQ(read_coords.size(), exp_read_coords.size());
    uint count = 0;
    for (const auto &r : read_coords) {
        if (exp_read_coords.size() <= count)
            break;
        EXPECT_EQ(r, exp_read_coords[count]);
        count++;
    }
}


TEST(ExtractReadsTest, read_coordinate_equals) {
    ReadCoordinate read_coord1 { 0, 6, 10, true };
    ReadCoordinate read_coord2 { 1, 16, 22, true };
    ReadCoordinate read_coord3 { 2, 3, 6, true };
    ReadCoordinate read_coord4 { 2, 8, 11, true };
    ReadCoordinate read_coord5 { 3, 2, 6, true };
    ReadCoordinate read_coord6 { 3, 6, 12, true };
    ReadCoordinate read_coord7 { 3, 2, 6, false };
    ReadCoordinate read_coord8 { 3, 6, 12, false };

    EXPECT_EQ(read_coord1, read_coord1);
    EXPECT_EQ(read_coord2, read_coord2);
    EXPECT_EQ(read_coord3, read_coord3);
    EXPECT_EQ(read_coord4, read_coord4);
    EXPECT_EQ(read_coord5, read_coord5);
    EXPECT_EQ(read_coord6, read_coord6);
    EXPECT_EQ(read_coord7, read_coord7);
    EXPECT_EQ(read_coord8, read_coord8);

    EXPECT_NE(read_coord1, read_coord2);
    EXPECT_NE(read_coord1, read_coord3);
    EXPECT_NE(read_coord1, read_coord4);
    EXPECT_NE(read_coord1, read_coord5);
    EXPECT_NE(read_coord1, read_coord6);
    EXPECT_NE(read_coord1, read_coord7);
    EXPECT_NE(read_coord1, read_coord8);

    EXPECT_NE(read_coord2, read_coord1);
    EXPECT_NE(read_coord2, read_coord3);
    EXPECT_NE(read_coord2, read_coord4);
    EXPECT_NE(read_coord2, read_coord5);
    EXPECT_NE(read_coord2, read_coord6);
    EXPECT_NE(read_coord2, read_coord7);
    EXPECT_NE(read_coord2, read_coord8);

    EXPECT_NE(read_coord3, read_coord2);
    EXPECT_NE(read_coord3, read_coord1);
    EXPECT_NE(read_coord3, read_coord4);
    EXPECT_NE(read_coord3, read_coord5);
    EXPECT_NE(read_coord3, read_coord6);
    EXPECT_NE(read_coord3, read_coord7);
    EXPECT_NE(read_coord3, read_coord8);

    EXPECT_NE(read_coord4, read_coord2);
    EXPECT_NE(read_coord4, read_coord3);
    EXPECT_NE(read_coord4, read_coord1);
    EXPECT_NE(read_coord4, read_coord5);
    EXPECT_NE(read_coord4, read_coord6);
    EXPECT_NE(read_coord4, read_coord7);
    EXPECT_NE(read_coord4, read_coord8);

    EXPECT_NE(read_coord5, read_coord2);
    EXPECT_NE(read_coord5, read_coord3);
    EXPECT_NE(read_coord5, read_coord4);
    EXPECT_NE(read_coord5, read_coord1);
    EXPECT_NE(read_coord5, read_coord6);
    EXPECT_NE(read_coord5, read_coord7);
    EXPECT_NE(read_coord5, read_coord8);

    EXPECT_NE(read_coord6, read_coord2);
    EXPECT_NE(read_coord6, read_coord3);
    EXPECT_NE(read_coord6, read_coord4);
    EXPECT_NE(read_coord6, read_coord5);
    EXPECT_NE(read_coord6, read_coord1);
    EXPECT_NE(read_coord6, read_coord7);
    EXPECT_NE(read_coord6, read_coord8);

    EXPECT_NE(read_coord7, read_coord2);
    EXPECT_NE(read_coord7, read_coord3);
    EXPECT_NE(read_coord7, read_coord4);
    EXPECT_NE(read_coord7, read_coord5);
    EXPECT_NE(read_coord7, read_coord6);
    EXPECT_NE(read_coord7, read_coord1);
    EXPECT_NE(read_coord7, read_coord8);

    EXPECT_NE(read_coord8, read_coord2);
    EXPECT_NE(read_coord8, read_coord3);
    EXPECT_NE(read_coord8, read_coord4);
    EXPECT_NE(read_coord8, read_coord5);
    EXPECT_NE(read_coord8, read_coord6);
    EXPECT_NE(read_coord8, read_coord7);
    EXPECT_NE(read_coord8, read_coord1);
}
