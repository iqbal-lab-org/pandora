#include <iostream>
#include <memory>
#include <set>
#include "gtest/gtest.h"
#include "../test_macro.cpp"
#include "denovo_discovery/extract_reads.h"
#include "interval.h"
#include "localPRG.h"
#include "minihit.h"
#include "minihits.h"
#include "prg/path.h"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"


using PanNodePtr =  std::shared_ptr<pangenome::Node>;
using PanReadPtr = std::shared_ptr<pangenome::Read>;
using std::make_pair;


TEST(PathComponentsEqivalenceOperatorTest, twoEqualPathComponents_returnTrue) {
    prg::Path flank_left;
    flank_left.initialize(Interval(0, 2));
    prg::Path slice;
    slice.initialize(Interval(3, 6));
    prg::Path flank_right;
    flank_right.initialize(Interval(7, 8));
    const PathComponents x {flank_left, slice, flank_right};
    const PathComponents y {flank_left, slice, flank_right};

    EXPECT_TRUE(x == y);
}

TEST(PathComponentsEqivalenceOperatorTest, twoDifferentPathComponents_returnFalse) {
    prg::Path flank_left_x;
    flank_left_x.initialize(Interval(0, 1));
    prg::Path flank_left_y;
    flank_left_x.initialize(Interval(0, 2));
    prg::Path slice;
    slice.initialize(Interval(3, 6));
    prg::Path flank_right;
    flank_right.initialize(Interval(7, 8));
    const PathComponents x {flank_left_x, slice, flank_right};
    const PathComponents y {flank_left_y, slice, flank_right};

    EXPECT_FALSE(x == y);
}

TEST(PathComponentsNonEqivalenceOperatorTest, twoEqualPathComponents_returnFalse) {
    prg::Path flank_left;
    flank_left.initialize(Interval(0, 2));
    prg::Path slice;
    slice.initialize(Interval(3, 6));
    prg::Path flank_right;
    flank_right.initialize(Interval(7, 8));
    const PathComponents x {flank_left, slice, flank_right};
    const PathComponents y {flank_left, slice, flank_right};

    EXPECT_FALSE(x != y);
}

TEST(PathComponentsNonEqivalenceOperatorTest, twoDifferentPathComponents_returnTrue) {
    prg::Path flank_left_x;
    flank_left_x.initialize(Interval(0, 1));
    prg::Path flank_left_y;
    flank_left_x.initialize(Interval(0, 2));
    prg::Path slice;
    slice.initialize(Interval(3, 6));
    prg::Path flank_right;
    flank_right.initialize(Interval(7, 8));
    const PathComponents x {flank_left_x, slice, flank_right};
    const PathComponents y {flank_left_y, slice, flank_right};

    EXPECT_TRUE(x != y);
}

TEST(IdentifyRegionsTest, identify_regions_null) {
    const auto min_len{5};
    const std::vector<uint32_t> covgs;
    auto low_covg_intervals{identify_regions(covgs, 0, min_len)};
    const std::vector<Interval> expected_intervals;
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, min_len);
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(IdentifyRegionsTest, identify_regions_1) {
    const std::vector<uint32_t> covgs({5, 6, 7, 5, 6, 6, 4, 4, 3, 5, 1, 1, 2, 3, 2, 4, 3});
    auto low_covg_intervals{identify_regions(covgs, 0, 0)};
    std::vector<Interval> expected_intervals;
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 0);
    expected_intervals = {Interval(10, 12)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 1);
    expected_intervals = {Interval(10, 12)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 2);
    expected_intervals = {Interval(10, 12)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 3);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(IdentifyRegionsTest, identify_regions_2) {
    std::vector<uint32_t> covgs({5, 6, 7, 5, 6, 6, 4, 4, 3, 5, 1, 1, 2, 3, 2, 4, 3});
    auto low_covg_intervals{identify_regions(covgs, 2, 0)};
    std::vector<Interval> expected_intervals{Interval(10, 13), Interval(14, 15)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 1);
    expected_intervals = {Interval(10, 13), Interval(14, 15)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 2);
    expected_intervals = {Interval(10, 13)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 3);
    expected_intervals = {Interval(10, 13)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 4);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(IdentifyRegionsTest, identify_regions_3) {
    const std::vector<uint32_t> covgs({5, 6, 7, 5, 6, 6, 4, 4, 3, 5, 1, 1, 2, 3, 2, 4, 3});
    auto low_covg_intervals{identify_regions(covgs, 3, 0)};
    std::vector<Interval> expected_intervals{Interval(8, 9), Interval(10, 15), Interval(16, 17)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 1);
    expected_intervals = {Interval(8, 9), Interval(10, 15), Interval(16, 17)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 2);
    expected_intervals = {Interval(10, 15)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 3);
    expected_intervals = {Interval(10, 15)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 4);
    expected_intervals = {Interval(10, 15)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 5);
    expected_intervals = {Interval(10, 15)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 6);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(IdentifyRegionsTest, identify_regions_4) {
    const std::vector<uint32_t> covgs({5, 6, 7, 5, 6, 6, 4, 4, 3, 5, 1, 1, 2, 3, 2, 4, 3});
    auto low_covg_intervals{identify_regions(covgs, 4, 0)};
    std::vector<Interval> expected_intervals{Interval(6, 9), Interval(10, 17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 1);
    expected_intervals = {Interval(6, 9), Interval(10, 17)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 2);
    expected_intervals = {Interval(6, 9), Interval(10, 17)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 3);
    expected_intervals = {Interval(6, 9), Interval(10, 17)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 4);
    expected_intervals = {Interval(10, 17)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 5);
    expected_intervals = {Interval(10, 17)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 6);
    expected_intervals = {Interval(10, 17)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 7);
    expected_intervals = {Interval(10, 17)};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 8);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(std::vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(FindIntervalInLocalPathTest, emptyInterval_returnEmptyResult) {
    const Interval interval; // fixme: need to implement an empty function on Interval
    LocalPRG local_prg{0, "test", "A"};
    const std::vector<LocalNodePtr> local_node_max_likelihood_path{local_prg.prg.nodes[0]};

    const auto actual {find_interval_in_localpath(interval, local_node_max_likelihood_path)};
    const PathComponents expected;

    EXPECT_EQ(actual, expected);
}

//TEST(FindIntervalInLocalPathTest, find_interval_in_localpath_minimal) {
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp{l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
//    // A G C T CGG  TAT
//    const auto found_components{find_interval_in_localpath(Interval(2, 3), lmp)};
//    const auto found_slice{l3.string_along_path(found_components.slice)};
//    const auto exp_slice{l3.prg.nodes[2]->seq};
//
//    EXPECT_EQ(exp_slice, found_slice);
//}
//
//TEST(ExtractReadsTest, find_interval_in_localpath_minimal_with_padding1) {
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp{l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
//    // A G C T CGG  TAT
//    const auto found_components{find_interval_in_localpath(Interval(2, 3), lmp)};
//    const auto found_slice{l3.string_along_path(found_components.slice)};
//    const auto exp_slice{LocalPRG::string_along_path({l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4]})};
//
//    EXPECT_EQ(exp_slice, found_slice);
//}
//
//TEST(ExtractReadsTest, find_interval_in_localpath_minimal_with_padding2) {
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp{l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
//    // A G C T CGG  TAT
//    const auto found_components{find_interval_in_localpath(Interval(2, 3), lmp)};
//    const auto found_slice{l3.string_along_path(found_components.slice)};
//    const auto exp_slice{LocalPRG::string_along_path({l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                     l3.prg.nodes[6]})};
//
//    EXPECT_EQ(exp_slice, found_slice);
//}
//
//TEST(ExtractReadsTest, find_interval_in_localpath_short) {
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp{l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
//    // A G C T CGG  TAT
//    const auto found_components{find_interval_in_localpath(Interval(1, 4), lmp)};
//    const auto found_slice{l3.string_along_path(found_components.slice)};
//    const auto exp_slice{LocalPRG::string_along_path({l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4]})};
//
//    EXPECT_EQ(exp_slice, found_slice);
//}
//
//TEST(ExtractReadsTest, find_interval_in_localpath_short_with_padding1) {
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp{l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
//    // A G C T CGG  TAT
//    const auto found_components{find_interval_in_localpath(Interval(1, 4), lmp)};
//    const auto found_slice{l3.string_along_path(found_components.slice)};
//    const auto exp_slice{LocalPRG::string_along_path({l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                     l3.prg.nodes[6]})};
//
//    EXPECT_EQ(exp_slice, found_slice);
//}
//
//TEST(ExtractReadsTest, find_interval_in_localpath_short_with_padding2) {
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp{l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
//    // A G C T CGG  TAT
//    const auto found_components{find_interval_in_localpath(Interval(1, 4), lmp)};
//    const auto found_slice{l3.string_along_path(found_components.slice)};
//    const auto exp_slice{LocalPRG::string_along_path({l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                     l3.prg.nodes[6]})};
//
//    EXPECT_EQ(exp_slice, found_slice);
//}
//
//TEST(ExtractReadsTest, find_interval_in_localpath_short_with_padding3) {
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp{l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
//    // A G C T CGG  TAT
//    const auto found_components{find_interval_in_localpath(Interval(1, 4), lmp)};
//    const auto found_slice{l3.string_along_path(found_components.slice)};
//    const auto exp_slice{LocalPRG::string_along_path({l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                     l3.prg.nodes[6]})};
//
//    EXPECT_EQ(exp_slice, found_slice);
//}
//
//TEST(ExtractReadsTest, find_interval_in_localpath_short_with_padding4) {
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp{l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
//    // A G C T CGG  TAT
//    const auto found_components{find_interval_in_localpath(Interval(1, 4), lmp)};
//    const auto found_slice{l3.string_along_path(found_components.slice)};
//    const auto exp_slice{LocalPRG::string_along_path({l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                     l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]})};
//
//    EXPECT_EQ(exp_slice, found_slice);
//}
//
//TEST(ExtractReadsTest, find_interval_in_localpath_multiple_sites) {
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp{l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
//                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
//    // A G C T CGG  TAT
//    const auto found_components{find_interval_in_localpath(Interval(2, 5), lmp)};
//    const auto found_slice{l3.string_along_path(found_components.slice)};
//    const auto exp_slice{LocalPRG::string_along_path({l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6]})};
//
//    EXPECT_EQ(exp_slice, found_slice);
//}

TEST(ExtractReadsTest, hits_inside_path) {
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    const std::vector<LocalNodePtr> lmp{//l3.prg.nodes[0],
            l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
            l3.prg.nodes[6], l3.prg.nodes[7]//, l3.prg.nodes[9]
    };
    // A G C T CGG  TAT
    std::set<MinimizerHitPtr, pComp_path> hits, expected_subset;
    uint32_t read_id = 0, prg_id = 3, knode_id = 0;
    const Interval read_interval(1, 4);
    const bool orientation(true);
    std::deque<Interval> d{Interval(7, 8), Interval(10, 12)};
    prg::Path prg_path;
    prg_path.initialize(d);

    // hit not on path
    MinimizerHitPtr mh(make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation));
    hits.insert(mh);

    // hits branching from path
    d = {Interval(7, 8), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31, 33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28, 30), Interval(33, 33), Interval(40, 41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    expected_subset.insert(mh);
    d = {Interval(8, 9), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    expected_subset.insert(mh);
    d = {Interval(16, 17), Interval(27, 29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    expected_subset.insert(mh);
    d = {Interval(27, 30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    expected_subset.insert(mh);

    // hit against different prg
    // hit different orientation
    prg::Path local_path;
    for (const auto &node : lmp) {
        local_path.add_end_interval(node->pos);
    }
    std::set<MinimizerHitPtr, pComp_path> found_subset{ find_hits_inside_path(hits, local_path)};
    EXPECT_EQ(expected_subset.size(), found_subset.size());
    auto jt = found_subset.begin();
    for (auto it = expected_subset.begin(); it != expected_subset.end() and jt != found_subset.end(); ++it) {
        EXPECT_EQ(**jt, **it);
        ++jt;
    }
}

TEST(ExtractReadsTest, get_read_overlap_coordinates) {
    //
    //  Read 0 has prg 3 sequence in interval (2,12] only
    //  Read 1 has prg 3 sequence in interval (6,16] as well as noise
    //  Read 2 has prg 3 sequence in interval (4,20] stretched out
    //  Read 3 has prg 3 sequence in interval (4,14] but is missing bits
    //  Read 4 doesn't have prg 3 sequence - on all hits are noise
    //
    uint32_t pnode_id = 3, prg_id = 3, read_id = 0, knode_id = 0;
    string pnode_name = "three";
    bool orientation(true);
    deque<Interval> d;
    prg::Path prg_path;
    MinimizerHitPtr mh;

    PanNodePtr pn = make_shared<pangenome::Node>(pnode_id, prg_id, pnode_name);
    PanReadPtr pr = make_shared<pangenome::Read>(read_id);

    set<MinimizerHitPtr, pComp> hits;


    // READ 0
    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(2, 5), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28, 30), Interval(33, 33), Interval(40, 41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(3, 6), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(8, 9), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16, 17), Interval(27, 29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(27, 30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 1
    read_id = 1;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(12, 15), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28, 30), Interval(33, 33), Interval(40, 41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(11, 14), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(8, 9), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16, 17), Interval(27, 29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(27, 30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = {Interval(7, 8), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31, 33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(78, 81)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 2
    read_id = 2;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(17, 20), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28, 30), Interval(33, 33), Interval(40, 41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(15, 18), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(8, 9), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16, 17), Interval(27, 29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(27, 30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = {Interval(7, 8), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31, 33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(78, 81)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 3
    read_id = 3;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28, 30), Interval(33, 33), Interval(40, 41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(8, 9), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16, 17), Interval(27, 29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);


    // noise
    d = {Interval(7, 8), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 4
    read_id = 4;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(17, 20), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = {Interval(7, 8), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31, 33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(78, 81)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // RUN GET_READ_OVERLAPS
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    const std::vector<LocalNodePtr> lmp{//l3.prg.nodes[0],
            l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
            l3.prg.nodes[6], l3.prg.nodes[7]//, l3.prg.nodes[9]
    };
    // A G C T CGG  TAT
    const std::set<ReadCoordinate> expected_overlaps{{0, 3, 9,  1},
                                                  {1, 7, 13, 1},
                                                  {2, 5, 13, 1},
                                                  {3, 6, 10, 1}};

    prg::Path local_path;
    for (const auto &node : lmp) {
    local_path.add_end_interval(node->pos);
    }
    const auto overlaps{get_read_overlap_coordinates(pn, local_path)};

    EXPECT_ITERABLE_EQ(std::set<ReadCoordinate>, expected_overlaps, overlaps);
}


TEST(ExtractReadsTest, get_read_overlap_coordinates_no_duplicates) {
    //
    //  Read 0 has prg 3 sequence in interval (2,12] only
    //  Read 1 has prg 3 sequence in interval (6,16] as well as noise
    //  Read 2 has prg 3 sequence in interval (4,20] stretched out
    //  Read 3 has prg 3 sequence in interval (4,14] but is missing bits
    //  Read 4 doesn't have prg 3 sequence - on all hits are noise
    //  Read 5 is a duplicate of Read 0
    //
    uint32_t pnode_id = 3, prg_id = 3, read_id = 0, knode_id = 0;
    string pnode_name = "three";
    bool orientation(true);
    deque<Interval> d;
    prg::Path prg_path;
    MinimizerHitPtr mh;

    PanNodePtr pn = make_shared<pangenome::Node>(pnode_id, prg_id, pnode_name);
    PanReadPtr pr = make_shared<pangenome::Read>(read_id);

    set<MinimizerHitPtr, pComp> hits;


    // READ 0
    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(2, 5), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28, 30), Interval(33, 33), Interval(40, 41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(3, 6), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(8, 9), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16, 17), Interval(27, 29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(27, 30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 1
    read_id = 1;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(12, 15), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28, 30), Interval(33, 33), Interval(40, 41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(11, 14), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(8, 9), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16, 17), Interval(27, 29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(27, 30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = {Interval(7, 8), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31, 33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(78, 81)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 2
    read_id = 2;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(17, 20), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28, 30), Interval(33, 33), Interval(40, 41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(15, 18), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(8, 9), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16, 17), Interval(27, 29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(27, 30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = {Interval(7, 8), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31, 33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(78, 81)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 3
    read_id = 3;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28, 30), Interval(33, 33), Interval(40, 41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(8, 9), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16, 17), Interval(27, 29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);


    // noise
    d = {Interval(7, 8), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 4
    read_id = 4;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(17, 20), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = {Interval(7, 8), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31, 33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(78, 81)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 5
    read_id = 0;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(2, 5), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(33, 33), Interval(40, 42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28, 30), Interval(33, 33), Interval(40, 41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(3, 6), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(8, 9), Interval(16, 17), Interval(27, 28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16, 17), Interval(27, 29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(27, 30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // RUN GET_READ_OVERLAPS
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    const std::vector<LocalNodePtr> lmp{//l3.prg.nodes[0],
            l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
            l3.prg.nodes[6], l3.prg.nodes[7]//, l3.prg.nodes[9]
    };
    // A G C T CGG  TAT
    const std::set<ReadCoordinate> expected_overlaps{{0, 3, 9,  1},
                                                  {1, 7, 13, 1},
                                                  {2, 5, 13, 1},
                                                  {3, 6, 10, 1}};

    prg::Path local_path;
    for (const auto &node : lmp) {
    local_path.add_end_interval(node->pos);
    }
    const auto overlaps{get_read_overlap_coordinates(pn, local_path)};

    EXPECT_ITERABLE_EQ(std::set<ReadCoordinate>, expected_overlaps, overlaps);
}


TEST(ExtractReadsTest, add_pnode_coordinate_pairs) {
    //
    //  Read 0 has both intervals
    //  Read 1 has first interval
    //  Read 2 has second interval
    //  Read 3 has both intervals, but only sparse hits which do not meet threshold
    //

    uint32_t pnode_id = 3, prg_id = 3, read_id = 0, knode_id = 0;
    string pnode_name = "three";
    bool orientation(true);
    MinimizerHitPtr mh;

    // define localPRG
    LocalPRG l3(prg_id, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    auto index = std::make_shared<Index>();
    auto w = 1, k = 3;
    l3.minimizer_sketch(index, w, k);

    // define localpath and kmerpath
    // corresponds to sequence (A) G C T CGG  (TAT)
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0],
                                l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]
    };

    vector<KmerNodePtr> kmp = {l3.kmer_prg.nodes[0],
                               l3.kmer_prg.nodes[1], l3.kmer_prg.nodes[4], l3.kmer_prg.nodes[8],
                               l3.kmer_prg.nodes[13], l3.kmer_prg.nodes[15], l3.kmer_prg.nodes[17],
                               l3.kmer_prg.nodes[19], l3.kmer_prg.nodes[20], l3.kmer_prg.nodes[21]
    };

    PanNodePtr pn = make_shared<pangenome::Node>(pnode_id, prg_id, pnode_name);
    pn->kmer_prg = l3.kmer_prg;

    bool strand = 0;
    uint32_t sample_id = 0;

    pn->kmer_prg.setup_coverages(1);

    pn->kmer_prg.nodes[0]->set_covg(10, strand, sample_id);
    pn->kmer_prg.nodes[1]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[4]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[8]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[13]->set_covg(10, strand, sample_id);
    pn->kmer_prg.nodes[15]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[17]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[19]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[20]->set_covg(10, strand, sample_id);
    pn->kmer_prg.nodes[21]->set_covg(10, strand, sample_id);

    PanReadPtr pr = make_shared<pangenome::Read>(read_id);
    set<MinimizerHitPtr, pComp> hits;

    // READ 0
    // covers whole path in interval 2,12
    mh = make_shared<MinimizerHit>(read_id, Interval(2, 5), prg_id, pn->kmer_prg.nodes[1]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(3, 6), prg_id, pn->kmer_prg.nodes[4]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, pn->kmer_prg.nodes[8]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, pn->kmer_prg.nodes[13]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, pn->kmer_prg.nodes[15]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, pn->kmer_prg.nodes[17]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, pn->kmer_prg.nodes[19]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, pn->kmer_prg.nodes[20]->path, knode_id,
                                   orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 1
    //covers first low covg site in interval [6,10]
    read_id = 1;
    pr = make_shared<pangenome::Read>(read_id);

    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, pn->kmer_prg.nodes[1]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, pn->kmer_prg.nodes[4]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, pn->kmer_prg.nodes[8]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, pn->kmer_prg.nodes[13]->path, knode_id,
                                   orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 2
    // covers second site in interval [16,22]
    read_id = 2;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    mh = make_shared<MinimizerHit>(read_id, Interval(15, 18), prg_id, pn->kmer_prg.nodes[13]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(16, 19), prg_id, pn->kmer_prg.nodes[15]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(17, 20), prg_id, pn->kmer_prg.nodes[17]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(18, 21), prg_id, pn->kmer_prg.nodes[19]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(19, 22), prg_id, pn->kmer_prg.nodes[20]->path, knode_id,
                                   orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 3
    // has sparse hits
    read_id = 3;
    pr = make_shared<pangenome::Read>(read_id);

    mh = make_shared<MinimizerHit>(read_id, Interval(3, 6), prg_id, pn->kmer_prg.nodes[4]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, pn->kmer_prg.nodes[8]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, pn->kmer_prg.nodes[19]->path, knode_id,
                                   orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    auto buff = 1, covg_thresh = 1, min_length = 1, min_num_hits = 2;

    std::vector<std::shared_ptr<LocalPRG>> prgs;
    prgs.emplace_back(std::make_shared<LocalPRG>(l3));
    std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> pairs;
    denovo_discovery::add_pnode_coordinate_pairs(prgs, pairs, pn, lmp, kmp, buff, covg_thresh, min_length);

    GeneIntervalInfo interval_info1{pn, Interval(1, 3), "AGCT"};
    GeneIntervalInfo interval_info2{pn, Interval(6, 7), "CGGTAT"};
    ReadCoordinate read_coord1{0, 2, 6, true};
    ReadCoordinate read_coord2{0, 6, 12, true};
    ReadCoordinate read_coord3{1, 6, 10, true};
    ReadCoordinate read_coord4{2, 16, 22, true};

    std::vector<std::pair<ReadCoordinate, GeneIntervalInfo>> expected_coords = {
            std::make_pair(read_coord1, interval_info1),
            std::make_pair(read_coord2, interval_info2),
            std::make_pair(read_coord3, interval_info1),
            std::make_pair(read_coord4, interval_info2)};
    EXPECT_EQ(expected_coords.size(), pairs.size());
    uint count = 0;
    for (const auto &p : pairs) {
        if (expected_coords.size() <= count)
            break;
        EXPECT_EQ(p.first, expected_coords[count].first);
        EXPECT_EQ(p.second, expected_coords[count].second);
        count++;
    }
}

TEST(ExtractReadsTest, add_pnode_coordinate_pairs_fewer_hits_needed) {
    //
    //  Read 0 has first interval
    //  Read 1 has second interval
    //  Read 2 has both intervals, but only sparse hits (meets lower threshold)
    //  Read 3 has both intervals
    //

    uint32_t pnode_id = 3, prg_id = 3, read_id = 0, knode_id = 0;
    string pnode_name = "three";
    bool orientation(true);
    MinimizerHitPtr mh;

    // define localPRG
    LocalPRG l3(prg_id, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    auto index = std::make_shared<Index>();
    auto w = 1, k = 3;
    l3.minimizer_sketch(index, w, k);

    // define localpath and kmerpath
    // corresponds to sequence (A) G C T CGG  (TAT)
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0],
                                l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]
    };

    vector<KmerNodePtr> kmp = {l3.kmer_prg.nodes[0],
                               l3.kmer_prg.nodes[1], l3.kmer_prg.nodes[4], l3.kmer_prg.nodes[8],
                               l3.kmer_prg.nodes[13], l3.kmer_prg.nodes[15], l3.kmer_prg.nodes[17],
                               l3.kmer_prg.nodes[19], l3.kmer_prg.nodes[20], l3.kmer_prg.nodes[21]
    };

    PanNodePtr pn = make_shared<pangenome::Node>(pnode_id, prg_id, pnode_name);
    pn->kmer_prg = l3.kmer_prg;

    uint32_t strand = 0;
    uint32_t sample_id = 0;
    pn->kmer_prg.setup_coverages(1);

    pn->kmer_prg.nodes[0]->set_covg(10, strand, sample_id);
    pn->kmer_prg.nodes[1]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[4]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[8]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[13]->set_covg(10, strand, sample_id);
    pn->kmer_prg.nodes[15]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[17]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[19]->set_covg(1, strand, sample_id);
    pn->kmer_prg.nodes[20]->set_covg(10, strand, sample_id);
    pn->kmer_prg.nodes[21]->set_covg(10, strand, sample_id);

    set<MinimizerHitPtr, pComp> hits;

    // READ 0
    //covers first low covg site in interval [6,10]
    read_id = 0;
    PanReadPtr pr = make_shared<pangenome::Read>(read_id);

    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, pn->kmer_prg.nodes[1]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, pn->kmer_prg.nodes[4]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, pn->kmer_prg.nodes[8]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, pn->kmer_prg.nodes[13]->path, knode_id,
                                   orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 1
    // covers second site in interval [16,22]
    read_id = 1;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    mh = make_shared<MinimizerHit>(read_id, Interval(15, 18), prg_id, pn->kmer_prg.nodes[13]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(16, 19), prg_id, pn->kmer_prg.nodes[15]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(17, 20), prg_id, pn->kmer_prg.nodes[17]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(18, 21), prg_id, pn->kmer_prg.nodes[19]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(19, 22), prg_id, pn->kmer_prg.nodes[20]->path, knode_id,
                                   orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 2
    // has sparse hits
    read_id = 2;
    pr = make_shared<pangenome::Read>(read_id);

    mh = make_shared<MinimizerHit>(read_id, Interval(3, 6), prg_id, pn->kmer_prg.nodes[4]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, pn->kmer_prg.nodes[8]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, pn->kmer_prg.nodes[19]->path, knode_id,
                                   orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 3
    // covers whole path in interval 2,12
    read_id = 3;
    pr = make_shared<pangenome::Read>(read_id);

    mh = make_shared<MinimizerHit>(read_id, Interval(2, 5), prg_id, pn->kmer_prg.nodes[1]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(3, 6), prg_id, pn->kmer_prg.nodes[4]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, pn->kmer_prg.nodes[8]->path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, pn->kmer_prg.nodes[13]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, pn->kmer_prg.nodes[15]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, pn->kmer_prg.nodes[17]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, pn->kmer_prg.nodes[19]->path, knode_id,
                                   orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, pn->kmer_prg.nodes[20]->path, knode_id,
                                   orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pn->reads.insert(pr);
    hits.clear();

    auto buff = 1, covg_thresh = 1, min_length = 1, min_num_hits = 1;

    std::vector<std::shared_ptr<LocalPRG>> prgs;
    prgs.emplace_back(std::make_shared<LocalPRG>(l3));
    std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> pairs;
    denovo_discovery::add_pnode_coordinate_pairs(prgs, pairs, pn, lmp, kmp, buff, covg_thresh, min_length);

    GeneIntervalInfo interval_info1{pn, Interval(1, 3), "AGCT"};
    GeneIntervalInfo interval_info2{pn, Interval(6, 7), "CGGTAT"};

    ReadCoordinate read_coord1{0, 6, 10, true};
    ReadCoordinate read_coord2{1, 16, 22, true};
    ReadCoordinate read_coord3{2, 3, 6, true};
    ReadCoordinate read_coord4{2, 8, 11, true};
    ReadCoordinate read_coord5{3, 2, 6, true};
    ReadCoordinate read_coord6{3, 6, 12, true};

    std::vector<std::pair<ReadCoordinate, GeneIntervalInfo>> expected_coords = {
            std::make_pair(read_coord1, interval_info1),
            std::make_pair(read_coord2, interval_info2),
            std::make_pair(read_coord3, interval_info1),
            std::make_pair(read_coord4, interval_info2),
            std::make_pair(read_coord5, interval_info1),
            std::make_pair(read_coord6, interval_info2)};
    EXPECT_EQ(expected_coords.size(), pairs.size());
    uint count = 0;
    for (const auto &p : pairs) {
        if (expected_coords.size() <= count)
            break;
        EXPECT_EQ(p.first, expected_coords[count].first);
        EXPECT_EQ(p.second, expected_coords[count].second);
        count++;
    }
}

TEST(ExtractReadsTest, read_coordinate_ordering) {
    ReadCoordinate read_coord1{0, 6, 10, true};
    ReadCoordinate read_coord2{1, 16, 22, true};
    ReadCoordinate read_coord3{2, 3, 6, true};
    ReadCoordinate read_coord4{2, 8, 11, true};
    ReadCoordinate read_coord5{3, 2, 6, true};
    ReadCoordinate read_coord6{3, 6, 12, true};
    ReadCoordinate read_coord7{3, 2, 6, false};
    ReadCoordinate read_coord8{3, 6, 12, false};

    set<ReadCoordinate> read_coords = {read_coord8, read_coord7, read_coord6, read_coord5, read_coord4,
                                       read_coord3, read_coord2, read_coord1};
    vector<ReadCoordinate> exp_read_coords = {read_coord1, read_coord2, read_coord3, read_coord4, read_coord5,
                                              read_coord7, read_coord6, read_coord8};

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
    ReadCoordinate read_coord1{0, 6, 10, true};
    ReadCoordinate read_coord2{1, 16, 22, true};
    ReadCoordinate read_coord3{2, 3, 6, true};
    ReadCoordinate read_coord4{2, 8, 11, true};
    ReadCoordinate read_coord5{3, 2, 6, true};
    ReadCoordinate read_coord6{3, 6, 12, true};
    ReadCoordinate read_coord7{3, 2, 6, false};
    ReadCoordinate read_coord8{3, 6, 12, false};

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


TEST(ExtractReadsTest, ordering_of_coordinate_pairs) {
    uint32_t pnode_id = 3, prg_id = 3, read_id = 0;
    string pnode_name = "three";

    PanNodePtr pn = make_shared<pangenome::Node>(pnode_id, prg_id, pnode_name);

    GeneIntervalInfo interval_info1{pn, Interval(1, 3), "AGCT"};
    GeneIntervalInfo interval_info2{pn, Interval(6, 7), "CGGTAT"};
    GeneIntervalInfo interval_info3{pn, Interval(6, 9), "CGGTAT"};

    ReadCoordinate read_coord1{0, 6, 10, true};
    ReadCoordinate read_coord2{1, 16, 22, true};
    ReadCoordinate read_coord3{2, 3, 6, true};
    ReadCoordinate read_coord4{2, 8, 11, true};
    ReadCoordinate read_coord5{3, 2, 6, true};
    ReadCoordinate read_coord6{3, 6, 12, true};
    ReadCoordinate read_coord7{3, 2, 6, false};
    ReadCoordinate read_coord8{3, 6, 12, false};

    std::vector<std::pair<ReadCoordinate, GeneIntervalInfo>> insert_order_pairs = {
            std::make_pair(read_coord1, interval_info1),
            std::make_pair(read_coord8, interval_info1),
            std::make_pair(read_coord4, interval_info1),
            std::make_pair(read_coord5, interval_info1),
            std::make_pair(read_coord7, interval_info1),
            std::make_pair(read_coord6, interval_info1),
            std::make_pair(read_coord2, interval_info1),
            std::make_pair(read_coord3, interval_info1),
            std::make_pair(read_coord7, interval_info3),
            std::make_pair(read_coord6, interval_info3),
            std::make_pair(read_coord7, interval_info2),
            std::make_pair(read_coord7, interval_info2),
            std::make_pair(read_coord6, interval_info2),
            std::make_pair(read_coord7, interval_info2),};

    std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> pairs(insert_order_pairs.begin(), insert_order_pairs.end());

    vector<uint> exp_order = {0, 6, 7, 2, 3, 4, 10, 8, 5, 12, 9, 1};
    EXPECT_EQ(exp_order.size(), pairs.size());
    uint count = 0;
    for (const auto &p : pairs) {
        if (exp_order.size() <= count)
            break;
        EXPECT_EQ(p.first, insert_order_pairs[exp_order[count]].first);
        EXPECT_EQ(p.second, insert_order_pairs[exp_order[count]].second);
        count++;
    }
}

TEST(ExtractReadsTest, collect_read_pileups) {
    uint32_t pnode_id = 3, prg_id = 3, read_id = 0;
    string pnode_name = "three";

    PanNodePtr pn = make_shared<pangenome::Node>(pnode_id, prg_id, pnode_name);

    GeneIntervalInfo interval_info1{pn, Interval(1, 3), "AGCT"};
    GeneIntervalInfo interval_info2{pn, Interval(6, 7), "CGGTAT"};
    GeneIntervalInfo interval_info3{pn, Interval(6, 9), "CGGTAT"};

    ReadCoordinate read_coord1{0, 6, 10, true}; //GCTA
    ReadCoordinate read_coord2{1, 16, 22, true}; //CGGTAT
    ReadCoordinate read_coord3{2, 3, 6, true}; //GTT
    ReadCoordinate read_coord4{2, 8, 11, true}; //TGT
    ReadCoordinate read_coord5{3, 2, 6, true}; //TCAC
    ReadCoordinate read_coord6{3, 6, 12, true}; //CTTAAC
    ReadCoordinate read_coord7{3, 2, 6, false}; //GTGA
    ReadCoordinate read_coord8{3, 6, 12, false}; //GTTAAG

    std::vector<std::pair<ReadCoordinate, GeneIntervalInfo>> insert_order_pairs = {
            std::make_pair(read_coord1, interval_info1),
            std::make_pair(read_coord8, interval_info1),
            std::make_pair(read_coord4, interval_info1),
            std::make_pair(read_coord5, interval_info1),
            std::make_pair(read_coord7, interval_info1),
            std::make_pair(read_coord6, interval_info1),
            std::make_pair(read_coord2, interval_info1),
            std::make_pair(read_coord3, interval_info1),
            std::make_pair(read_coord7, interval_info3),
            std::make_pair(read_coord6, interval_info3),
            std::make_pair(read_coord7, interval_info2),
            std::make_pair(read_coord7, interval_info2),
            std::make_pair(read_coord6, interval_info2),
            std::make_pair(read_coord7, interval_info2),};

    std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> pairs(insert_order_pairs.begin(), insert_order_pairs.end());

    boost::filesystem::path fpath("../../test/test_cases/reads_for_pileups.fa");
    auto padding_size = 0;
    auto pileups = denovo_discovery::collect_read_pileups(pairs, fpath);

    std::map<GeneIntervalInfo, ReadPileup> exp_pileups;
    exp_pileups[interval_info1] = {"GCTA", "CGGTAT", "GTT", "TGT", "TCAC", "GTGA", "CTTAAC", "GTTAAG"};
    exp_pileups[interval_info2] = {"GTGA", "CTTAAC"};
    exp_pileups[interval_info3] = {"GTGA", "CTTAAC"};

    EXPECT_ITERABLE_EQ(ReadPileup, exp_pileups[interval_info1], pileups[interval_info1]);
    EXPECT_ITERABLE_EQ(ReadPileup, exp_pileups[interval_info2], pileups[interval_info2]);
    EXPECT_ITERABLE_EQ(ReadPileup, exp_pileups[interval_info3], pileups[interval_info3]);

}
