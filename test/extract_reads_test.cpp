#include <iostream>
#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "extract_reads.h"
#include "interval.h"
#include "localPRG.h"

using namespace std;

TEST(ExtractReadsTest, identify_regions_null) {
    vector<uint32_t> covgs;
    auto low_covg_intervals = identify_regions(covgs, 0);
    vector<Interval> expected_intervals;
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2);
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(ExtractReadsTest, identify_regions_1) {
    vector<uint32_t> covgs({5, 6, 7, 5, 6, 6, 4, 4, 3, 5, 1, 1, 2, 3, 2, 4, 3});
    auto low_covg_intervals = identify_regions(covgs, 0, 0);
    vector<Interval> expected_intervals;
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 0);
    expected_intervals = {Interval(10, 12)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 1);
    expected_intervals = {Interval(10, 12)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 2);
    expected_intervals = {Interval(10, 12)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 3);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(ExtractReadsTest, identify_regions_2) {
    vector<uint32_t> covgs({5,6,7,5,6,6,4,4,3,5,1,1,2,3,2,4,3});
    auto low_covg_intervals = identify_regions(covgs, 2, 0);
    vector<Interval> expected_intervals = {Interval(10,13),Interval(14,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 1);
    expected_intervals = {Interval(10,13),Interval(14,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 2);
    expected_intervals = {Interval(10,13)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 3);
    expected_intervals = {Interval(10,13)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 4);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(ExtractReadsTest, identify_regions_3) {
    vector<uint32_t> covgs({5,6,7,5,6,6,4,4,3,5,1,1,2,3,2,4,3});
    auto low_covg_intervals = identify_regions(covgs, 3, 0);
    vector<Interval> expected_intervals = {Interval(8,9),Interval(10,15), Interval(16,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 1);
    expected_intervals = {Interval(8,9),Interval(10,15), Interval(16,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 2);
    expected_intervals = {Interval(10,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 3);
    expected_intervals = {Interval(10,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 4);
    expected_intervals = {Interval(10,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 5);
    expected_intervals = {Interval(10,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 6);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(ExtractReadsTest, identify_regions_4) {
    vector<uint32_t> covgs({5,6,7,5,6,6,4,4,3,5,1,1,2,3,2,4,3});
    auto low_covg_intervals = identify_regions(covgs, 4, 0);
    vector<Interval> expected_intervals = {Interval(6,9),Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 1);
    expected_intervals = {Interval(6,9),Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 2);
    expected_intervals = {Interval(6,9),Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 3);
    expected_intervals = {Interval(6,9),Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 4);
    expected_intervals = {Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 5);
    expected_intervals = {Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 6);
    expected_intervals = {Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 7);
    expected_intervals = {Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 8);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(ExtractReadsTest, find_interval_in_localpath_minimal) {
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
    // A G C T CGG  TAT
    auto found_path = find_interval_in_localpath(Interval(2,3), lmp);

    vector<LocalNodePtr> exp_path = {l3.prg.nodes[2]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, exp_path, found_path);
}

TEST(ExtractReadsTest, find_interval_in_localpath_short) {
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
    // A G C T CGG  TAT
    auto found_path = find_interval_in_localpath(Interval(1,4), lmp);
    for (auto n : found_path){
        cout << n->pos << " ";
    }
    cout << endl;

    vector<LocalNodePtr> exp_path = {l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, exp_path, found_path);
}

TEST(ExtractReadsTest, find_interval_in_localpath_uneven) {
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
    // A G C T CGG  TAT
    auto found_path = find_interval_in_localpath(Interval(2,4), lmp);

    vector<LocalNodePtr> exp_path = {l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, exp_path, found_path);

    found_path = find_interval_in_localpath(Interval(1,3), lmp);

    exp_path = {l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, exp_path, found_path);
}

TEST(ExtractReadsTest, find_interval_in_localpath_multiple_sites) {
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
    // A G C T CGG  TAT
    auto found_path = find_interval_in_localpath(Interval(2,5), lmp);

    vector<LocalNodePtr> exp_path = {l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                     l3.prg.nodes[6], l3.prg.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, exp_path, found_path);

}
