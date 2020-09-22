#include "gtest/gtest.h"
#include <cstdint>
#include <iostream>
#include "interval.h"

TEST(IntervalTest, create)
{
    Interval i(0, 0);
    uint32_t j = 0;
    EXPECT_EQ(i.start, j);
    EXPECT_EQ(i.get_end(), j);
    EXPECT_EQ(i.length, j);

    i = Interval(1, 9);
    j = 1;
    EXPECT_EQ(i.start, j);
    j = 9;
    EXPECT_EQ(i.get_end(), j);
    j = 8;
    EXPECT_EQ(i.length, j);

    // should fail if end is before start
    EXPECT_DEATH(Interval(9, 1), "");
    // input should be non-negative
    EXPECT_DEATH(Interval(-1, 10), "");
}

TEST(IntervalTest, write)
{
    Interval i(1, 5);
    std::stringstream out;
    out << i;
    EXPECT_EQ(out.str(), "[1, 5)");
}

TEST(IntervalTest, read)
{
    Interval i(1, 5);
    std::stringstream out;
    out << i;
    Interval j;
    out >> j;
    EXPECT_EQ(i, j);
}

TEST(IntervalTest, equals)
{
    Interval i(1, 5);
    Interval j(1, 5);
    EXPECT_EQ(i, j);
    EXPECT_EQ(j, i);

    Interval k(0, 4);
    EXPECT_EQ((i == k), false);
    EXPECT_EQ((k == j), false);

    i = Interval(0, 0);
    j = Interval(1, 1);
    EXPECT_EQ((i == i), true);
    EXPECT_EQ((j == j), true);
    EXPECT_EQ((i == j), false);
    EXPECT_EQ(i, i);
    EXPECT_EQ(j, j);
}

TEST(IntervalTest, notequals)
{
    Interval i(1, 5);
    Interval j(1, 5);
    EXPECT_EQ((i != j), false);

    Interval k(0, 4);
    EXPECT_EQ((i != k), true);
    EXPECT_EQ((j != k), true);
    EXPECT_NE(i, k);
    EXPECT_NE(k, j);

    i = Interval(0, 0);
    j = Interval(1, 1);
    EXPECT_EQ((i != i), false);
    EXPECT_EQ((j != j), false);
    EXPECT_EQ((i != j), true);
    EXPECT_NE(i, j);
    EXPECT_NE(j, i);
}

TEST(IntervalTest, lessthan)
{
    Interval i(1, 5);
    Interval j(2, 5);
    Interval k(0, 4);
    Interval l(0, 7);

    EXPECT_EQ((i < i), false);

    EXPECT_EQ((i < j), true);
    EXPECT_EQ((j < i), false);

    EXPECT_EQ((i < k), false);
    EXPECT_EQ((k < i), true);
    EXPECT_EQ((k < j), true);
    EXPECT_EQ((j < k), false);

    EXPECT_EQ((k < l), true);
    EXPECT_EQ((l < i), true);
    EXPECT_EQ((l < j), true);
    EXPECT_EQ((l < k), false);
    EXPECT_EQ((i < l), false);
    EXPECT_EQ((j < l), false);
}

TEST(intervalEmptyTest, emptyIntervalReturnsTrue)
{
    const Interval empty_interval {};

    EXPECT_TRUE(empty_interval.empty());
}

TEST(intervalEmptyTest, nonEmptyIntervalReturnsFalse)
{
    const Interval non_empty_interval { 1, 4 };

    EXPECT_FALSE(non_empty_interval.empty());
}

class MergeIntervalsTest : public ::testing::Test {
protected:
    void SetUp() override
    {
        single_.insert(single_.end(), { Interval(0, 1) });
        unsorted_disjoint_.insert(
            unsorted_disjoint_.end(), { Interval(4, 6), Interval(0, 2) });
        contained_.insert(contained_.end(), { Interval(0, 10), Interval(4, 6) });
        overlap_.insert(overlap_.end(), { Interval(0, 3), Interval(2, 5) });
        overlap_and_disjoint_.insert(overlap_and_disjoint_.end(),
            { Interval(0, 3), Interval(2, 5), Interval(7, 10) });
        disjoint_overlap_disjoint_.insert(disjoint_overlap_disjoint_.end(),
            { Interval(0, 3), Interval(6, 9), Interval(8, 11), Interval(16, 20),
                Interval(30, 34) });
    }
    std::vector<Interval> empty_;
    std::vector<Interval> single_;
    std::vector<Interval> unsorted_disjoint_;
    std::vector<Interval> contained_;
    std::vector<Interval> overlap_;
    std::vector<Interval> overlap_and_disjoint_;
    std::vector<Interval> disjoint_overlap_disjoint_;
};

TEST_F(MergeIntervalsTest, HandlesEmptyInput)
{
    const auto dist { 0 };

    merge_intervals_within(empty_, dist);
    const auto actual { empty_ };
    const std::vector<Interval> expected;

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, SingleElementDoesNothing)
{
    const auto dist { 0 };

    const auto expected = single_;
    merge_intervals_within(single_, dist);
    const auto actual = single_;

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, UnsortedDisjointIntervalsReturnsSortedIntervals)
{
    const auto dist { 0 };

    std::vector<Interval> expected = unsorted_disjoint_;
    std::sort(expected.begin(), expected.end());
    merge_intervals_within(unsorted_disjoint_, dist);
    const auto actual = unsorted_disjoint_;

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, DisjointWithDistEdgeCase)
{
    const auto dist { 2 };

    std::vector<Interval> expected = unsorted_disjoint_;
    std::sort(expected.begin(), expected.end());
    merge_intervals_within(unsorted_disjoint_, dist);
    const auto actual = unsorted_disjoint_;

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, DisjointWithDistCausingOverlap)
{
    const auto dist { 3 };

    merge_intervals_within(unsorted_disjoint_, dist);
    const auto actual = unsorted_disjoint_;
    const std::vector<Interval> expected = { Interval(0, 6) };

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, ContainedIntervalIsMerged)
{
    const auto dist { 0 };

    merge_intervals_within(contained_, dist);
    const auto actual = contained_;
    const std::vector<Interval> expected = { Interval(0, 10) };

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, OverlappingIsMerged)
{
    const auto dist { 0 };

    merge_intervals_within(overlap_, dist);
    const auto actual = overlap_;
    const std::vector<Interval> expected = { Interval(0, 5) };

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, OverlapAndDisjointNoDistOnlyOverlapMerged)
{
    const auto dist { 0 };

    merge_intervals_within(overlap_and_disjoint_, dist);
    const auto actual = overlap_and_disjoint_;
    const std::vector<Interval> expected = { Interval(0, 5), Interval(7, 10) };

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, OverlapAndDisjointEdgeDistOnlyOverlapMerged)
{
    const auto dist { 2 };

    merge_intervals_within(overlap_and_disjoint_, dist);
    const auto actual = overlap_and_disjoint_;
    const std::vector<Interval> expected = { Interval(0, 5), Interval(7, 10) };

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, OverlapAndDisjointOverlapDistMergeAll)
{
    const auto dist { 3 };

    merge_intervals_within(overlap_and_disjoint_, dist);
    const auto actual = overlap_and_disjoint_;
    const std::vector<Interval> expected = { Interval(0, 10) };

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, ManyDisjointAndOverlappingNoDistOnlyMergeOverlaps)
{
    const auto dist { 0 };

    merge_intervals_within(disjoint_overlap_disjoint_, dist);
    const auto actual = disjoint_overlap_disjoint_;
    const std::vector<Interval> expected
        = { Interval(0, 3), Interval(6, 11), Interval(16, 20), Interval(30, 34) };

    EXPECT_EQ(actual, expected);
}

TEST_F(
    MergeIntervalsTest, ManyDisjointAndOverlappingSmallDistMergeOverlapsAndOneDisjoint)
{
    const auto dist { 5 };

    merge_intervals_within(disjoint_overlap_disjoint_, dist);
    const auto actual = disjoint_overlap_disjoint_;
    const std::vector<Interval> expected
        = { Interval(0, 11), Interval(16, 20), Interval(30, 34) };

    EXPECT_EQ(actual, expected);
}

TEST_F(
    MergeIntervalsTest, ManyDisjointAndOverlappingMediumDistMergeOverlapsAndTwoDisjoint)
{
    const auto dist { 9 };

    merge_intervals_within(disjoint_overlap_disjoint_, dist);
    const auto actual = disjoint_overlap_disjoint_;
    const std::vector<Interval> expected = { Interval(0, 20), Interval(30, 34) };

    EXPECT_EQ(actual, expected);
}

TEST_F(MergeIntervalsTest, ManyDisjointAndOverlappingLargeDistMergeAll)
{
    const auto dist { 100 };

    merge_intervals_within(disjoint_overlap_disjoint_, dist);
    const auto actual = disjoint_overlap_disjoint_;
    const std::vector<Interval> expected = { Interval(0, 34) };

    EXPECT_EQ(actual, expected);
}
