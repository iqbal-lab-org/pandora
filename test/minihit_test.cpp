#include "gtest/gtest.h"
#include "minihit.h"
#include "minimizer.h"
#include "minirecord.h"
#include "interval.h"
#include "prg/path.h"
#include "inthash.h"
#include <stdint.h>
#include <iostream>
#include <algorithm>

using namespace std;

TEST(MinimizerHitTest, create)
{
    KmerHash hash;
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    Minimizer m(min(kh.first, kh.second), 0, 5, 0);
    deque<Interval> d = { Interval(7, 8), Interval(10, 14) };
    prg::Path p;
    p.initialize(d);
    MiniRecord mr(0, p, 0, 0);
    MinimizerHit mh(1, m, mr);
    uint32_t j(1);
    EXPECT_EQ(j, mh.get_read_id());
    EXPECT_EQ((uint)0, mh.get_read_start_position());
    j = 0;
    EXPECT_EQ(j, mh.get_prg_id());
    EXPECT_EQ(p, mh.get_prg_path());
    bool b = true;
    EXPECT_EQ(b, mh.is_forward());

    kh = hash.kmerhash("hell", 4);
    m = Minimizer(min(kh.first, kh.second), 1, 5, 0);
    EXPECT_DEATH(MinimizerHit(1, m, mr), "");
    // TEST SECOND CONSTRUCTOR!!
}

TEST(MinimizerHitTest, checkStrand)
{
    KmerHash hash;
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    deque<Interval> d = { Interval(7, 8), Interval(10, 14) };
    prg::Path p;
    p.initialize(d);

    {
        Minimizer m(min(kh.first, kh.second), 0, 5, 0);
        ;
        MiniRecord mr(0, p, 0, 0);
        MinimizerHit mh(1, m, mr);
        EXPECT_EQ(mh.is_forward(), true);
    }

    {
        Minimizer m(min(kh.first, kh.second), 0, 5, 1);
        MiniRecord mr(0, p, 0, 0);
        MinimizerHit mh(1, m, mr);
        EXPECT_EQ(mh.is_forward(), false);
    }

    {
        Minimizer m(min(kh.first, kh.second), 0, 5, 0);
        MiniRecord mr(0, p, 0, 1);
        MinimizerHit mh(1, m, mr);
        EXPECT_EQ(mh.is_forward(), false);
    }
}

TEST(MinimizerHitTest, equals)
{
    KmerHash hash;
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    Minimizer m1(min(kh.first, kh.second), 0, 5, 0);
    deque<Interval> d = { Interval(7, 8), Interval(10, 14) };
    prg::Path p;
    p.initialize(d);
    MiniRecord mr1(0, p, 0, 0);
    MinimizerHit mh1(1, m1, mr1);

    Minimizer m2(min(kh.first, kh.second), 0, 5, 0);
    d = { Interval(7, 9), Interval(11, 14) };
    p.initialize(d);
    MiniRecord mr2(0, p, 0, 0);
    MinimizerHit mh2(1, m2, mr2);

    EXPECT_EQ(mh1, mh1);
    EXPECT_EQ(mh2, mh2);
    EXPECT_EQ((mh1 == mh2), false);
}

TEST(MinimizerHitTest, compare)
{
    set<MinimizerHit> hits;
    KmerHash hash;

    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    Minimizer m12 = Minimizer(min(kh.first, kh.second), 1, 6, 0);
    deque<Interval> d = { Interval(7, 8), Interval(10, 14) };
    prg::Path p;
    p.initialize(d);
    MiniRecord mr12 = MiniRecord(0, p, 0, 0);
    MinimizerHit mh1(1, m12, mr12);
    MinimizerHit mh2(0, m12, mr12);

    Minimizer m3 = Minimizer(min(kh.first, kh.second), 0, 5, 0);
    d = { Interval(6, 10), Interval(11, 12) };
    p.initialize(d);
    MiniRecord mr3 = MiniRecord(0, p, 0, 0);
    MinimizerHit mh3(1, m3, mr3);

    Minimizer m4 = Minimizer(min(kh.first, kh.second), 0, 5, 0);
    d = { Interval(6, 10), Interval(12, 13) };
    p.initialize(d);
    MiniRecord mr4 = MiniRecord(0, p, 0, 0);
    MinimizerHit mh4(1, m4, mr4);

    Minimizer m5 = Minimizer(min(kh.first, kh.second), 0, 5, 0);
    d = { Interval(6, 10), Interval(13, 13), Interval(14, 15) };
    p.initialize(d);
    MiniRecord mr5 = MiniRecord(0, p, 0, 0);
    MinimizerHit mh5(1, m5, mr5);

    Minimizer m6 = Minimizer(min(kh.first, kh.second), 0, 5, 0);
    d = { Interval(6, 10), Interval(14, 14), Interval(14, 15) };
    p.initialize(d);
    MiniRecord mr6 = MiniRecord(0, p, 0, 0);
    MinimizerHit mh6(1, m6, mr6);

    hits.insert(mh1);
    hits.insert(mh2);
    hits.insert(mh3);
    hits.insert(mh4);
    hits.insert(mh5);
    hits.insert(mh6);

    vector<MinimizerHit> expected;
    expected.push_back(mh1);
    expected.push_back(mh2);
    expected.push_back(mh3);
    expected.push_back(mh4);
    expected.push_back(mh5);
    expected.push_back(mh6);

    uint32_t j(1);
    for (set<MinimizerHit>::iterator it = hits.begin(); it != --hits.end(); ++it) {
        EXPECT_EQ(expected[j], *it);
        j++;
    }
    EXPECT_EQ(expected[0], *(--hits.end()));
}
