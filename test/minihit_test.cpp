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

TEST(MinimizerHitTest, create) {
    Minimizer m;
    KmerHash hash;
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = Minimizer(min(kh.first, kh.second), 0, 5, 0);
    deque<Interval> d = {Interval(7, 8), Interval(10, 14)};
    prg::Path p;
    p.initialize(d);
    MiniRecord *mr;
    mr = new MiniRecord(0, p, 0, 0);
    MinimizerHit mh(1, m, mr);
    uint32_t j(1);
    EXPECT_EQ(j, mh.read_id);
    EXPECT_EQ((uint) 0, mh.read_start_position);
    j = 0;
    EXPECT_EQ(j, mh.prg_id);
    EXPECT_EQ(p, mh.prg_path);
    bool b = true;
    EXPECT_EQ(b, mh.strand);

    kh = hash.kmerhash("hell", 4);
    m = Minimizer(min(kh.first, kh.second), 1, 5, 0);
    EXPECT_DEATH(MinimizerHit(1, m, mr), "");

    delete mr;
    //TEST SECOND CONSTRUCTOR!!
}

TEST(MinimizerHitTest, checkStrand) {
    Minimizer m;
    KmerHash hash;
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = Minimizer(min(kh.first, kh.second), 0, 5, 0);
    deque<Interval> d = {Interval(7, 8), Interval(10, 14)};
    prg::Path p;
    p.initialize(d);
    MiniRecord *mr;
    mr = new MiniRecord(0, p, 0, 0);
    MinimizerHit mh(1, m, mr);
    EXPECT_EQ(mh.strand, true);

    delete mr;
    m = Minimizer(min(kh.first, kh.second), 0, 5, 1);
    mr = new MiniRecord(0, p, 0, 1);
    MinimizerHit mh1(1, m, mr);
    EXPECT_EQ(mh1.strand, true);

    delete mr;
    m = Minimizer(min(kh.first, kh.second), 0, 5, 1);
    mr = new MiniRecord(0, p, 0, 0);
    MinimizerHit mh2(1, m, mr);
    EXPECT_EQ(mh2.strand, false);

    delete mr;
    m = Minimizer(min(kh.first, kh.second), 0, 5, 0);
    mr = new MiniRecord(0, p, 0, 1);
    MinimizerHit mh3(1, m, mr);
    EXPECT_EQ(mh3.strand, false);

    delete mr;
}

TEST(MinimizerHitTest, equals) {
    Minimizer m;
    KmerHash hash;
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = Minimizer(min(kh.first, kh.second), 0, 5, 0);
    deque<Interval> d = {Interval(7, 8), Interval(10, 14)};
    prg::Path p;
    p.initialize(d);
    MiniRecord *mr;
    mr = new MiniRecord(0, p, 0, 0);
    MinimizerHit mh1(1, m, mr);

    m = Minimizer(min(kh.first, kh.second), 0, 5, 0);
    d = {Interval(7, 9), Interval(11, 14)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    MinimizerHit mh2(1, m, mr);

    EXPECT_EQ(mh1, mh1);
    EXPECT_EQ(mh2, mh2);
    EXPECT_EQ((mh1 == mh2), false);
    delete mr;
}

TEST(MinimizerHitTest, compare) {
    set<MinimizerHit> hits;
    KmerHash hash;

    Minimizer m;
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = Minimizer(min(kh.first, kh.second), 1, 6, 0);
    deque<Interval> d = {Interval(7, 8), Interval(10, 14)};
    prg::Path p;
    p.initialize(d);
    MiniRecord *mr;
    mr = new MiniRecord(0, p, 0, 0);
    MinimizerHit mh1(1, m, mr);

    MinimizerHit mh2(0, m, mr);

    m = Minimizer(min(kh.first, kh.second), 0, 5, 0);

    d = {Interval(6, 10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    MinimizerHit mh3(1, m, mr);

    d = {Interval(6, 10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    MinimizerHit mh4(1, m, mr);

    d = {Interval(6, 10), Interval(13, 13), Interval(14, 15)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    MinimizerHit mh5(1, m, mr);

    d = {Interval(6, 10), Interval(14, 14), Interval(14, 15)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    MinimizerHit mh6(1, m, mr);

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
    delete mr;
}

