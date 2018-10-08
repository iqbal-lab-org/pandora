#include "gtest/gtest.h"
#include "minirecord.h"
#include "prg/path.h"
#include <stdint.h>
#include <iostream>


using namespace std;

TEST(MiniRecordTest, create) {
    deque<Interval> v1 = {Interval(0, 5)};
    deque<Interval> v2 = {Interval(1, 4), Interval(15, 17)};
    deque<Interval> v3 = {Interval(1, 6)};
    deque<Interval> v4 = {Interval(0, 3), Interval(16, 18)};

    Path p;
    p.initialize(v1);
    MiniRecord m1(1, p, 0, 0);
    uint32_t j = 1;
    EXPECT_EQ(j, m1.prg_id);
    EXPECT_EQ(p, m1.path);
    p.initialize(v2);
    MiniRecord m2(2, p, 0, 0);
    j = 2;
    EXPECT_EQ(j, m2.prg_id);
    EXPECT_EQ(p, m2.path);
    p.initialize(v3);
    MiniRecord m3(3, p, 0, 0);
    j = 3;
    EXPECT_EQ(j, m3.prg_id);
    EXPECT_EQ(p, m3.path);
    p.initialize(v4);
    MiniRecord m4(4, p, 0, 0);
    j = 4;
    EXPECT_EQ(j, m4.prg_id);
    EXPECT_EQ(p, m4.path);
}

TEST(MiniRecordTest, equals) {
    deque<Interval> v1 = {Interval(0, 5)};
    deque<Interval> v2 = {Interval(1, 4), Interval(15, 17)};
    deque<Interval> v3 = {Interval(1, 6)};
    deque<Interval> v4 = {Interval(0, 3), Interval(16, 18)};

    Path p;
    p.initialize(v1);
    MiniRecord m1(1, p, 0, 0);
    p.initialize(v2);
    MiniRecord m2(2, p, 0, 0);
    p.initialize(v3);
    MiniRecord m3(3, p, 0, 0);
    p.initialize(v4);
    MiniRecord m4(4, p, 0, 0);

    EXPECT_EQ(m1, m1);
    EXPECT_EQ(m2, m2);
    EXPECT_EQ(m3, m3);
    EXPECT_EQ(m4, m4);
    EXPECT_EQ((m1 == m2), false);
    EXPECT_EQ((m3 == m2), false);
    EXPECT_EQ((m1 == m4), false);
    EXPECT_EQ((m3 == m4), false);
}

TEST(MiniRecordTest, write) {
    deque<Interval> d;
    d = {Interval(1, 3), Interval(4, 5), Interval(6, 6), Interval(9, 40)};
    Path p;
    p.initialize(d);
    MiniRecord mr(1, p, 0, 0);

    stringstream out;
    out << mr;
    EXPECT_EQ(out.str(), "(1, 4{[1, 3)[4, 5)[6, 6)[9, 40)}, 0, 0)");
}

TEST(MiniRecordTest, read) {
    deque<Interval> d;
    d = {Interval(1, 3), Interval(4, 5), Interval(6, 6), Interval(9, 40)};
    Path p;
    p.initialize(d);
    MiniRecord mr1(1, p, 0, 0), mr2;

    stringstream out;
    out << mr1;

    out >> mr2;
    EXPECT_EQ(mr1, mr2);
}
