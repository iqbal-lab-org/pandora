#include "gtest/gtest.h"
#include "kmernode.h"
#include "interval.h"
#include "prg/path.h"
#include <stdint.h>
#include <iostream>


using namespace std;

TEST(KmerNodeTest, create) {

    deque<Interval> d = {Interval(0, 4)};
    Path p;
    p.initialize(d);
    KmerNode kn(0, p);

    uint j = 0;
    EXPECT_EQ(j, kn.id);
    j = 1;
    EXPECT_EQ(j, kn.covg_new.size());
    j = 0;
    EXPECT_EQ(j, kn.get_covg(0, 0));
    EXPECT_EQ(j, kn.get_covg(0, 0));
    EXPECT_EQ(j, kn.num_AT);
    EXPECT_EQ(p, kn.path);
}

TEST(KmerNodeTest, assign) {

    deque<Interval> d = {Interval(0, 4)};
    Path p;
    p.initialize(d);
    KmerNode kn(0, p);
    kn.increment_covg(0,0);
    kn.increment_covg(1,0);
    kn.increment_covg(1,0);

    EXPECT_EQ((uint)0, kn.id);
    EXPECT_EQ((uint)1, kn.covg_new.size());
    EXPECT_EQ((uint)1, kn.get_covg(0, 0));
    EXPECT_EQ((uint)2, kn.get_covg(1, 0));
    EXPECT_EQ((uint)0, kn.num_AT);
    EXPECT_EQ(p, kn.path);

    KmerNode kn_prime = kn;

    EXPECT_EQ((uint)0, kn_prime.id);
    EXPECT_EQ((uint)1, kn_prime.covg_new.size());
    EXPECT_EQ((uint)1, kn_prime.get_covg(0, 0));
    EXPECT_EQ((uint)2, kn_prime.get_covg(1, 0));
    EXPECT_EQ((uint)0, kn_prime.num_AT);
    EXPECT_EQ(p, kn.path);
}

TEST(KmerNodeTest, equals) {

    deque<Interval> d = {Interval(0, 4)};
    Path p1, p2;
    p1.initialize(d);
    d = {Interval(2, 6)};
    p2.initialize(d);
    KmerNode kn1(0, p1);
    KmerNode kn2(3, p1);
    KmerNode kn3(0, p2);
    KmerNode kn4(3, p2);

    // simple cases
    EXPECT_EQ(kn1, kn1);
    EXPECT_EQ(kn2, kn2);
    EXPECT_EQ(kn3, kn3);
    EXPECT_EQ(kn4, kn4);
    EXPECT_EQ((kn1 == kn2), true); //id doesn't affect ==
    EXPECT_EQ((kn1 == kn3), false);
    EXPECT_EQ((kn1 == kn4), false);
    EXPECT_EQ((kn2 == kn1), true);
    EXPECT_EQ((kn2 == kn3), false);
    EXPECT_EQ((kn2 == kn4), false);
    EXPECT_EQ((kn3 == kn1), false);
    EXPECT_EQ((kn3 == kn2), false);
    EXPECT_EQ((kn3 == kn4), true); //id doesn't affect ==

    // covg doesn't affect whether equal
    KmerNode kn5(0, p1);
    kn5.set_covg(6, 0, 0);
    EXPECT_EQ(kn5, kn5);
    EXPECT_EQ(kn1, kn5);
    EXPECT_EQ(kn5, kn1);

}

