#include "interval.h"
#include "localnode.h"
#include "gtest/gtest.h"
#include <iostream>
#include <stdint.h>

using namespace std;

TEST(LocalNodeTest, create)
{

    LocalNode ln("ACGTA", Interval(0, 5), 0);

    uint32_t j = 1;
    EXPECT_EQ("ACGTA", ln.seq);
    EXPECT_EQ(Interval(0, 5), ln.pos);
    j = 0;
    EXPECT_EQ(j, ln.id);
}

TEST(LocalNodeTest, equals)
{
    LocalNode ln1("ACGTA", Interval(0, 5), 0);
    LocalNode ln2("AGCTA", Interval(0, 5), 0);
    LocalNode ln3("ACGTA", Interval(0, 4), 0);
    LocalNode ln4("ACGTA", Interval(0, 5), 1);
    // can't compare outNodes bit outside of localGraph
    EXPECT_EQ(ln1, ln1);
    EXPECT_EQ(ln2, ln2);
    EXPECT_EQ(ln3, ln3);
    EXPECT_EQ(ln4, ln4);
    EXPECT_EQ((ln1 == ln2), false);
    // EXPECT_EQ((ln1==ln3), false); //now interval associated with a node does not
    // matter
    EXPECT_EQ((ln1 == ln4), false);
    EXPECT_EQ((ln2 == ln3), false);
    EXPECT_EQ((ln2 == ln4), false);
    EXPECT_EQ((ln3 == ln4), false);
}
