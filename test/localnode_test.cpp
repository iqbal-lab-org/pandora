#include "gtest/gtest.h"
#include "localnode.h"
#include "interval.h"
#include <stdint.h>
#include <iostream>

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

TEST(LocalNodeTest, to_string) {
    LocalNode local_node("ACGT", Interval(10, 14), 0);

    std::string actual = local_node.to_string();
    std::string expected {"10 14 ACGT"};

    EXPECT_EQ(actual, expected);
}


TEST(LocalNodeTest, to_string_vector) {
    LocalNode ln1("AAAAA", Interval(0, 5), 0);
    LocalNode ln2("CCCCC", Interval(10, 15), 1);
    LocalNode ln3("GGGGG", Interval(20, 25), 2);
    LocalNode ln4("TTTTT", Interval(30, 35), 3);
    std::vector<LocalNodePtr> local_nodes {
        std::make_shared<LocalNode>(ln1),
        std::make_shared<LocalNode>(ln2),
        std::make_shared<LocalNode>(ln3),
        std::make_shared<LocalNode>(ln4),
    };

    std::string actual = LocalNode::to_string_vector(local_nodes);
    std::string expected = "4 nodes\n(0 [0, 5) AAAAA)\n(1 [10, 15) CCCCC)\n(2 [20, 25) GGGGG)\n(3 [30, 35) TTTTT)";

    EXPECT_EQ(actual, expected);
}
