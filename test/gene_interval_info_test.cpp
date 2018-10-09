#include "gtest/gtest.h"
#include "gene_interval_info.h"
#include "pangenome/pannode.h"

TEST(GeneIntervalInfoTest,create){
    auto prg_id = 4;
    auto node_id = 3;
    auto name = "dummy";
    auto pn_ptr = make_shared<pangenome::Node>(prg_id, node_id, name);

    GeneIntervalInfo gi{pn_ptr, Interval(0,6), "hello"};

    EXPECT_EQ(gi.pnode, pn_ptr);
    EXPECT_EQ(gi.interval, Interval(0,6));
    EXPECT_EQ(gi.seq, "hello");

    EXPECT_EQ(gi.pnode->prg_id, prg_id);
    EXPECT_EQ(gi.pnode->node_id, node_id);
    EXPECT_EQ(gi.pnode->name, name);
}

TEST(GeneIntervalInfoTest,less_than){
    auto prg_id = 4;
    auto node_id = 3;
    auto name = "dummy";

    auto pn1_ptr = make_shared<pangenome::Node>(prg_id, node_id, name);
    PanNodePtr pn2_ptr = nullptr;

    prg_id = 5;
    node_id = 2;
    auto pn3_ptr = make_shared<pangenome::Node>(prg_id, node_id, name);

    prg_id = 4;
    auto pn4_ptr = make_shared<pangenome::Node>(prg_id, node_id, name);

    GeneIntervalInfo gi1{pn1_ptr, Interval(1,6), "hello"};
    GeneIntervalInfo gi2{pn1_ptr, Interval(1,9), "hello"};
    GeneIntervalInfo gi3{pn1_ptr, Interval(0,9), "hello"};
    GeneIntervalInfo gi4{pn2_ptr, Interval(1,6), "hello"};
    GeneIntervalInfo gi5{pn3_ptr, Interval(1,6), "hello"};
    GeneIntervalInfo gi6{pn4_ptr, Interval(1,6), "hello"};
    
    EXPECT_TRUE(gi1<gi2);
    EXPECT_TRUE(gi3<gi1);
    EXPECT_TRUE(gi4<gi1);
    EXPECT_TRUE(gi1<gi5);
    EXPECT_TRUE(gi6<gi1);

    EXPECT_FALSE(gi2<gi1);
    EXPECT_FALSE(gi1<gi3);
    EXPECT_FALSE(gi1<gi4);
    EXPECT_FALSE(gi5<gi1);
    EXPECT_FALSE(gi1<gi6);

    EXPECT_FALSE(gi1<gi1);
    EXPECT_FALSE(gi2<gi2);
    EXPECT_FALSE(gi3<gi3);
    EXPECT_FALSE(gi4<gi4);
    EXPECT_FALSE(gi5<gi5);
    EXPECT_FALSE(gi6<gi6);
}

TEST(GeneIntervalInfoTest,equals){
    auto prg_id = 4;
    auto node_id = 3;
    auto name = "dummy";

    auto pn1_ptr = make_shared<pangenome::Node>(prg_id, node_id, name);
    PanNodePtr pn2_ptr = nullptr;

    prg_id = 5;
    node_id = 2;
    auto pn3_ptr = make_shared<pangenome::Node>(prg_id, node_id, name);

    prg_id = 4;
    auto pn4_ptr = make_shared<pangenome::Node>(prg_id, node_id, name);

    GeneIntervalInfo gi1{pn1_ptr, Interval(1,6), "hello"};
    GeneIntervalInfo gi2{pn1_ptr, Interval(1,9), "hello"};
    GeneIntervalInfo gi3{pn1_ptr, Interval(0,9), "hello"};
    GeneIntervalInfo gi4{pn2_ptr, Interval(1,6), "hello"};
    GeneIntervalInfo gi5{pn3_ptr, Interval(1,6), "hello"};
    GeneIntervalInfo gi6{pn4_ptr, Interval(1,6), "hello"};

    cout << "defined gii" << endl;

    EXPECT_EQ(pn1_ptr, pn1_ptr);
    EXPECT_EQ(pn2_ptr, pn2_ptr);
    EXPECT_EQ(pn3_ptr, pn3_ptr);
    EXPECT_EQ(pn4_ptr, pn4_ptr);

    EXPECT_NE(pn1_ptr, pn2_ptr);
    EXPECT_NE(pn1_ptr, pn3_ptr);
    EXPECT_NE(pn1_ptr, pn4_ptr);
    EXPECT_NE(pn2_ptr, pn1_ptr);
    EXPECT_NE(pn2_ptr, pn3_ptr);
    EXPECT_NE(pn2_ptr, pn4_ptr);
    EXPECT_NE(pn3_ptr, pn1_ptr);
    EXPECT_NE(pn3_ptr, pn2_ptr);
    EXPECT_NE(pn3_ptr, pn4_ptr);
    EXPECT_NE(pn4_ptr, pn1_ptr);
    EXPECT_NE(pn4_ptr, pn2_ptr);
    EXPECT_NE(pn4_ptr, pn3_ptr);
}
