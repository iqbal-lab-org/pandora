#include "minihit.h"
#include "noise_filtering.h"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "pangenome_graph_class.h"
#include "test_macro.cpp"
#include "gtest/gtest.h"

using namespace std;

TEST(NoiseFilteringTest, node_plus_orientation_to_num)
{
    EXPECT_EQ((uint_least32_t)0, node_plus_orientation_to_num(0, false));
    EXPECT_EQ((uint_least32_t)1, node_plus_orientation_to_num(0, true));
    EXPECT_EQ((uint_least32_t)2, node_plus_orientation_to_num(1, false));
    EXPECT_EQ((uint_least32_t)3, node_plus_orientation_to_num(1, true));
}

TEST(NoiseFilteringTest, num_to_node_plus_orientation)
{
    uint_least32_t node_id;
    bool node_orientation;

    num_to_node_plus_orientation(node_id, node_orientation, 0);
    EXPECT_EQ((uint_least32_t)0, node_id);
    EXPECT_EQ(node_orientation, false);

    num_to_node_plus_orientation(node_id, node_orientation, 1);
    EXPECT_EQ((uint_least32_t)0, node_id);
    EXPECT_EQ(node_orientation, true);

    num_to_node_plus_orientation(node_id, node_orientation, 2);
    EXPECT_EQ((uint_least32_t)1, node_id);
    EXPECT_EQ(node_orientation, false);

    num_to_node_plus_orientation(node_id, node_orientation, 3);
    EXPECT_EQ((uint_least32_t)1, node_id);
    EXPECT_EQ(node_orientation, true);
}

TEST(NoiseFilteringTest, rc_num)
{
    EXPECT_EQ(node_plus_orientation_to_num(0, false),
        rc_num(node_plus_orientation_to_num(0, true)));
    EXPECT_EQ(node_plus_orientation_to_num(0, true),
        rc_num(node_plus_orientation_to_num(0, false)));
    EXPECT_EQ(node_plus_orientation_to_num(1, false),
        rc_num(node_plus_orientation_to_num(1, true)));
    EXPECT_EQ(node_plus_orientation_to_num(1, true),
        rc_num(node_plus_orientation_to_num(1, false)));
    EXPECT_EQ(node_plus_orientation_to_num(2, false),
        rc_num(node_plus_orientation_to_num(2, true)));
    EXPECT_EQ(node_plus_orientation_to_num(2, true),
        rc_num(node_plus_orientation_to_num(2, false)));
}

TEST(NoiseFilteringTest, hashed_node_ids_to_ids_and_orientations)
{
    std::deque<uint_least32_t> d = { 3, 1, 2, 0 };
    vector<uint_least32_t> v;
    vector<uint_least32_t> v_exp = { 1, 0, 1, 0 };
    vector<bool> b;
    vector<bool> b_exp = { true, true, false, false };

    hashed_node_ids_to_ids_and_orientations(d, v, b);
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, v_exp, v);
    EXPECT_ITERABLE_EQ(vector<bool>, b_exp, b);
}

TEST(NoiseFilteringOverlapForwards, SimpleCase_True)
{
    std::deque<uint_least32_t> d1 = { 0, 1, 2 };
    std::deque<uint_least32_t> d2 = { 1, 2, 3 };
    bool result = overlap_forwards(d1, d2);
    EXPECT_TRUE(result);
    result = overlap_forwards(d2, d1);
    EXPECT_FALSE(result);
}

TEST(NoiseFilteringOverlapForwards, SimpleCase_False)
{
    std::deque<uint_least32_t> d1 = { 0, 1, 2 };
    std::deque<uint_least32_t> d2 = { 1, 2, 3 };
    bool result = overlap_forwards(d2, d1);
    EXPECT_FALSE(result);
}

TEST(NoiseFilteringOverlapForwards, FirstLongerThanSecond_True)
{
    std::deque<uint_least32_t> d1 = { 0, 4, 6, 2, 5, 4, 0, 1, 2 };
    std::deque<uint_least32_t> d2 = { 1, 2, 3 };
    bool result = overlap_forwards(d1, d2);
    EXPECT_TRUE(result);
}

TEST(NoiseFilteringOverlapForwards, OverlapShiftedMoreThanOne_False)
{
    std::deque<uint_least32_t> d1 = { 0, 4, 6, 2, 5, 4, 0, 1, 2 };
    std::deque<uint_least32_t> d2 = { 1, 2, 3, 4 };
    bool result = overlap_forwards(d1, d2);
    EXPECT_FALSE(result);
}

TEST(NoiseFilteringOverlapForwards, SecondLongerThanFirst_Death)
{
    std::deque<uint_least32_t> d1 = { 0, 4, 6, 2, 5, 4, 0, 1, 2 };
    std::deque<uint_least32_t> d2 = { 0, 4, 6, 2, 5, 4, 0, 1, 2, 3 };
    EXPECT_DEATH(overlap_forwards(d1, d2), "");
}

TEST(NoiseFilteringTest, overlap_backwards)
{
    std::deque<uint_least32_t> d1 = { 0, 1, 2 };
    std::deque<uint_least32_t> d2 = { 1, 2, 3 };
    EXPECT_EQ(overlap_backwards(d2, d1), true);
    EXPECT_EQ(overlap_backwards(d1, d2), false);
    EXPECT_EQ(overlap_backwards(d1, d1), false);
    EXPECT_EQ(overlap_backwards(d2, d2), false);

    // works when d1 longer than d2
    d1 = { 0, 4, 6, 2, 5, 4, 0, 1, 2 };
    d2 = { 1, 0, 4 };
    EXPECT_EQ(overlap_backwards(d1, d2), true);
    EXPECT_EQ(overlap_backwards(d1, d1), false);
    EXPECT_EQ(overlap_backwards(d2, d2), false);

    // if overlap > 1 then is false
    d2 = { 1, 2, 0, 4 };
    EXPECT_EQ(overlap_backwards(d1, d2), false);
    EXPECT_EQ(overlap_backwards(d2, d2), false);

    // if d2 longer than d1, should still work?
    d2 = { 3, 0, 4, 6, 2, 5, 4, 0, 1, 2, 4, 6 };
    EXPECT_EQ(overlap_backwards(d1, d2), true);
}

TEST(NoiseFilteringTest, rc_hashed_node_ids)
{
    std::deque<uint_least32_t> d1 = { 0, 1, 2, 5, 9, 10, 46, 322, 6779 };
    std::deque<uint_least32_t> d2 = { 6778, 323, 47, 11, 8, 4, 3, 0, 1 };
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, d1, rc_hashed_node_ids(d2));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, d2, rc_hashed_node_ids(d1));
}

TEST(NoiseFilteringDbgNodeIdsToIdsAndOrientations,
    AllNodesOverlapForward_ConvertCorrectly)
{
    uint32_t read_id = 0;

    debruijn::Graph dbg(3);
    std::deque<uint_least32_t> d = { 0, 2, 6 };
    dbg.add_node(d, read_id);
    d = { 2, 6, 11 };
    dbg.add_node(d, read_id);
    d = { 6, 11, 12 };
    dbg.add_node(d, read_id);
    d = { 11, 12, 198 };
    dbg.add_node(d, read_id);
    d = { 12, 198, 60 };
    dbg.add_node(d, read_id);
    d = { 198, 60, 6 };
    dbg.add_node(d, read_id);

    std::deque<uint32_t> tig = { 0, 1, 2, 3, 4, 5 };
    vector<uint_least32_t> node_ids;
    vector<bool> node_orients;
    dbg_node_ids_to_ids_and_orientations(dbg, tig, node_ids, node_orients);

    vector<uint_least32_t> exp = { 0, 1, 3, 5, 6, 99, 30, 3 };
    vector<bool> exp_o = { 0, 0, 0, 1, 0, 0, 0, 0 };
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, exp, node_ids);
    EXPECT_ITERABLE_EQ(vector<bool>, exp_o, node_orients);
}

TEST(NoiseFilteringDbgNodeIdsToIdsAndOrientations,
    AllNodesOverlapBackward_ConvertCorrectly)
{
    uint32_t read_id = 0;

    debruijn::Graph dbg(3);
    std::deque<uint_least32_t> d = { 0, 2, 6 };
    dbg.add_node(d, read_id);
    d = { 2, 6, 11 };
    dbg.add_node(d, read_id);
    d = { 6, 11, 12 };
    dbg.add_node(d, read_id);
    d = { 11, 12, 198 };
    dbg.add_node(d, read_id);
    d = { 12, 198, 60 };
    dbg.add_node(d, read_id);
    d = { 198, 60, 6 };
    dbg.add_node(d, read_id);

    std::deque<uint32_t> tig = { 5, 4, 3, 2, 1, 0 };
    vector<uint_least32_t> node_ids;
    vector<bool> node_orients;
    dbg_node_ids_to_ids_and_orientations(dbg, tig, node_ids, node_orients);

    vector<uint_least32_t> exp = { 0, 1, 3, 5, 6, 99, 30, 3 };
    vector<bool> exp_o = { 0, 0, 0, 1, 0, 0, 0, 0 };
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, exp, node_ids);
    EXPECT_ITERABLE_EQ(vector<bool>, exp_o, node_orients);
}

TEST(NoiseFilteringDbgNodeIdsToIdsAndOrientations,
    OverlapForwardsSomeReverseComplement_ConvertCorrectly)
{

    debruijn::Graph dbg(3);

    // read 0 0->1->3->5->6
    //        0  0  0  1  0
    uint32_t read_id = 0;
    std::deque<uint_least32_t> d = { 0, 2, 6 };
    dbg.add_node(d, read_id);
    d = { 2, 6, 11 };
    dbg.add_node(d, read_id);
    d = { 6, 11, 12 };
    dbg.add_node(d, read_id);

    // read 1 3->30->99->6->5->3
    //        0  0   1   1  0  1
    read_id = 1;
    d = { 6, 60, 199 };
    dbg.add_node(d, read_id);
    d = { 60, 199, 13 };
    dbg.add_node(d, read_id);
    d = { 199, 13, 10 };
    dbg.add_node(d, read_id);
    d = { 13, 10, 7 };
    dbg.add_node(d, read_id);

    std::deque<uint32_t> tig = { 0, 1, 2, 5, 4, 3 };
    vector<uint_least32_t> node_ids;
    vector<bool> node_orients;
    dbg_node_ids_to_ids_and_orientations(dbg, tig, node_ids, node_orients);

    vector<uint_least32_t> exp = { 0, 1, 3, 5, 6, 99, 30, 3 };
    vector<bool> exp_o = { 0, 0, 0, 1, 0, 0, 1, 1 };
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, exp, node_ids);
    EXPECT_ITERABLE_EQ(vector<bool>, exp_o, node_orients);
}

TEST(NoiseFilteringDbgNodeIdsToIdsAndOrientations,
    OverlapBackwardsSomeReverseComplement_ConvertCorrectly)
{

    debruijn::Graph dbg(3);

    // read 0 0->1->3->5->6
    //        0  0  0  1  0
    uint32_t read_id = 0;
    std::deque<uint_least32_t> d = { 0, 2, 6 };
    dbg.add_node(d, read_id);
    d = { 2, 6, 11 };
    dbg.add_node(d, read_id);
    d = { 6, 11, 12 };
    dbg.add_node(d, read_id);

    // read 1 3->30->99->6->5->3
    //        0  0   1   1  0  1
    read_id = 1;
    d = { 6, 60, 199 };
    dbg.add_node(d, read_id);
    d = { 60, 199, 13 };
    dbg.add_node(d, read_id);
    d = { 199, 13, 10 };
    dbg.add_node(d, read_id);
    d = { 13, 10, 7 };
    dbg.add_node(d, read_id);

    std::deque<uint32_t> tig = { 3, 4, 5, 2, 1, 0 };
    vector<uint_least32_t> node_ids;
    vector<bool> node_orients;
    dbg_node_ids_to_ids_and_orientations(dbg, tig, node_ids, node_orients);

    vector<uint_least32_t> exp = { 3, 30, 99, 6, 5, 3, 1, 0 };
    vector<bool> exp_o = { 0, 0, 1, 1, 0, 1, 1, 1 };
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, exp, node_ids);
    EXPECT_ITERABLE_EQ(vector<bool>, exp_o, node_orients);
}

TEST(NoiseFilteringTest, construct_debruijn_graph)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // overlaps to create loop
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    // starts correct and deviates
    pangraph->add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l7, 2, dummy_cluster);

    // all disjoint, short
    pangraph->add_hits_between_PRG_and_read(l0, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l6, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 3, dummy_cluster);

    // deviates in middle
    pangraph->add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l6, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    // all disjoint, long
    pangraph->add_hits_between_PRG_and_read(l6, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l6, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 5, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);

    debruijn::Graph dbg_exp(3);
    std::deque<uint_least32_t> d = { 0, 2, 4 };
    debruijn::OrientedNodePtr n1 = dbg_exp.add_node(d, 0);
    d = { 2, 4, 6 };
    debruijn::OrientedNodePtr n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1, n2);
    d = { 4, 6, 8 };
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2, n1);
    d = { 6, 8, 10 };
    n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1, n2);

    d = { 6, 8, 10 };
    n2 = dbg_exp.add_node(d, 1);
    d = { 8, 10, 0 };
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2, n1);
    d = { 10, 0, 2 };
    n2 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n1, n2);
    d = { 0, 2, 4 };
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2, n1);

    d = { 2, 4, 6 };
    n1 = dbg_exp.add_node(d, 2);
    d = { 4, 6, 14 };
    n2 = dbg_exp.add_node(d, 2);
    dbg_exp.add_edge(n1, n2);

    d = { 0, 12, 6 };
    n1 = dbg_exp.add_node(d, 3);
    d = { 12, 6, 8 };
    n2 = dbg_exp.add_node(d, 3);
    dbg_exp.add_edge(n1, n2);

    d = { 0, 2, 4 };
    n1 = dbg_exp.add_node(d, 4);
    d = { 2, 4, 12 };
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1, n2);
    d = { 4, 12, 6 };
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);
    d = { 12, 6, 8 };
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1, n2);
    d = { 6, 8, 10 };
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);

    d = { 12, 2, 4 };
    n1 = dbg_exp.add_node(d, 5);
    d = { 2, 4, 12 };
    n2 = dbg_exp.add_node(d, 5);
    dbg_exp.add_edge(n1, n2);
    d = { 4, 12, 6 };
    n1 = dbg_exp.add_node(d, 5);
    dbg_exp.add_edge(n2, n1);

    EXPECT_EQ(dbg_exp, dbg);
}

TEST(NoiseFilteringRemoveLeaves, OneDBGNode_RemovedFromPanGraph)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");

    // first an example where only one read giving 1 dbg node, want no segfaults
    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    remove_leaves(pangraph, dbg);

    pangenome::Graph pg_exp;
    EXPECT_EQ(pg_exp, *pangraph);
}

TEST(NoiseFilteringRemoveLeaves, OneDBGNode_RemovedFromDBGraph)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");

    // first an example where only one read giving 1 dbg node, want no segfaults
    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    remove_leaves(pangraph, dbg);
    debruijn::Graph dbg_exp(3);
    EXPECT_EQ(dbg_exp, dbg);
}

TEST(NoiseFilteringRemoveLeaves, OneLoop_NoLeavesRemoved)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // overlaps to create loop
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    uint pg_size = pangraph->nodes.size();
    uint dbg_size = dbg.nodes.size();
    remove_leaves(pangraph, dbg);

    EXPECT_EQ(pangraph->nodes.size(), pg_size);
    EXPECT_EQ(dbg.nodes.size(), dbg_size);
}

TEST(NoiseFilteringRemoveLeaves, OneLoopAndDeviantPath_OneLeafRemoved)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // overlaps to create loop
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    // starts correct and deviates
    pangraph->add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l7, 2, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    uint pg_size = pangraph->nodes.size();
    uint dbg_size = dbg.nodes.size();
    remove_leaves(pangraph, dbg);

    EXPECT_EQ(pangraph->nodes.size(), pg_size - 1);
    EXPECT_TRUE(pangraph->nodes.find(7) == pangraph->nodes.end());
    EXPECT_EQ(dbg.nodes.size(), dbg_size - 1);
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{ 4, 6, 14 }]) == dbg.nodes.end());
}

TEST(NoiseFilteringRemoveLeaves, OneLoopAndIncorrectPath_TwoLeavesRemoved)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // overlaps to create loop
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    // incorrect short
    pangraph->add_hits_between_PRG_and_read(l0, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 3, dummy_cluster); // 6
    pangraph->add_hits_between_PRG_and_read(l3, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 3, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    uint pg_size = pangraph->nodes.size();
    uint dbg_size = dbg.nodes.size();
    remove_leaves(pangraph, dbg);

    EXPECT_EQ(pangraph->nodes.size(), pg_size);
    EXPECT_EQ(dbg.nodes.size(), dbg_size - 2);
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{ 0, 10, 6 }]) == dbg.nodes.end());
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{ 10, 6, 8 }]) == dbg.nodes.end());
}

TEST(NoiseFilteringRemoveLeaves, OneLoopAndDeviatesInMiddle_NoLeavesRemoved)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // overlaps to create loop
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    // deviates in middle
    pangraph->add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l6, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    uint pg_size = pangraph->nodes.size();
    uint dbg_size = dbg.nodes.size();
    remove_leaves(pangraph, dbg);

    EXPECT_EQ(pangraph->nodes.size(), pg_size);
    EXPECT_EQ(dbg.nodes.size(), dbg_size);
}

TEST(NoiseFilteringRemoveLeaves, OneLoopAndLongerWrongPath_LeavesRemoved)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // overlaps to create loop
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    // incorrect longer
    pangraph->add_hits_between_PRG_and_read(l6, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l7, 5, dummy_cluster); // 2
    pangraph->add_hits_between_PRG_and_read(l6, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 5, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    uint pg_size = pangraph->nodes.size();
    uint dbg_size = dbg.nodes.size();
    remove_leaves(pangraph, dbg);

    EXPECT_EQ(pangraph->nodes.size(), pg_size - 2);
    EXPECT_TRUE(pangraph->nodes.find(6) == pangraph->nodes.end());
    EXPECT_TRUE(pangraph->nodes.find(7) == pangraph->nodes.end());
    EXPECT_EQ(dbg.nodes.size(), dbg_size - 3);
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{ 12, 2, 14 }]) == dbg.nodes.end());
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{ 2, 14, 12 }]) == dbg.nodes.end());
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{ 14, 12, 6 }]) == dbg.nodes.end());
}

TEST(NoiseFilteringRemoveLeaves, AllTogether_GraphsLookCorrect)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    // first an example where only one read giving 1 dbg node, want no segfaults
    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    remove_leaves(pangraph, dbg);
    debruijn::Graph dbg_exp(3);
    EXPECT_EQ(dbg_exp, dbg);
    pangenome::Graph pg_exp;
    EXPECT_EQ(pg_exp, *pangraph);

    // and now a full example
    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // overlapping in loop
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    // starts correct and deviates
    pangraph->add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l7, 2, dummy_cluster);

    // incorrect short
    pangraph->add_hits_between_PRG_and_read(l0, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 3, dummy_cluster); // 6
    pangraph->add_hits_between_PRG_and_read(l3, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 3, dummy_cluster);

    // deviates in middle
    pangraph->add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l6, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    // incorrect longer
    pangraph->add_hits_between_PRG_and_read(l6, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 5, dummy_cluster); // 2
    pangraph->add_hits_between_PRG_and_read(l6, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 5, dummy_cluster);

    construct_debruijn_graph(pangraph, dbg);
    remove_leaves(pangraph, dbg);

    std::deque<uint_least32_t> d = { 0, 2, 4 };
    debruijn::OrientedNodePtr n1 = dbg_exp.add_node(d, 0);
    d = { 2, 4, 6 };
    debruijn::OrientedNodePtr n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1, n2);
    d = { 4, 6, 8 };
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2, n1);
    d = { 6, 8, 10 };
    n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1, n2);

    d = { 6, 8, 10 };
    n2 = dbg_exp.add_node(d, 1);
    d = { 8, 10, 0 };
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2, n1);
    d = { 10, 0, 2 };
    n2 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n1, n2);
    d = { 0, 2, 4 };
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2, n1);

    d = { 2, 4, 6 };
    n1 = dbg_exp.add_node(d, 2);

    d = { 0, 2, 4 };
    n1 = dbg_exp.add_node(d, 4);
    d = { 2, 4, 12 };
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1, n2);
    d = { 4, 12, 6 };
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);
    d = { 12, 6, 8 };
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1, n2);
    d = { 6, 8, 10 };
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);

    d = { 2, 4, 12 };
    n2 = dbg_exp.add_node(d, 4);
    d = { 4, 12, 6 };
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);

    EXPECT_EQ(dbg_exp, dbg);

    pg_exp.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    pg_exp.add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    pg_exp.add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 2, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l3, 2, dummy_cluster);

    pg_exp.add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l6, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    EXPECT_EQ(pg_exp, *pangraph);
}

TEST(NoiseFilteringFilterUnitigs, SimpleCaseNothingToDo_ReadsUnchanged)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);

    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    filter_unitigs(pangraph, dbg, 1);

    pangenome::Graph* pg_exp;
    pg_exp = new pangenome::Graph();
    pg_exp->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);

    pg_exp->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);

    // check contents of read 0 are correct
    EXPECT_EQ(
        pg_exp->reads[0]->get_nodes().size(), pangraph->reads[0]->get_nodes().size());
    for (uint i = 0; i < min(pangraph->reads[0]->get_nodes().size(),
                         pg_exp->reads[0]->get_nodes().size());
         ++i) {
        EXPECT_EQ(*pangraph->reads[0]->get_nodes()[i].lock(),
            *pg_exp->reads[0]->get_nodes()[i].lock());
    }
    // check contents of read 1 are correct
    EXPECT_EQ(
        pg_exp->reads[1]->get_nodes().size(), pangraph->reads[1]->get_nodes().size());
    for (uint i = 0; i < min(pangraph->reads[1]->get_nodes().size(),
                         pg_exp->reads[1]->get_nodes().size());
         ++i) {
        EXPECT_EQ(*pangraph->reads[1]->get_nodes()[i].lock(),
            *pg_exp->reads[1]->get_nodes()[i].lock());
    }

    delete pg_exp;
}

TEST(NoiseFilteringFilterUnitigs, SimpleCaseNothingToDoCycle_ReadsUnchanged)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);

    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    filter_unitigs(pangraph, dbg, 1);

    pangenome::Graph* pg_exp;
    pg_exp = new pangenome::Graph();
    pg_exp->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);

    pg_exp->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);

    // check contents of read 0 are correct
    EXPECT_EQ(
        pg_exp->reads[0]->get_nodes().size(), pangraph->reads[0]->get_nodes().size());
    for (uint i = 0; i < min(pangraph->reads[0]->get_nodes().size(),
                         pg_exp->reads[0]->get_nodes().size());
         ++i) {
        EXPECT_EQ(*pangraph->reads[0]->get_nodes()[i].lock(),
            *pg_exp->reads[0]->get_nodes()[i].lock());
    }
    // check contents of read 1 are correct
    EXPECT_EQ(
        pg_exp->reads[1]->get_nodes().size(), pangraph->reads[1]->get_nodes().size());
    for (uint i = 0; i < min(pangraph->reads[1]->get_nodes().size(),
                         pg_exp->reads[1]->get_nodes().size());
         ++i) {
        EXPECT_EQ(*pangraph->reads[1]->get_nodes()[i].lock(),
            *pg_exp->reads[1]->get_nodes()[i].lock());
    }

    delete pg_exp;
}

/*TEST(NoiseFilteringFilterUnitigs,FilterUnitigsReadStartsRightAndDeviates_ReadPruned)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0,"0",0, mhs);
    pg->add_node(1,"1",0, mhs);
    pg->add_node(2,"2",0, mhs);
    pg->add_node(3,"3",0, mhs);
    pg->add_node(4,"4",0, mhs);
    pg->add_node(5,"5",0, mhs);
    pg->add_node(0,"0",0, mhs);

    pg->add_node(0,"0",1, mhs);
    pg->add_node(1,"1",1, mhs);
    pg->add_node(2,"2",1, mhs);
    pg->add_node(3,"3",1, mhs);
    pg->add_node(4,"4",1, mhs);
    pg->add_node(5,"5",1, mhs);
    pg->add_node(0,"0",1, mhs);

    // starts correct and deviates
    pg->add_node(1,"1",2, mhs);
    pg->add_node(2,"2",2, mhs);
    pg->add_node(3,"3",2, mhs);
    pg->add_node(7,"7",2, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph_from_pangraph(pg, dbg);
    filter_unitigs(pg, dbg, 1);

    pangenome::Graph *pg_exp;
    pg_exp = new pangenome::Graph();
    pg_exp->add_node(0,"0",0, mhs);
    pg_exp->add_node(1,"1",0, mhs);
    pg_exp->add_node(2,"2",0, mhs);
    pg_exp->add_node(3,"3",0, mhs);
    pg_exp->add_node(4,"4",0, mhs);
    pg_exp->add_node(5,"5",0, mhs);
    pg_exp->add_node(0,"0",0, mhs);

    pg_exp->add_node(0,"0",1, mhs);
    pg_exp->add_node(1,"1",1, mhs);
    pg_exp->add_node(2,"2",1, mhs);
    pg_exp->add_node(3,"3",1, mhs);
    pg_exp->add_node(4,"4",1, mhs);
    pg_exp->add_node(5,"5",1, mhs);
    pg_exp->add_node(0,"0",1, mhs);

    pg_exp->add_node(1,"1",2, mhs);
    pg_exp->add_node(2,"2",2, mhs);
    pg_exp->add_node(3,"3",2, mhs);

    // check contents of read 0 are correct
    EXPECT_EQ(pg_exp->reads[0]->get_nodes().size(), pg->reads[0]->get_nodes().size());
    for (uint i=0; i < min(pg->reads[0]->get_nodes().size(),
pg_exp->reads[0]->get_nodes().size()); ++i)
    {
        EXPECT_EQ(*pg->reads[0]->get_nodes()[i], *pg_exp->reads[0]->get_nodes()[i]);
    }
    // check contents of read 1 are correct
    EXPECT_EQ(pg_exp->reads[1]->get_nodes().size(), pg->reads[1]->get_nodes().size());
    for (uint i=0; i <
min(pg->reads[1]->get_nodes().size(),pg_exp->reads[1]->get_nodes().size()) ; ++i)
    {
        EXPECT_EQ(*pg->reads[1]->get_nodes()[i], *pg_exp->reads[1]->get_nodes()[i]);
    }
    // check contents of read 2 are correct
    EXPECT_EQ(pg_exp->reads[2]->get_nodes().size(), pg->reads[2]->get_nodes().size());
    for (uint i=0; i <
min(pg->reads[2]->get_nodes().size(),pg_exp->reads[2]->get_nodes().size()) ; ++i)
    {
        EXPECT_EQ(*pg->reads[2]->get_nodes()[i], *pg_exp->reads[2]->get_nodes()[i]);
    }

    delete pg;
    delete pg_exp;
}

TEST(NoiseFilteringFilterUnitigs,FilterUnitigsReadStartsRightAndDeviates_DbgAlsoPruned)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0,"0",0, mhs);
    pg->add_node(1,"1",0, mhs);
    pg->add_node(2,"2",0, mhs);
    pg->add_node(3,"3",0, mhs);
    pg->add_node(4,"4",0, mhs);
    pg->add_node(5,"5",0, mhs);
    pg->add_node(0,"0",0, mhs);

    pg->add_node(0,"0",1, mhs);
    pg->add_node(1,"1",1, mhs);
    pg->add_node(2,"2",1, mhs);
    pg->add_node(3,"3",1, mhs);
    pg->add_node(4,"4",1, mhs);
    pg->add_node(5,"5",1, mhs);
    pg->add_node(0,"0",1, mhs);

    // starts correct and deviates
    pg->add_node(1,"1",2, mhs);
    pg->add_node(2,"2",2, mhs);
    pg->add_node(3,"3",2, mhs);
    pg->add_node(7,"7",2, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    uint num_nodes = dbg.nodes.size();
    filter_unitigs(pg, dbg, 1);

    EXPECT_EQ(num_nodes-1,dbg.nodes.size());
    bool result = dbg.nodes.find(5) == dbg.nodes.end(); // node <2,3,7> will be 5th in
dbg EXPECT_TRUE(result);

    delete pg;
}*/

/*TEST(NoiseFilteringFilterUnitigs,FilterUnitigsReadWrong_ReadRemoved)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0,"0",0, mhs);
    pg->add_node(1,"1",0, mhs);
    pg->add_node(2,"2",0, mhs);
    pg->add_node(3,"3",0, mhs);
    pg->add_node(4,"4",0, mhs);
    pg->add_node(5,"5",0, mhs);
    pg->add_node(0,"0",0, mhs);

    pg->add_node(0,"0",1, mhs);
    pg->add_node(1,"1",1, mhs);
    pg->add_node(2,"2",1, mhs);
    pg->add_node(3,"3",1, mhs);
    pg->add_node(4,"4",1, mhs);
    pg->add_node(5,"5",1, mhs);
    pg->add_node(0,"0",1, mhs);

    // short incorrect
    pg->add_node(0,"0",3, mhs);
    pg->add_node(5,"5",3, mhs);
    pg->add_node(3,"3",3, mhs);
    pg->add_node(4,"4",3, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph_from_pangraph(pg, dbg);
    filter_unitigs(pg, dbg, 1);

    bool found = pg->reads.find(3) != pg->reads.end();

    EXPECT_FALSE(found);
    delete pg;
}

TEST(NoiseFilteringFilterUnitigs,FilterUnitigsReadWrong_DbgAlsoPruned)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0,"0",0, mhs);
    pg->add_node(1,"1",0, mhs);
    pg->add_node(2,"2",0, mhs);
    pg->add_node(3,"3",0, mhs);
    pg->add_node(4,"4",0, mhs);
    pg->add_node(5,"5",0, mhs);
    pg->add_node(0,"0",0, mhs);

    pg->add_node(0,"0",1, mhs);
    pg->add_node(1,"1",1, mhs);
    pg->add_node(2,"2",1, mhs);
    pg->add_node(3,"3",1, mhs);
    pg->add_node(4,"4",1, mhs);
    pg->add_node(5,"5",1, mhs);
    pg->add_node(0,"0",1, mhs);

    // short incorrect
    pg->add_node(0,"0",3, mhs);
    pg->add_node(5,"5",3, mhs);
    pg->add_node(3,"3",3, mhs);
    pg->add_node(4,"4",3, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    uint num_nodes = dbg.nodes.size();
    filter_unitigs(pg, dbg, 1);

    EXPECT_EQ(num_nodes-2,dbg.nodes.size());
    bool result = dbg.nodes.find(5) == dbg.nodes.end(); // node <0,5,3> will be 5th in
dbg EXPECT_TRUE(result); result = dbg.nodes.find(6) == dbg.nodes.end(); // node <5,3,4>
will be 6th in dbg EXPECT_TRUE(result);

    delete pg;
}*/

TEST(NoiseFilteringFilterUnitigs, FilterUnitigsReadDeviatesInMiddle_ReadPruned)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);

    // deviates in middle
    pangraph->add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l6, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    filter_unitigs(pangraph, dbg, 1);

    pangenome::Graph* pg_exp;
    pg_exp = new pangenome::Graph();
    pg_exp->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    pg_exp->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);

    pg_exp->add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    // check contents of read 4 are correct
    EXPECT_EQ(
        pg_exp->reads[4]->get_nodes().size(), pangraph->reads[4]->get_nodes().size());
    for (uint i = 0; i < min(pangraph->reads[4]->get_nodes().size(),
                         pg_exp->reads[4]->get_nodes().size());
         ++i) {
        EXPECT_EQ(*pangraph->reads[4]->get_nodes()[i].lock(),
            *pg_exp->reads[4]->get_nodes()[i].lock());
    }
    delete pg_exp;
}

/*TEST(NoiseFilteringFilterUnitigs,FilterUnitigsDeviatesInMiddle_DbgAlsoPruned)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0,"0",0, mhs);
    pg->add_node(1,"1",0, mhs);
    pg->add_node(2,"2",0, mhs);
    pg->add_node(3,"3",0, mhs);
    pg->add_node(4,"4",0, mhs);
    pg->add_node(5,"5",0, mhs);
    pg->add_node(0,"0",0, mhs);

    pg->add_node(0,"0",1, mhs);
    pg->add_node(1,"1",1, mhs);
    pg->add_node(2,"2",1, mhs);
    pg->add_node(3,"3",1, mhs);
    pg->add_node(4,"4",1, mhs);
    pg->add_node(5,"5",1, mhs);
    pg->add_node(0,"0",1, mhs);

    // deviates in middle
    pg->add_node(0,"0",4, mhs);
    pg->add_node(1,"1",4, mhs);
    pg->add_node(2,"2",4, mhs);
    pg->add_node(6,"6",4, mhs);
    pg->add_node(3,"3",4, mhs);
    pg->add_node(4,"4",4, mhs);
    pg->add_node(5,"5",4, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    uint num_nodes = dbg.nodes.size();
    filter_unitigs(pg, dbg, 1);

    EXPECT_EQ(num_nodes-3,dbg.nodes.size());
    bool result = dbg.nodes.find(5) == dbg.nodes.end(); // node <1,2,6> will be 5th in
dbg EXPECT_TRUE(result); result = dbg.nodes.find(6) == dbg.nodes.end(); // node <2,6,3>
will be 6th in dbg EXPECT_TRUE(result); result = dbg.nodes.find(7) == dbg.nodes.end();
// node <6,3,4> will be 7th in dbg EXPECT_TRUE(result);

    delete pg;
}*/

TEST(NoiseFilteringFilterUnitigs, FilterUnitigsReadDeviatesLongerInMiddle_ReadPruned)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");
    auto l8 = std::make_shared<LocalPRG>(8, "8", "");
    auto l9 = std::make_shared<LocalPRG>(9, "9", "");
    auto l10 = std::make_shared<LocalPRG>(10, "10", "");
    auto l11 = std::make_shared<LocalPRG>(11, "11", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);

    // deviates in middle longer
    pangraph->add_hits_between_PRG_and_read(l0, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l9, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l10, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l11, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 5, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    filter_unitigs(pangraph, dbg, 1);

    pangenome::Graph* pg_exp;
    pg_exp = new pangenome::Graph();
    pg_exp->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    pg_exp->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);

    pg_exp->add_hits_between_PRG_and_read(l0, 5, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l1, 5, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l2, 5, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l3, 5, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l4, 5, dummy_cluster);
    pg_exp->add_hits_between_PRG_and_read(l5, 5, dummy_cluster);

    // check contents of read 4 are correct
    EXPECT_EQ(
        pg_exp->reads[5]->get_nodes().size(), pangraph->reads[5]->get_nodes().size());
    for (uint i = 0; i < min(pangraph->reads[5]->get_nodes().size(),
                         pg_exp->reads[5]->get_nodes().size());
         ++i) {
        EXPECT_EQ(*pangraph->reads[5]->get_nodes()[i].lock(),
            *pg_exp->reads[5]->get_nodes()[i].lock());
    }

    delete pg_exp;
}

/*TEST(NoiseFilteringFilterUnitigs,FilterUnitigsDeviatesLongerInMiddle_DbgAlsoPruned)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0,"0",0, mhs);
    pg->add_node(1,"1",0, mhs);
    pg->add_node(2,"2",0, mhs);
    pg->add_node(3,"3",0, mhs);
    pg->add_node(4,"4",0, mhs);
    pg->add_node(5,"5",0, mhs);
    pg->add_node(0,"0",0, mhs);

    pg->add_node(0,"0",1, mhs);
    pg->add_node(1,"1",1, mhs);
    pg->add_node(2,"2",1, mhs);
    pg->add_node(3,"3",1, mhs);
    pg->add_node(4,"4",1, mhs);
    pg->add_node(5,"5",1, mhs);
    pg->add_node(0,"0",1, mhs);

    // deviates in middle longer
    pg->add_node(0,"0",5, mhs);
    pg->add_node(1,"1",5, mhs);
    pg->add_node(2,"2",5, mhs);
    pg->add_node(9,"9",5, mhs);
    pg->add_node(10,"10",5, mhs);
    pg->add_node(11,"11",5, mhs);
    pg->add_node(3,"3",5, mhs);
    pg->add_node(4,"4",5, mhs);
    pg->add_node(5,"5",5, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    uint num_nodes = dbg.nodes.size();
    filter_unitigs(pg, dbg, 1);

    EXPECT_EQ(num_nodes-5,dbg.nodes.size());
    bool result = dbg.nodes.find(5) == dbg.nodes.end(); // node <1,2,9> will be 5th in
dbg EXPECT_TRUE(result); result = dbg.nodes.find(6) == dbg.nodes.end(); // node <2,9,10>
will be 6th in dbg EXPECT_TRUE(result); result = dbg.nodes.find(7) == dbg.nodes.end();
// node <9,10,11> will be 7th in dbg EXPECT_TRUE(result); result = dbg.nodes.find(8) ==
dbg.nodes.end(); // node <10,11,3> will be 8th in dbg EXPECT_TRUE(result); result =
dbg.nodes.find(9) == dbg.nodes.end(); // node <11,3,4> will be 9th in dbg
    EXPECT_TRUE(result);

    delete pg;
}*/

TEST(NoiseFilteringFilterUnitigs, AllTogether_PanGraphIsAsExpected)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");
    auto l8 = std::make_shared<LocalPRG>(8, "8", "");
    auto l9 = std::make_shared<LocalPRG>(9, "9", "");
    auto l10 = std::make_shared<LocalPRG>(10, "10", "");
    auto l11 = std::make_shared<LocalPRG>(11, "11", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // starts correct and deviates
    pangraph->add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l7, 2, dummy_cluster);

    // incorrect short
    pangraph->add_hits_between_PRG_and_read(l0, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 3, dummy_cluster); // 6
    pangraph->add_hits_between_PRG_and_read(l3, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 3, dummy_cluster);

    // deviates in middle
    pangraph->add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l6, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    // deviates in middle longer
    pangraph->add_hits_between_PRG_and_read(l0, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l9, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l10, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l11, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 5, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 5, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    filter_unitigs(pangraph, dbg, 1);

    pangenome::Graph pg_exp;
    pg_exp.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // starts correct and deviates
    pg_exp.add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 2, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l7, 2, dummy_cluster);

    // incorrect short
    pg_exp.add_hits_between_PRG_and_read(l0, 3, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l5, 3, dummy_cluster); // 6
    pg_exp.add_hits_between_PRG_and_read(l3, 3, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l4, 3, dummy_cluster);

    // deviates in middle
    pg_exp.add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    // deviates in middle longer
    pg_exp.add_hits_between_PRG_and_read(l0, 5, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l1, 5, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 5, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l3, 5, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l4, 5, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l5, 5, dummy_cluster);

    EXPECT_EQ(pg_exp, *pangraph);
}

/*TEST(NoiseFilteringFilterUnitigs,AllTogether_DbgIsAsExpected)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0,"0",0, mhs);
    pg->add_node(1,"1",0, mhs);
    pg->add_node(2,"2",0, mhs);
    pg->add_node(3,"3",0, mhs);
    pg->add_node(4,"4",0, mhs);
    pg->add_node(5,"5",0, mhs);
    pg->add_node(0,"0",0, mhs);

    // overlapping in loop
    pg->add_node(3,"3",1, mhs);
    pg->add_node(4,"4",1, mhs);
    pg->add_node(5,"5",1, mhs);
    pg->add_node(0,"0",1, mhs);
    pg->add_node(1,"1",1, mhs);
    pg->add_node(2,"2",1, mhs);

    // starts correct and deviates
    pg->add_node(1,"1",2, mhs);
    pg->add_node(2,"2",2, mhs);
    pg->add_node(3,"3",2, mhs);
    pg->add_node(7,"7",2, mhs);

    // incorrect short
    pg->add_node(0,"0",3, mhs);
    pg->add_node(5,"5",3, mhs);//6
    pg->add_node(3,"3",3, mhs);
    pg->add_node(4,"4",3, mhs);

    // deviates in middle
    pg->add_node(0,"0",4, mhs);
    pg->add_node(1,"1",4, mhs);
    pg->add_node(2,"2",4, mhs);
    pg->add_node(6,"6",4, mhs);
    pg->add_node(3,"3",4, mhs);
    pg->add_node(4,"4",4, mhs);
    pg->add_node(5,"5",4, mhs);

    // deviates in middle longer
    pg->add_node(0,"0",5, mhs);
    pg->add_node(1,"1",5, mhs);
    pg->add_node(2,"2",5, mhs);
    pg->add_node(9,"9",5, mhs);
    pg->add_node(10,"10",5, mhs);
    pg->add_node(11,"11",5, mhs);
    pg->add_node(3,"3",5, mhs);
    pg->add_node(4,"4",5, mhs);
    pg->add_node(5,"5",5, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    filter_unitigs(pg, dbg, 1);

    debruijn::Graph dbg_exp(3);
    std::deque<uint_least32_t> d = {0,2,4};
    debruijn::OrientedNodePtr n1 = dbg_exp.add_node(d, 0);
    d = {2,4,6};
    debruijn::OrientedNodePtr n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1,n2);
    d = {4,6,8};
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2,n1);
    d = {6,8,10};
    n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1,n2);
    d = {8,10,0};
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2,n1);

    d = {6,8,10};
    n2 = dbg_exp.add_node(d, 1);
    d = {8,10,0};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2,n1);
    d = {10,0,2};
    n2 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n1,n2);
    d = {0,2,4};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2,n1);

    d = {2,4,6};
    n1 = dbg_exp.add_node(d, 2);
    d = {4,6,14};
    n2 = dbg_exp.add_node(d, 2);
    dbg_exp.add_edge(n1,n2);

    d = {0,10,6};
    n1 = dbg_exp.add_node(d, 3);
    d = {10,6,8};
    n2 = dbg_exp.add_node(d, 3);
    dbg_exp.add_edge(n1,n2);

    d = {0,2,4};
    n1 = dbg_exp.add_node(d, 4);
    d = {6,8,10};
    n2 = dbg_exp.add_node(d, 4);

    EXPECT_EQ(dbg_exp, dbg);

    delete pg;
}*/

TEST(NoiseFilteringTest, detangle_pangraph_with_debruijn_graph)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);

    // overlapping in loop
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    // starts correct and deviates
    pangraph->add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l7, 2, dummy_cluster);

    // incorrect short
    pangraph->add_hits_between_PRG_and_read(l0, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 3, dummy_cluster); // 6
    pangraph->add_hits_between_PRG_and_read(l3, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 3, dummy_cluster);

    // deviates in middle
    pangraph->add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l6, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pangraph, dbg);
    // detangle_pangraph_with_debruijn_graph(pg, dbg);

    pangenome::Graph pg_exp;
    pangenome::ReadPtr r;
    pangenome::NodePtr n;
    n = make_shared<pangenome::Node>(l0, 0);
    pg_exp.nodes[0] = n;
    n = make_shared<pangenome::Node>(l1, 1);
    pg_exp.nodes[1] = n;
    n = make_shared<pangenome::Node>(l2, 2);
    pg_exp.nodes[2] = n;
    n = make_shared<pangenome::Node>(l3, 8);
    pg_exp.nodes[8] = n;
    n = make_shared<pangenome::Node>(l4, 9);
    pg_exp.nodes[9] = n;
    n = make_shared<pangenome::Node>(l5, 10);
    pg_exp.nodes[10] = n;
    n = make_shared<pangenome::Node>(l0, 11);
    pg_exp.nodes[11] = n;
    r = make_shared<pangenome::Read>(0);
    pg_exp.reads[0] = r;
    r->set_nodes({ pg_exp.nodes[0], pg_exp.nodes[1], pg_exp.nodes[2], pg_exp.nodes[8],
        pg_exp.nodes[9], pg_exp.nodes[10], pg_exp.nodes[11] });

    n = make_shared<pangenome::Node>(l1, 12);
    pg_exp.nodes[12] = n;
    n = make_shared<pangenome::Node>(l2, 13);
    pg_exp.nodes[13] = n;
    r = make_shared<pangenome::Read>(1);
    pg_exp.reads[1] = r;
    r->set_nodes({ pg_exp.nodes[8], pg_exp.nodes[9], pg_exp.nodes[10], pg_exp.nodes[11],
        pg_exp.nodes[12], pg_exp.nodes[13] });

    n = make_shared<pangenome::Node>(l1, 20);
    pg_exp.nodes[20] = n;
    n = make_shared<pangenome::Node>(l2, 21);
    pg_exp.nodes[21] = n;
    n = make_shared<pangenome::Node>(l3, 22);
    pg_exp.nodes[22] = n;
    n = make_shared<pangenome::Node>(l7, 7);
    pg_exp.nodes[7] = n;
    r = make_shared<pangenome::Read>(2);
    pg_exp.reads[2] = r;
    r->set_nodes(
        { pg_exp.nodes[20], pg_exp.nodes[21], pg_exp.nodes[22], pg_exp.nodes[7] });

    n = make_shared<pangenome::Node>(l0, 23);
    pg_exp.nodes[23] = n;
    n = make_shared<pangenome::Node>(l5, 5);
    pg_exp.nodes[5] = n;
    n = make_shared<pangenome::Node>(l3, 3);
    pg_exp.nodes[3] = n;
    n = make_shared<pangenome::Node>(l4, 4);
    pg_exp.nodes[4] = n;
    r = make_shared<pangenome::Read>(3);
    pg_exp.reads[3] = r;
    r->set_nodes(
        { pg_exp.nodes[23], pg_exp.nodes[5], pg_exp.nodes[3], pg_exp.nodes[4] });

    n = make_shared<pangenome::Node>(l0, 14);
    pg_exp.nodes[14] = n;
    n = make_shared<pangenome::Node>(l1, 15);
    pg_exp.nodes[15] = n;
    n = make_shared<pangenome::Node>(l2, 16);
    pg_exp.nodes[16] = n;
    n = make_shared<pangenome::Node>(l6, 6);
    pg_exp.nodes[6] = n;
    n = make_shared<pangenome::Node>(l3, 17);
    pg_exp.nodes[17] = n;
    n = make_shared<pangenome::Node>(l4, 18);
    pg_exp.nodes[18] = n;
    n = make_shared<pangenome::Node>(l5, 19);
    pg_exp.nodes[19] = n;
    r = make_shared<pangenome::Read>(4);
    pg_exp.reads[4] = r;
    r->set_nodes({ pg_exp.nodes[14], pg_exp.nodes[15], pg_exp.nodes[16],
        pg_exp.nodes[6], pg_exp.nodes[17], pg_exp.nodes[18], pg_exp.nodes[19] });

    // EXPECT_EQ(pg_exp, *pg);
}

TEST(NoiseFilteringTest, clean_pangraph_with_debruijn_graph)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // starts correct and deviates
    pangraph->add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l7, 2, dummy_cluster);

    // incorrect short
    pangraph->add_hits_between_PRG_and_read(l0, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 3, dummy_cluster); // 6
    pangraph->add_hits_between_PRG_and_read(l3, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 3, dummy_cluster);

    // deviates in middle
    pangraph->add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l6, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    // clean_pangraph_with_debruijn_graph(pg, 3, 1);

    pangenome::Graph pg_exp;
    pg_exp.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l5, 0, dummy_cluster);

    // starts correct and deviates
    pg_exp.add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 2, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l3, 2, dummy_cluster);

    // deviates in middle
    pg_exp.add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pg_exp.add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    // EXPECT_EQ(pg_exp, *pg);
}

TEST(NoiseFilteringTest, write_pangraph_gfa)
{
    set<MinimizerHitPtr, pComp> dummy_cluster;
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");
    auto l5 = std::make_shared<LocalPRG>(5, "5", "");
    auto l6 = std::make_shared<LocalPRG>(6, "6", "");
    auto l7 = std::make_shared<LocalPRG>(7, "7", "");

    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 0, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 0, dummy_cluster);

    // overlapping in loop
    pangraph->add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    // starts correct and deviates
    pangraph->add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l7, 2, dummy_cluster);

    // incorrect short
    pangraph->add_hits_between_PRG_and_read(l0, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 3, dummy_cluster); // 6
    pangraph->add_hits_between_PRG_and_read(l3, 3, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 3, dummy_cluster);

    // deviates in middle
    pangraph->add_hits_between_PRG_and_read(l0, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l1, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l2, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l6, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l3, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l4, 4, dummy_cluster);
    pangraph->add_hits_between_PRG_and_read(l5, 4, dummy_cluster);

    write_pangraph_gfa("../test/test_cases/noisefiltering_test.pangraph.gfa", pangraph);
}
