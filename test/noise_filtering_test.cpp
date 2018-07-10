#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "noise_filtering.h"
#include "pangenome_graph_class.h"
#include "pangenome/panread.h"
#include "pangenome/pannode.h"
#include "minihit.h"

using namespace std;

TEST(NoiseFilteringTest, node_plus_orientation_to_num) {
    EXPECT_EQ((uint_least32_t) 0, node_plus_orientation_to_num(0, false));
    EXPECT_EQ((uint_least32_t) 1, node_plus_orientation_to_num(0, true));
    EXPECT_EQ((uint_least32_t) 2, node_plus_orientation_to_num(1, false));
    EXPECT_EQ((uint_least32_t) 3, node_plus_orientation_to_num(1, true));
}

TEST(NoiseFilteringTest, num_to_node_plus_orientation) {
    uint_least32_t node_id;
    bool node_orientation;

    num_to_node_plus_orientation(node_id, node_orientation, 0);
    EXPECT_EQ((uint_least32_t) 0, node_id);
    EXPECT_EQ(node_orientation, false);

    num_to_node_plus_orientation(node_id, node_orientation, 1);
    EXPECT_EQ((uint_least32_t) 0, node_id);
    EXPECT_EQ(node_orientation, true);

    num_to_node_plus_orientation(node_id, node_orientation, 2);
    EXPECT_EQ((uint_least32_t) 1, node_id);
    EXPECT_EQ(node_orientation, false);

    num_to_node_plus_orientation(node_id, node_orientation, 3);
    EXPECT_EQ((uint_least32_t) 1, node_id);
    EXPECT_EQ(node_orientation, true);
}

TEST(NoiseFilteringTest, rc_num) {
    EXPECT_EQ(node_plus_orientation_to_num(0, false), rc_num(node_plus_orientation_to_num(0, true)));
    EXPECT_EQ(node_plus_orientation_to_num(0, true), rc_num(node_plus_orientation_to_num(0, false)));
    EXPECT_EQ(node_plus_orientation_to_num(1, false), rc_num(node_plus_orientation_to_num(1, true)));
    EXPECT_EQ(node_plus_orientation_to_num(1, true), rc_num(node_plus_orientation_to_num(1, false)));
    EXPECT_EQ(node_plus_orientation_to_num(2, false), rc_num(node_plus_orientation_to_num(2, true)));
    EXPECT_EQ(node_plus_orientation_to_num(2, true), rc_num(node_plus_orientation_to_num(2, false)));
}

TEST(NoiseFilteringTest, hashed_node_ids_to_ids_and_orientations) {
    std::deque<uint_least32_t> d = {3, 1, 2, 0};
    vector<uint_least32_t> v;
    vector<uint_least32_t> v_exp = {1, 0, 1, 0};
    vector<bool> b;
    vector<bool> b_exp = {true, true, false, false};

    hashed_node_ids_to_ids_and_orientations(d, v, b);
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, v_exp, v);
    EXPECT_ITERABLE_EQ(vector<bool>, b_exp, b);
}

TEST(NoiseFilteringOverlapForwards, SimpleCase_True) {
    std::deque<uint_least32_t> d1 = {0, 1, 2};
    std::deque<uint_least32_t> d2 = {1, 2, 3};
    bool result = overlap_forwards(d1, d2);
    EXPECT_TRUE(result);
    result = overlap_forwards(d2, d1);
    EXPECT_FALSE(result);
}

TEST(NoiseFilteringOverlapForwards, SimpleCase_False) {
    std::deque<uint_least32_t> d1 = {0, 1, 2};
    std::deque<uint_least32_t> d2 = {1, 2, 3};
    bool result = overlap_forwards(d2, d1);
    EXPECT_FALSE(result);
}

TEST(NoiseFilteringOverlapForwards, FirstLongerThanSecond_True) {
    std::deque<uint_least32_t> d1 = {0, 4, 6, 2, 5, 4, 0, 1, 2};
    std::deque<uint_least32_t> d2 = {1, 2, 3};
    bool result = overlap_forwards(d1, d2);
    EXPECT_TRUE(result);
}

TEST(NoiseFilteringOverlapForwards, OverlapShiftedMoreThanOne_False) {
    std::deque<uint_least32_t> d1 = {0, 4, 6, 2, 5, 4, 0, 1, 2};
    std::deque<uint_least32_t> d2 = {1, 2, 3, 4};
    bool result = overlap_forwards(d1, d2);
    EXPECT_FALSE(result);
}

TEST(NoiseFilteringOverlapForwards, SecondLongerThanFirst_Death) {
    std::deque<uint_least32_t> d1 = {0, 4, 6, 2, 5, 4, 0, 1, 2};
    std::deque<uint_least32_t> d2 = {0, 4, 6, 2, 5, 4, 0, 1, 2, 3};
    EXPECT_DEATH(overlap_forwards(d1, d2), "");
}

TEST(NoiseFilteringTest, overlap_backwards) {
    std::deque<uint_least32_t> d1 = {0, 1, 2};
    std::deque<uint_least32_t> d2 = {1, 2, 3};
    EXPECT_EQ(overlap_backwards(d2, d1), true);
    EXPECT_EQ(overlap_backwards(d1, d2), false);
    EXPECT_EQ(overlap_backwards(d1, d1), false);
    EXPECT_EQ(overlap_backwards(d2, d2), false);

    // works when d1 longer than d2
    d1 = {0, 4, 6, 2, 5, 4, 0, 1, 2};
    d2 = {1, 0, 4};
    EXPECT_EQ(overlap_backwards(d1, d2), true);
    EXPECT_EQ(overlap_backwards(d1, d1), false);
    EXPECT_EQ(overlap_backwards(d2, d2), false);


    // if overlap > 1 then is false
    d2 = {1, 2, 0, 4};
    EXPECT_EQ(overlap_backwards(d1, d2), false);
    EXPECT_EQ(overlap_backwards(d2, d2), false);

    // if d2 longer than d1, should still work?
    d2 = {3, 0, 4, 6, 2, 5, 4, 0, 1, 2, 4, 6};
    EXPECT_EQ(overlap_backwards(d1, d2), true);
}

TEST(NoiseFilteringTest, rc_hashed_node_ids) {
    std::deque<uint_least32_t> d1 = {0, 1, 2, 5, 9, 10, 46, 322, 6779};
    std::deque<uint_least32_t> d2 = {6778, 323, 47, 11, 8, 4, 3, 0, 1};
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, d1, rc_hashed_node_ids(d2));
    EXPECT_ITERABLE_EQ(std::deque<uint_least32_t>, d2, rc_hashed_node_ids(d1));
}

TEST(NoiseFilteringDbgNodeIdsToIdsAndOrientations, AllNodesOverlapForward_ConvertCorrectly) {
    uint32_t read_id = 0;

    debruijn::Graph dbg(3);
    std::deque<uint_least32_t> d = {0, 2, 6};
    dbg.add_node(d, read_id);
    d = {2, 6, 11};
    dbg.add_node(d, read_id);
    d = {6, 11, 12};
    dbg.add_node(d, read_id);
    d = {11, 12, 198};
    dbg.add_node(d, read_id);
    d = {12, 198, 60};
    dbg.add_node(d, read_id);
    d = {198, 60, 6};
    dbg.add_node(d, read_id);

    std::deque<uint32_t> tig = {0, 1, 2, 3, 4, 5};
    vector<uint_least32_t> node_ids;
    vector<bool> node_orients;
    dbg_node_ids_to_ids_and_orientations(dbg, tig, node_ids, node_orients);

    vector<uint_least32_t> exp = {0, 1, 3, 5, 6, 99, 30, 3};
    vector<bool> exp_o = {0, 0, 0, 1, 0, 0, 0, 0};
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, exp, node_ids);
    EXPECT_ITERABLE_EQ(vector<bool>, exp_o, node_orients);
}

TEST(NoiseFilteringDbgNodeIdsToIdsAndOrientations, AllNodesOverlapBackward_ConvertCorrectly) {
    uint32_t read_id = 0;

    debruijn::Graph dbg(3);
    std::deque<uint_least32_t> d = {0, 2, 6};
    dbg.add_node(d, read_id);
    d = {2, 6, 11};
    dbg.add_node(d, read_id);
    d = {6, 11, 12};
    dbg.add_node(d, read_id);
    d = {11, 12, 198};
    dbg.add_node(d, read_id);
    d = {12, 198, 60};
    dbg.add_node(d, read_id);
    d = {198, 60, 6};
    dbg.add_node(d, read_id);

    std::deque<uint32_t> tig = {5, 4, 3, 2, 1, 0};
    vector<uint_least32_t> node_ids;
    vector<bool> node_orients;
    dbg_node_ids_to_ids_and_orientations(dbg, tig, node_ids, node_orients);

    vector<uint_least32_t> exp = {0, 1, 3, 5, 6, 99, 30, 3};
    vector<bool> exp_o = {0, 0, 0, 1, 0, 0, 0, 0};
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, exp, node_ids);
    EXPECT_ITERABLE_EQ(vector<bool>, exp_o, node_orients);
}

TEST(NoiseFilteringDbgNodeIdsToIdsAndOrientations, OverlapForwardsSomeReverseComplement_ConvertCorrectly) {

    debruijn::Graph dbg(3);

    // read 0 0->1->3->5->6
    //        0  0  0  1  0
    uint32_t read_id = 0;
    std::deque<uint_least32_t> d = {0, 2, 6};
    dbg.add_node(d, read_id);
    d = {2, 6, 11};
    dbg.add_node(d, read_id);
    d = {6, 11, 12};
    dbg.add_node(d, read_id);

    // read 1 3->30->99->6->5->3
    //        0  0   1   1  0  1
    read_id = 1;
    d = {6, 60, 199};
    dbg.add_node(d, read_id);
    d = {60, 199, 13};
    dbg.add_node(d, read_id);
    d = {199, 13, 10};
    dbg.add_node(d, read_id);
    d = {13, 10, 7};
    dbg.add_node(d, read_id);

    std::deque<uint32_t> tig = {0, 1, 2, 5, 4, 3};
    vector<uint_least32_t> node_ids;
    vector<bool> node_orients;
    dbg_node_ids_to_ids_and_orientations(dbg, tig, node_ids, node_orients);

    vector<uint_least32_t> exp = {0, 1, 3, 5, 6, 99, 30, 3};
    vector<bool> exp_o = {0, 0, 0, 1, 0, 0, 1, 1};
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, exp, node_ids);
    EXPECT_ITERABLE_EQ(vector<bool>, exp_o, node_orients);
}

TEST(NoiseFilteringDbgNodeIdsToIdsAndOrientations, OverlapBackwardsSomeReverseComplement_ConvertCorrectly) {

    debruijn::Graph dbg(3);

    // read 0 0->1->3->5->6
    //        0  0  0  1  0
    uint32_t read_id = 0;
    std::deque<uint_least32_t> d = {0, 2, 6};
    dbg.add_node(d, read_id);
    d = {2, 6, 11};
    dbg.add_node(d, read_id);
    d = {6, 11, 12};
    dbg.add_node(d, read_id);

    // read 1 3->30->99->6->5->3
    //        0  0   1   1  0  1
    read_id = 1;
    d = {6, 60, 199};
    dbg.add_node(d, read_id);
    d = {60, 199, 13};
    dbg.add_node(d, read_id);
    d = {199, 13, 10};
    dbg.add_node(d, read_id);
    d = {13, 10, 7};
    dbg.add_node(d, read_id);

    std::deque<uint32_t> tig = {3, 4, 5, 2, 1, 0};
    vector<uint_least32_t> node_ids;
    vector<bool> node_orients;
    dbg_node_ids_to_ids_and_orientations(dbg, tig, node_ids, node_orients);

    vector<uint_least32_t> exp = {3, 30, 99, 6, 5, 3, 1, 0};
    vector<bool> exp_o = {0, 0, 1, 1, 0, 1, 1, 1};
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, exp, node_ids);
    EXPECT_ITERABLE_EQ(vector<bool>, exp_o, node_orients);
}

TEST(NoiseFilteringTest, construct_debruijn_graph) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    // overlaps to create loop
    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);

    // starts correct and deviates
    pg->add_node(1, "1", 2, mhs);
    pg->add_node(2, "2", 2, mhs);
    pg->add_node(3, "3", 2, mhs);
    pg->add_node(7, "7", 2, mhs);

    // all disjoint, short
    pg->add_node(0, "0", 3, mhs);
    pg->add_node(6, "6", 3, mhs);
    pg->add_node(3, "3", 3, mhs);
    pg->add_node(4, "4", 3, mhs);

    // deviates in middle
    pg->add_node(0, "0", 4, mhs);
    pg->add_node(1, "1", 4, mhs);
    pg->add_node(2, "2", 4, mhs);
    pg->add_node(6, "6", 4, mhs);
    pg->add_node(3, "3", 4, mhs);
    pg->add_node(4, "4", 4, mhs);
    pg->add_node(5, "5", 4, mhs);

    // all disjoint, long
    pg->add_node(6, "6", 5, mhs);
    pg->add_node(1, "1", 5, mhs);
    pg->add_node(2, "2", 5, mhs);
    pg->add_node(6, "6", 5, mhs);
    pg->add_node(3, "3", 5, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);

    debruijn::Graph dbg_exp(3);
    std::deque<uint_least32_t> d = {0, 2, 4};
    debruijn::OrientedNodePtr n1 = dbg_exp.add_node(d, 0);
    d = {2, 4, 6};
    debruijn::OrientedNodePtr n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1, n2);
    d = {4, 6, 8};
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2, n1);
    d = {6, 8, 10};
    n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1, n2);

    d = {6, 8, 10};
    n2 = dbg_exp.add_node(d, 1);
    d = {8, 10, 0};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2, n1);
    d = {10, 0, 2};
    n2 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n1, n2);
    d = {0, 2, 4};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2, n1);

    d = {2, 4, 6};
    n1 = dbg_exp.add_node(d, 2);
    d = {4, 6, 14};
    n2 = dbg_exp.add_node(d, 2);
    dbg_exp.add_edge(n1, n2);

    d = {0, 12, 6};
    n1 = dbg_exp.add_node(d, 3);
    d = {12, 6, 8};
    n2 = dbg_exp.add_node(d, 3);
    dbg_exp.add_edge(n1, n2);

    d = {0, 2, 4};
    n1 = dbg_exp.add_node(d, 4);
    d = {2, 4, 12};
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1, n2);
    d = {4, 12, 6};
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);
    d = {12, 6, 8};
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1, n2);
    d = {6, 8, 10};
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);

    d = {12, 2, 4};
    n1 = dbg_exp.add_node(d, 5);
    d = {2, 4, 12};
    n2 = dbg_exp.add_node(d, 5);
    dbg_exp.add_edge(n1, n2);
    d = {4, 12, 6};
    n1 = dbg_exp.add_node(d, 5);
    dbg_exp.add_edge(n2, n1);

    EXPECT_EQ(dbg_exp, dbg);
    delete pg;
}

TEST(NoiseFilteringRemoveLeaves, OneDBGNode_RemovedFromPanGraph) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();

    // first an example where only one read giving 1 dbg node, want no segfaults
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    remove_leaves(pg, dbg);

    pangenome::Graph pg_exp;
    EXPECT_EQ(pg_exp, *pg);
    delete pg;
}

TEST(NoiseFilteringRemoveLeaves, OneDBGNode_RemovedFromDBGraph) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();

    // first an example where only one read giving 1 dbg node, want no segfaults
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    remove_leaves(pg, dbg);
    debruijn::Graph dbg_exp(3);
    EXPECT_EQ(dbg_exp, dbg);
    delete pg;
}

TEST(NoiseFilteringRemoveLeaves, OneLoop_NoLeavesRemoved) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();

    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    uint pg_size = pg->nodes.size();
    uint dbg_size = dbg.nodes.size();
    remove_leaves(pg, dbg);

    EXPECT_EQ(pg->nodes.size(), pg_size);
    EXPECT_EQ(dbg.nodes.size(), dbg_size);

    delete pg;
}

TEST(NoiseFilteringRemoveLeaves, OneLoopAndDeviantPath_OneLeafRemoved) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();

    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);

    // starts correct and deviates
    pg->add_node(1, "1", 2, mhs);
    pg->add_node(2, "2", 2, mhs);
    pg->add_node(3, "3", 2, mhs);
    pg->add_node(7, "7", 2, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    uint pg_size = pg->nodes.size();
    uint dbg_size = dbg.nodes.size();
    remove_leaves(pg, dbg);

    EXPECT_EQ(pg->nodes.size(), pg_size - 1);
    EXPECT_TRUE(pg->nodes.find(7) == pg->nodes.end());
    EXPECT_EQ(dbg.nodes.size(), dbg_size - 1);
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{4, 6, 14}]) == dbg.nodes.end());

    delete pg;
}

TEST(NoiseFilteringRemoveLeaves, OneLoopAndIncorrectPath_TwoLeavesRemoved) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();

    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);

    // incorrect short
    pg->add_node(0, "0", 3, mhs);
    pg->add_node(5, "5", 3, mhs);//6
    pg->add_node(3, "3", 3, mhs);
    pg->add_node(4, "4", 3, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    uint pg_size = pg->nodes.size();
    uint dbg_size = dbg.nodes.size();
    remove_leaves(pg, dbg);

    EXPECT_EQ(pg->nodes.size(), pg_size);
    EXPECT_EQ(dbg.nodes.size(), dbg_size - 2);
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{0, 10, 6}]) == dbg.nodes.end());
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{10, 6, 8}]) == dbg.nodes.end());

    delete pg;
}

TEST(NoiseFilteringRemoveLeaves, OneLoopAndDeviatesInMiddle_NoLeavesRemoved) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();

    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);

    // deviates in middle
    pg->add_node(0, "0", 4, mhs);
    pg->add_node(1, "1", 4, mhs);
    pg->add_node(2, "2", 4, mhs);
    pg->add_node(6, "6", 4, mhs);
    pg->add_node(3, "3", 4, mhs);
    pg->add_node(4, "4", 4, mhs);
    pg->add_node(5, "5", 4, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    uint pg_size = pg->nodes.size();
    uint dbg_size = dbg.nodes.size();
    remove_leaves(pg, dbg);

    EXPECT_EQ(pg->nodes.size(), pg_size);
    EXPECT_EQ(dbg.nodes.size(), dbg_size);

    delete pg;
}

TEST(NoiseFilteringRemoveLeaves, OneLoopAndLongerWrongPath_LeavesRemoved) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();

    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);

    // incorrect longer
    pg->add_node(6, "6", 5, mhs);
    pg->add_node(1, "1", 5, mhs);
    pg->add_node(7, "7", 5, mhs);//2
    pg->add_node(6, "6", 5, mhs);
    pg->add_node(3, "3", 5, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    uint pg_size = pg->nodes.size();
    uint dbg_size = dbg.nodes.size();
    remove_leaves(pg, dbg);

    EXPECT_EQ(pg->nodes.size(), pg_size - 2);
    EXPECT_TRUE(pg->nodes.find(6) == pg->nodes.end());
    EXPECT_TRUE(pg->nodes.find(7) == pg->nodes.end());
    EXPECT_EQ(dbg.nodes.size(), dbg_size - 3);
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{12, 2, 14}]) == dbg.nodes.end());
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{2, 14, 12}]) == dbg.nodes.end());
    EXPECT_TRUE(dbg.nodes.find(dbg.node_hash[{14, 12, 6}]) == dbg.nodes.end());

    delete pg;
}

TEST(NoiseFilteringRemoveLeaves, AllTogether_GraphsLookCorrect) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();

    // first an example where only one read giving 1 dbg node, want no segfaults
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    remove_leaves(pg, dbg);
    debruijn::Graph dbg_exp(3);
    EXPECT_EQ(dbg_exp, dbg);
    pangenome::Graph pg_exp;
    EXPECT_EQ(pg_exp, *pg);

    // and now a full example
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    // overlapping in loop
    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);

    // starts correct and deviates
    pg->add_node(1, "1", 2, mhs);
    pg->add_node(2, "2", 2, mhs);
    pg->add_node(3, "3", 2, mhs);
    pg->add_node(7, "7", 2, mhs);

    // incorrect short
    pg->add_node(0, "0", 3, mhs);
    pg->add_node(5, "5", 3, mhs);//6
    pg->add_node(3, "3", 3, mhs);
    pg->add_node(4, "4", 3, mhs);

    // deviates in middle
    pg->add_node(0, "0", 4, mhs);
    pg->add_node(1, "1", 4, mhs);
    pg->add_node(2, "2", 4, mhs);
    pg->add_node(6, "6", 4, mhs);
    pg->add_node(3, "3", 4, mhs);
    pg->add_node(4, "4", 4, mhs);
    pg->add_node(5, "5", 4, mhs);

    // incorrect longer
    pg->add_node(6, "6", 5, mhs);
    pg->add_node(1, "1", 5, mhs);
    pg->add_node(1, "1", 5, mhs);//2
    pg->add_node(6, "6", 5, mhs);
    pg->add_node(3, "3", 5, mhs);

    cout << "pg is now: " << endl << *pg << endl;

    construct_debruijn_graph(pg, dbg);
    remove_leaves(pg, dbg);

    std::deque<uint_least32_t> d = {0, 2, 4};
    debruijn::OrientedNodePtr n1 = dbg_exp.add_node(d, 0);
    d = {2, 4, 6};
    debruijn::OrientedNodePtr n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1, n2);
    d = {4, 6, 8};
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2, n1);
    d = {6, 8, 10};
    n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1, n2);

    d = {6, 8, 10};
    n2 = dbg_exp.add_node(d, 1);
    d = {8, 10, 0};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2, n1);
    d = {10, 0, 2};
    n2 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n1, n2);
    d = {0, 2, 4};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2, n1);

    d = {2, 4, 6};
    n1 = dbg_exp.add_node(d, 2);

    d = {0, 2, 4};
    n1 = dbg_exp.add_node(d, 4);
    d = {2, 4, 12};
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1, n2);
    d = {4, 12, 6};
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);
    d = {12, 6, 8};
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1, n2);
    d = {6, 8, 10};
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);

    d = {2, 4, 12};
    n2 = dbg_exp.add_node(d, 4);
    d = {4, 12, 6};
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2, n1);

    EXPECT_EQ(dbg_exp, dbg);

    pg_exp.add_node(0, "0", 0, mhs);
    pg_exp.add_node(1, "1", 0, mhs);
    pg_exp.add_node(2, "2", 0, mhs);
    pg_exp.add_node(3, "3", 0, mhs);
    pg_exp.add_node(4, "4", 0, mhs);
    pg_exp.add_node(5, "5", 0, mhs);

    pg_exp.add_node(3, "3", 1, mhs);
    pg_exp.add_node(4, "4", 1, mhs);
    pg_exp.add_node(5, "5", 1, mhs);
    pg_exp.add_node(0, "0", 1, mhs);
    pg_exp.add_node(1, "1", 1, mhs);
    pg_exp.add_node(2, "2", 1, mhs);

    pg_exp.add_node(1, "1", 2, mhs);
    pg_exp.add_node(2, "2", 2, mhs);
    pg_exp.add_node(3, "3", 2, mhs);

    pg_exp.add_node(0, "0", 4, mhs);
    pg_exp.add_node(1, "1", 4, mhs);
    pg_exp.add_node(2, "2", 4, mhs);
    pg_exp.add_node(6, "6", 4, mhs);
    pg_exp.add_node(3, "3", 4, mhs);
    pg_exp.add_node(4, "4", 4, mhs);
    pg_exp.add_node(5, "5", 4, mhs);

    EXPECT_EQ(pg_exp, *pg);
    delete pg;
}

TEST(NoiseFilteringFilterUnitigs, SimpleCaseNothingToDo_ReadsUnchanged) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);
    pg->add_node(0, "0", 0, mhs);

    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);
    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    filter_unitigs(pg, dbg, 1);

    pangenome::Graph *pg_exp;
    pg_exp = new pangenome::Graph();
    pg_exp->add_node(0, "0", 0, mhs);
    pg_exp->add_node(1, "1", 0, mhs);
    pg_exp->add_node(2, "2", 0, mhs);
    pg_exp->add_node(3, "3", 0, mhs);
    pg_exp->add_node(4, "4", 0, mhs);
    pg_exp->add_node(5, "5", 0, mhs);
    pg_exp->add_node(0, "0", 0, mhs);

    pg_exp->add_node(0, "0", 1, mhs);
    pg_exp->add_node(1, "1", 1, mhs);
    pg_exp->add_node(2, "2", 1, mhs);
    pg_exp->add_node(3, "3", 1, mhs);
    pg_exp->add_node(4, "4", 1, mhs);
    pg_exp->add_node(5, "5", 1, mhs);
    pg_exp->add_node(0, "0", 1, mhs);

    // check contents of read 0 are correct
    EXPECT_EQ(pg_exp->reads[0]->nodes.size(), pg->reads[0]->nodes.size());
    for (uint i = 0; i < min(pg->reads[0]->nodes.size(), pg_exp->reads[0]->nodes.size()); ++i) {
        EXPECT_EQ(*pg->reads[0]->nodes[i], *pg_exp->reads[0]->nodes[i]);
    }
    // check contents of read 1 are correct
    EXPECT_EQ(pg_exp->reads[1]->nodes.size(), pg->reads[1]->nodes.size());
    for (uint i = 0; i < min(pg->reads[1]->nodes.size(), pg_exp->reads[1]->nodes.size()); ++i) {
        EXPECT_EQ(*pg->reads[1]->nodes[i], *pg_exp->reads[1]->nodes[i]);
    }

    delete pg;
    delete pg_exp;
}

TEST(NoiseFilteringFilterUnitigs, SimpleCaseNothingToDoCycle_ReadsUnchanged) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);
    pg->add_node(0, "0", 0, mhs);

    pg->add_node(2, "2", 1, mhs);
    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    filter_unitigs(pg, dbg, 1);

    pangenome::Graph *pg_exp;
    pg_exp = new pangenome::Graph();
    pg_exp->add_node(0, "0", 0, mhs);
    pg_exp->add_node(1, "1", 0, mhs);
    pg_exp->add_node(2, "2", 0, mhs);
    pg_exp->add_node(3, "3", 0, mhs);
    pg_exp->add_node(4, "4", 0, mhs);
    pg_exp->add_node(5, "5", 0, mhs);
    pg_exp->add_node(0, "0", 0, mhs);

    pg_exp->add_node(2, "2", 1, mhs);
    pg_exp->add_node(3, "3", 1, mhs);
    pg_exp->add_node(4, "4", 1, mhs);
    pg_exp->add_node(5, "5", 1, mhs);
    pg_exp->add_node(0, "0", 1, mhs);
    pg_exp->add_node(1, "1", 1, mhs);


    // check contents of read 0 are correct
    EXPECT_EQ(pg_exp->reads[0]->nodes.size(), pg->reads[0]->nodes.size());
    for (uint i = 0; i < min(pg->reads[0]->nodes.size(), pg_exp->reads[0]->nodes.size()); ++i) {
        EXPECT_EQ(*pg->reads[0]->nodes[i], *pg_exp->reads[0]->nodes[i]);
    }
    // check contents of read 1 are correct
    EXPECT_EQ(pg_exp->reads[1]->nodes.size(), pg->reads[1]->nodes.size());
    for (uint i = 0; i < min(pg->reads[1]->nodes.size(), pg_exp->reads[1]->nodes.size()); ++i) {
        EXPECT_EQ(*pg->reads[1]->nodes[i], *pg_exp->reads[1]->nodes[i]);
    }

    delete pg;
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
    EXPECT_EQ(pg_exp->reads[0]->nodes.size(), pg->reads[0]->nodes.size());
    for (uint i=0; i < min(pg->reads[0]->nodes.size(), pg_exp->reads[0]->nodes.size()); ++i)
    {
        EXPECT_EQ(*pg->reads[0]->nodes[i], *pg_exp->reads[0]->nodes[i]);
    }
    // check contents of read 1 are correct
    EXPECT_EQ(pg_exp->reads[1]->nodes.size(), pg->reads[1]->nodes.size());
    for (uint i=0; i < min(pg->reads[1]->nodes.size(),pg_exp->reads[1]->nodes.size()) ; ++i)
    {
        EXPECT_EQ(*pg->reads[1]->nodes[i], *pg_exp->reads[1]->nodes[i]);
    }
    // check contents of read 2 are correct
    EXPECT_EQ(pg_exp->reads[2]->nodes.size(), pg->reads[2]->nodes.size());
    for (uint i=0; i < min(pg->reads[2]->nodes.size(),pg_exp->reads[2]->nodes.size()) ; ++i)
    {
        EXPECT_EQ(*pg->reads[2]->nodes[i], *pg_exp->reads[2]->nodes[i]);
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
    bool result = dbg.nodes.find(5) == dbg.nodes.end(); // node <2,3,7> will be 5th in dbg
    EXPECT_TRUE(result);

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
    bool result = dbg.nodes.find(5) == dbg.nodes.end(); // node <0,5,3> will be 5th in dbg
    EXPECT_TRUE(result);
    result = dbg.nodes.find(6) == dbg.nodes.end(); // node <5,3,4> will be 6th in dbg
    EXPECT_TRUE(result);

    delete pg;
}*/

TEST(NoiseFilteringFilterUnitigs, FilterUnitigsReadDeviatesInMiddle_ReadPruned) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);
    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);

    // deviates in middle
    pg->add_node(0, "0", 4, mhs);
    pg->add_node(1, "1", 4, mhs);
    pg->add_node(2, "2", 4, mhs);
    pg->add_node(6, "6", 4, mhs);
    pg->add_node(3, "3", 4, mhs);
    pg->add_node(4, "4", 4, mhs);
    pg->add_node(5, "5", 4, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    filter_unitigs(pg, dbg, 1);

    pangenome::Graph *pg_exp;
    pg_exp = new pangenome::Graph();
    pg_exp->add_node(0, "0", 0, mhs);
    pg_exp->add_node(1, "1", 0, mhs);
    pg_exp->add_node(2, "2", 0, mhs);
    pg_exp->add_node(3, "3", 0, mhs);
    pg_exp->add_node(4, "4", 0, mhs);
    pg_exp->add_node(5, "5", 0, mhs);

    pg_exp->add_node(0, "0", 1, mhs);
    pg_exp->add_node(1, "1", 1, mhs);
    pg_exp->add_node(2, "2", 1, mhs);
    pg_exp->add_node(3, "3", 1, mhs);
    pg_exp->add_node(4, "4", 1, mhs);
    pg_exp->add_node(5, "5", 1, mhs);

    pg_exp->add_node(0, "0", 4, mhs);
    pg_exp->add_node(1, "1", 4, mhs);
    pg_exp->add_node(2, "2", 4, mhs);
    pg_exp->add_node(3, "3", 4, mhs);
    pg_exp->add_node(4, "4", 4, mhs);
    pg_exp->add_node(5, "5", 4, mhs);

    // check contents of read 4 are correct
    EXPECT_EQ(pg_exp->reads[4]->nodes.size(), pg->reads[4]->nodes.size());
    for (uint i = 0; i < min(pg->reads[4]->nodes.size(), pg_exp->reads[4]->nodes.size()); ++i) {
        EXPECT_EQ(*pg->reads[4]->nodes[i], *pg_exp->reads[4]->nodes[i]);
    }
    delete pg;
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
    bool result = dbg.nodes.find(5) == dbg.nodes.end(); // node <1,2,6> will be 5th in dbg
    EXPECT_TRUE(result);
    result = dbg.nodes.find(6) == dbg.nodes.end(); // node <2,6,3> will be 6th in dbg
    EXPECT_TRUE(result);
    result = dbg.nodes.find(7) == dbg.nodes.end(); // node <6,3,4> will be 7th in dbg
    EXPECT_TRUE(result);

    delete pg;
}*/

TEST(NoiseFilteringFilterUnitigs, FilterUnitigsReadDeviatesLongerInMiddle_ReadPruned) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);
    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);

    // deviates in middle longer
    pg->add_node(0, "0", 5, mhs);
    pg->add_node(1, "1", 5, mhs);
    pg->add_node(2, "2", 5, mhs);
    pg->add_node(9, "9", 5, mhs);
    pg->add_node(10, "10", 5, mhs);
    pg->add_node(11, "11", 5, mhs);
    pg->add_node(3, "3", 5, mhs);
    pg->add_node(4, "4", 5, mhs);
    pg->add_node(5, "5", 5, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    filter_unitigs(pg, dbg, 1);

    pangenome::Graph *pg_exp;
    pg_exp = new pangenome::Graph();
    pg_exp->add_node(0, "0", 0, mhs);
    pg_exp->add_node(1, "1", 0, mhs);
    pg_exp->add_node(2, "2", 0, mhs);
    pg_exp->add_node(3, "3", 0, mhs);
    pg_exp->add_node(4, "4", 0, mhs);
    pg_exp->add_node(5, "5", 0, mhs);

    pg_exp->add_node(0, "0", 1, mhs);
    pg_exp->add_node(1, "1", 1, mhs);
    pg_exp->add_node(2, "2", 1, mhs);
    pg_exp->add_node(3, "3", 1, mhs);
    pg_exp->add_node(4, "4", 1, mhs);
    pg_exp->add_node(5, "5", 1, mhs);

    pg_exp->add_node(0, "0", 5, mhs);
    pg_exp->add_node(1, "1", 5, mhs);
    pg_exp->add_node(2, "2", 5, mhs);
    pg_exp->add_node(3, "3", 5, mhs);
    pg_exp->add_node(4, "4", 5, mhs);
    pg_exp->add_node(5, "5", 5, mhs);

    // check contents of read 4 are correct
    EXPECT_EQ(pg_exp->reads[5]->nodes.size(), pg->reads[5]->nodes.size());
    for (uint i = 0; i < min(pg->reads[5]->nodes.size(), pg_exp->reads[5]->nodes.size()); ++i) {
        EXPECT_EQ(*pg->reads[5]->nodes[i], *pg_exp->reads[5]->nodes[i]);
    }

    delete pg;
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
    bool result = dbg.nodes.find(5) == dbg.nodes.end(); // node <1,2,9> will be 5th in dbg
    EXPECT_TRUE(result);
    result = dbg.nodes.find(6) == dbg.nodes.end(); // node <2,9,10> will be 6th in dbg
    EXPECT_TRUE(result);
    result = dbg.nodes.find(7) == dbg.nodes.end(); // node <9,10,11> will be 7th in dbg
    EXPECT_TRUE(result);
    result = dbg.nodes.find(8) == dbg.nodes.end(); // node <10,11,3> will be 8th in dbg
    EXPECT_TRUE(result);
    result = dbg.nodes.find(9) == dbg.nodes.end(); // node <11,3,4> will be 9th in dbg
    EXPECT_TRUE(result);

    delete pg;
}*/

TEST(NoiseFilteringFilterUnitigs, AllTogether_PanGraphIsAsExpected) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    // starts correct and deviates
    pg->add_node(1, "1", 2, mhs);
    pg->add_node(2, "2", 2, mhs);
    pg->add_node(3, "3", 2, mhs);
    pg->add_node(7, "7", 2, mhs);

    // incorrect short
    pg->add_node(0, "0", 3, mhs);
    pg->add_node(5, "5", 3, mhs);//6
    pg->add_node(3, "3", 3, mhs);
    pg->add_node(4, "4", 3, mhs);

    // deviates in middle
    pg->add_node(0, "0", 4, mhs);
    pg->add_node(1, "1", 4, mhs);
    pg->add_node(2, "2", 4, mhs);
    pg->add_node(6, "6", 4, mhs);
    pg->add_node(3, "3", 4, mhs);
    pg->add_node(4, "4", 4, mhs);
    pg->add_node(5, "5", 4, mhs);

    // deviates in middle longer
    pg->add_node(0, "0", 5, mhs);
    pg->add_node(1, "1", 5, mhs);
    pg->add_node(2, "2", 5, mhs);
    pg->add_node(9, "9", 5, mhs);
    pg->add_node(10, "10", 5, mhs);
    pg->add_node(11, "11", 5, mhs);
    pg->add_node(3, "3", 5, mhs);
    pg->add_node(4, "4", 5, mhs);
    pg->add_node(5, "5", 5, mhs);

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    filter_unitigs(pg, dbg, 1);

    pangenome::Graph pg_exp;
    pg_exp.add_node(0, "0", 0, mhs);
    pg_exp.add_node(1, "1", 0, mhs);
    pg_exp.add_node(2, "2", 0, mhs);
    pg_exp.add_node(3, "3", 0, mhs);
    pg_exp.add_node(4, "4", 0, mhs);
    pg_exp.add_node(5, "5", 0, mhs);


    // starts correct and deviates
    pg_exp.add_node(1, "1", 2, mhs);
    pg_exp.add_node(2, "2", 2, mhs);
    pg_exp.add_node(3, "3", 2, mhs);
    pg_exp.add_node(7, "7", 2, mhs);

    // incorrect short
    pg_exp.add_node(0, "0", 3, mhs);
    pg_exp.add_node(5, "5", 3, mhs);//6
    pg_exp.add_node(3, "3", 3, mhs);
    pg_exp.add_node(4, "4", 3, mhs);

    // deviates in middle
    pg_exp.add_node(0, "0", 4, mhs);
    pg_exp.add_node(1, "1", 4, mhs);
    pg_exp.add_node(2, "2", 4, mhs);
    pg_exp.add_node(3, "3", 4, mhs);
    pg_exp.add_node(4, "4", 4, mhs);
    pg_exp.add_node(5, "5", 4, mhs);

    // deviates in middle longer
    pg_exp.add_node(0, "0", 5, mhs);
    pg_exp.add_node(1, "1", 5, mhs);
    pg_exp.add_node(2, "2", 5, mhs);
    pg_exp.add_node(3, "3", 5, mhs);
    pg_exp.add_node(4, "4", 5, mhs);
    pg_exp.add_node(5, "5", 5, mhs);

    EXPECT_EQ(pg_exp, *pg);
    delete pg;
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

TEST(NoiseFilteringTest, detangle_pangraph_with_debruijn_graph) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();

    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);
    pg->add_node(0, "0", 0, mhs);

    // overlapping in loop
    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);

    // starts correct and deviates
    pg->add_node(1, "1", 2, mhs);
    pg->add_node(2, "2", 2, mhs);
    pg->add_node(3, "3", 2, mhs);
    pg->add_node(7, "7", 2, mhs);

    // incorrect short
    pg->add_node(0, "0", 3, mhs);
    pg->add_node(5, "5", 3, mhs);//6
    pg->add_node(3, "3", 3, mhs);
    pg->add_node(4, "4", 3, mhs);

    // deviates in middle
    pg->add_node(0, "0", 4, mhs);
    pg->add_node(1, "1", 4, mhs);
    pg->add_node(2, "2", 4, mhs);
    pg->add_node(6, "6", 4, mhs);
    pg->add_node(3, "3", 4, mhs);
    pg->add_node(4, "4", 4, mhs);
    pg->add_node(5, "5", 4, mhs);

    cout << "original pg is: " << endl << *pg << endl;

    debruijn::Graph dbg(3);
    construct_debruijn_graph(pg, dbg);
    //detangle_pangraph_with_debruijn_graph(pg, dbg);

    pangenome::Graph pg_exp;
    pangenome::ReadPtr r;
    pangenome::NodePtr n;
    n = make_shared<pangenome::Node>(0, 0, "0");
    pg_exp.nodes[0] = n;
    n = make_shared<pangenome::Node>(1, 1, "1");
    pg_exp.nodes[1] = n;
    n = make_shared<pangenome::Node>(2, 2, "2");
    pg_exp.nodes[2] = n;
    n = make_shared<pangenome::Node>(3, 8, "3");
    pg_exp.nodes[8] = n;
    n = make_shared<pangenome::Node>(4, 9, "4");
    pg_exp.nodes[9] = n;
    n = make_shared<pangenome::Node>(5, 10, "5");
    pg_exp.nodes[10] = n;
    n = make_shared<pangenome::Node>(0, 11, "0");
    pg_exp.nodes[11] = n;
    r = make_shared<pangenome::Read>(0);
    pg_exp.reads[0] = r;
    r->nodes = {pg_exp.nodes[0], pg_exp.nodes[1], pg_exp.nodes[2], pg_exp.nodes[8],
                pg_exp.nodes[9], pg_exp.nodes[10], pg_exp.nodes[11]};

    n = make_shared<pangenome::Node>(1, 12, "1");
    pg_exp.nodes[12] = n;
    n = make_shared<pangenome::Node>(2, 13, "2");
    pg_exp.nodes[13] = n;
    r = make_shared<pangenome::Read>(1);
    pg_exp.reads[1] = r;
    r->nodes = {pg_exp.nodes[8], pg_exp.nodes[9], pg_exp.nodes[10], pg_exp.nodes[11],
                pg_exp.nodes[12], pg_exp.nodes[13]};

    n = make_shared<pangenome::Node>(1, 20, "1");
    pg_exp.nodes[20] = n;
    n = make_shared<pangenome::Node>(2, 21, "2");
    pg_exp.nodes[21] = n;
    n = make_shared<pangenome::Node>(3, 22, "3");
    pg_exp.nodes[22] = n;
    n = make_shared<pangenome::Node>(7, 7, "7");
    pg_exp.nodes[7] = n;
    r = make_shared<pangenome::Read>(2);
    pg_exp.reads[2] = r;
    r->nodes = {pg_exp.nodes[20], pg_exp.nodes[21], pg_exp.nodes[22], pg_exp.nodes[7]};

    n = make_shared<pangenome::Node>(0, 23, "0");
    pg_exp.nodes[23] = n;
    n = make_shared<pangenome::Node>(5, 5, "5");
    pg_exp.nodes[5] = n;
    n = make_shared<pangenome::Node>(3, 3, "3");
    pg_exp.nodes[3] = n;
    n = make_shared<pangenome::Node>(4, 4, "4");
    pg_exp.nodes[4] = n;
    r = make_shared<pangenome::Read>(3);
    pg_exp.reads[3] = r;
    r->nodes = {pg_exp.nodes[23], pg_exp.nodes[5], pg_exp.nodes[3], pg_exp.nodes[4]};

    n = make_shared<pangenome::Node>(0, 14, "0");
    pg_exp.nodes[14] = n;
    n = make_shared<pangenome::Node>(1, 15, "1");
    pg_exp.nodes[15] = n;
    n = make_shared<pangenome::Node>(2, 16, "2");
    pg_exp.nodes[16] = n;
    n = make_shared<pangenome::Node>(6, 6, "6");
    pg_exp.nodes[6] = n;
    n = make_shared<pangenome::Node>(3, 17, "3");
    pg_exp.nodes[17] = n;
    n = make_shared<pangenome::Node>(4, 18, "4");
    pg_exp.nodes[18] = n;
    n = make_shared<pangenome::Node>(5, 19, "5");
    pg_exp.nodes[19] = n;
    r = make_shared<pangenome::Read>(4);
    pg_exp.reads[4] = r;
    r->nodes = {pg_exp.nodes[14], pg_exp.nodes[15], pg_exp.nodes[16], pg_exp.nodes[6],
                pg_exp.nodes[17], pg_exp.nodes[18], pg_exp.nodes[19]};

    //EXPECT_EQ(pg_exp, *pg);
    delete pg;
}

TEST(NoiseFilteringTest, clean_pangraph_with_debruijn_graph) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);

    // starts correct and deviates
    pg->add_node(1, "1", 2, mhs);
    pg->add_node(2, "2", 2, mhs);
    pg->add_node(3, "3", 2, mhs);
    pg->add_node(7, "7", 2, mhs);

    // incorrect short
    pg->add_node(0, "0", 3, mhs);
    pg->add_node(5, "5", 3, mhs);//6
    pg->add_node(3, "3", 3, mhs);
    pg->add_node(4, "4", 3, mhs);

    // deviates in middle
    pg->add_node(0, "0", 4, mhs);
    pg->add_node(1, "1", 4, mhs);
    pg->add_node(2, "2", 4, mhs);
    pg->add_node(6, "6", 4, mhs);
    pg->add_node(3, "3", 4, mhs);
    pg->add_node(4, "4", 4, mhs);
    pg->add_node(5, "5", 4, mhs);

    //clean_pangraph_with_debruijn_graph(pg, 3, 1);

    pangenome::Graph pg_exp;
    pg_exp.add_node(0, "0", 0, mhs);
    pg_exp.add_node(1, "1", 0, mhs);
    pg_exp.add_node(2, "2", 0, mhs);
    pg_exp.add_node(3, "3", 0, mhs);
    pg_exp.add_node(4, "4", 0, mhs);
    pg_exp.add_node(5, "5", 0, mhs);

    // starts correct and deviates
    pg_exp.add_node(1, "1", 2, mhs);
    pg_exp.add_node(2, "2", 2, mhs);
    pg_exp.add_node(3, "3", 2, mhs);

    // deviates in middle
    pg_exp.add_node(0, "0", 4, mhs);
    pg_exp.add_node(1, "1", 4, mhs);
    pg_exp.add_node(2, "2", 4, mhs);
    pg_exp.add_node(3, "3", 4, mhs);
    pg_exp.add_node(4, "4", 4, mhs);
    pg_exp.add_node(5, "5", 4, mhs);

    //EXPECT_EQ(pg_exp, *pg);
    delete pg;
}

TEST(NoiseFilteringTest, write_pangraph_gfa) {
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph *pg;
    pg = new pangenome::Graph();
    pg->add_node(0, "0", 0, mhs);
    pg->add_node(1, "1", 0, mhs);
    pg->add_node(2, "2", 0, mhs);
    pg->add_node(3, "3", 0, mhs);
    pg->add_node(4, "4", 0, mhs);
    pg->add_node(5, "5", 0, mhs);
    pg->add_node(0, "0", 0, mhs);

    // overlapping in loop
    pg->add_node(3, "3", 1, mhs);
    pg->add_node(4, "4", 1, mhs);
    pg->add_node(5, "5", 1, mhs);
    pg->add_node(0, "0", 1, mhs);
    pg->add_node(1, "1", 1, mhs);
    pg->add_node(2, "2", 1, mhs);

    // starts correct and deviates
    pg->add_node(1, "1", 2, mhs);
    pg->add_node(2, "2", 2, mhs);
    pg->add_node(3, "3", 2, mhs);
    pg->add_node(7, "7", 2, mhs);

    // incorrect short
    pg->add_node(0, "0", 3, mhs);
    pg->add_node(5, "5", 3, mhs);//6
    pg->add_node(3, "3", 3, mhs);
    pg->add_node(4, "4", 3, mhs);

    // deviates in middle
    pg->add_node(0, "0", 4, mhs);
    pg->add_node(1, "1", 4, mhs);
    pg->add_node(2, "2", 4, mhs);
    pg->add_node(6, "6", 4, mhs);
    pg->add_node(3, "3", 4, mhs);
    pg->add_node(4, "4", 4, mhs);
    pg->add_node(5, "5", 4, mhs);

    write_pangraph_gfa("../test/test_cases/noisefiltering_test.pangraph.gfa", pg);
}
