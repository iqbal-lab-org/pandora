#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "pangenome/ns.cpp"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "pangenome_graph_class.h"
#include "test_helpers.h"
#include "minihit.h"
#include <stdint.h>
#include <iostream>
#include <utility>
#include <functional>
#include <set>

using namespace pangenome;

TEST(PangenomeReadTest, create)
{

    Read pr(3);
    EXPECT_EQ((uint)3, pr.id);
    EXPECT_EQ((uint)0, pr.get_nodes().size());
    EXPECT_EQ((uint)0, pr.node_orientations.size());
    EXPECT_EQ((uint)0, pr.get_hits_as_unordered_map().size());
}

TEST(ReadAddHits, AddOneEmptyClusterToHits_ReadHitsSizeOne)
{
    uint32_t read_id = 1;
    Read read(read_id);
    std::set<MinimizerHitPtr, pComp> cluster;
    uint32_t prg_id = 4;
    auto local_prg_ptr { std::make_shared<LocalPRG>(prg_id, "four", "") };
    PanNodePtr pan_node = make_shared<pangenome::Node>(local_prg_ptr);
    read.add_hits(pan_node, cluster);

    uint result = read.get_hits_as_unordered_map().size();
    uint expect = 1;
    EXPECT_EQ(result, expect);
}

TEST(ReadAddHits, AddOneEmptyClusterToHits_ReadHitsMapContainsCorrectPrgId)
{
    uint32_t read_id = 1;
    Read read(read_id);
    std::set<MinimizerHitPtr, pComp> cluster;
    uint32_t prg_id = 4;
    auto local_prg_ptr { std::make_shared<LocalPRG>(prg_id, "four", "") };
    PanNodePtr pan_node = make_shared<pangenome::Node>(local_prg_ptr);
    read.add_hits(pan_node, cluster);

    auto hits = read.get_hits_as_unordered_map();
    auto result = hits.find(prg_id) != hits.end();
    EXPECT_TRUE(result);
}

TEST(ReadAddHits, AddClusterSecondTime_DeathAndReadHitsNotChanged)
{
    uint32_t read_id = 1;
    Read read(read_id);
    std::set<MinimizerHitPtr, pComp> cluster;
    uint32_t prg_id = 4;

    Interval interval(0, 5);
    std::deque<Interval> raw_path = { Interval(7, 8), Interval(10, 14) };
    prg::Path path;
    path.initialize(raw_path);
    Minimizer m1(0, interval.start, interval.get_end(), 0); // kmer, start, end, strand
    MiniRecord mr1(prg_id, path, 0, 0);
    MinimizerHitPtr minimizer_hit(std::make_shared<MinimizerHit>(read_id, m1, mr1));

    cluster.insert(minimizer_hit);
    auto local_prg_ptr { std::make_shared<LocalPRG>(prg_id, "four", "") };
    PanNodePtr pan_node = make_shared<pangenome::Node>(local_prg_ptr);
    read.add_hits(pan_node, cluster);
    ASSERT_EXCEPTION(read.add_hits(pan_node, cluster), FatalRuntimeError,
        "Error when adding hits to Pangraph read");
    EXPECT_EQ((uint)1, read.get_hits_as_unordered_map()[prg_id].size());
}

TEST(ReadAddHits, AddSecondCluster_ReadHitsMapContainsCorrectPrgIds)
{
    uint32_t read_id = 1;
    Read read(read_id);
    std::set<MinimizerHitPtr, pComp> cluster;
    uint32_t prg_id = 4;
    auto local_prg_ptr_4 { std::make_shared<LocalPRG>(prg_id, "four", "") };
    PanNodePtr pan_node_4 = make_shared<pangenome::Node>(local_prg_ptr_4);
    read.add_hits(pan_node_4, cluster);

    prg_id = 5;
    Interval interval(0, 5);
    std::deque<Interval> raw_path = { Interval(7, 8), Interval(10, 14) };
    prg::Path path;
    path.initialize(raw_path);
    Minimizer m1(0, interval.start, interval.get_end(), 0); // kmer, start, end, strand
    MiniRecord mr1(prg_id, path, 0, 0);
    MinimizerHitPtr minimizer_hit(std::make_shared<MinimizerHit>(read_id, m1, mr1));
    cluster.insert(minimizer_hit);
    auto local_prg_ptr_5 { std::make_shared<LocalPRG>(prg_id, "five", "") };
    PanNodePtr pan_node_5 = make_shared<pangenome::Node>(local_prg_ptr_5);
    read.add_hits(pan_node_5, cluster);
    auto hits = read.get_hits_as_unordered_map();
    auto result = hits.find(prg_id) != hits.end();
    EXPECT_TRUE(result);
}

TEST(PangenomeReadTest, find_position)
{
    std::set<MinimizerHitPtr, pComp> dummy_cluster;
    PGraphTester pg;

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

    // read 0: 0->1->2->3->5->0->7->2->3->5->9
    pg.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l5, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l7, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l5, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l9, 0, dummy_cluster);

    // read 1: 0->1->2
    pg.add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l2, 1, dummy_cluster);

    pg.reads[0]->node_orientations[6] = 1;

    std::vector<uint_least32_t> v = { 2, 3, 5 };
    std::vector<bool> b = { 0, 0, 0 };
    std::pair<uint, uint> p = pg.reads[0]->find_position(v, b);
    std::pair<uint, uint> truth = std::make_pair(2, 4);

    EXPECT_EQ(p.first, truth.first);
    EXPECT_EQ(p.second, truth.second);

    // one at the end of the string
    v = { 3, 5, 9 };
    p = pg.reads[0]->find_position(v, b);
    truth = std::make_pair(8, 10);
    EXPECT_EQ(p.first, truth.first);
    EXPECT_EQ(p.second, truth.second);

    // one in reverse
    v = { 0, 5, 3 };
    b = { 1, 1, 1 };
    p = pg.reads[0]->find_position(v, b);
    truth = std::make_pair(3, 5);
    EXPECT_EQ(p.first, truth.first);
    EXPECT_EQ(p.second, truth.second);

    // one overlapping start
    v = { 9, 0, 1 };
    b = { 0, 0, 0 };
    p = pg.reads[0]->find_position(v, b);
    truth = std::make_pair(0, 1);
    EXPECT_EQ(p.first, truth.first);
    EXPECT_EQ(p.second, truth.second);
    // one in reverse overlapping start
    v = { 1, 0, 9 };
    b = { 1, 1, 1 };
    p = pg.reads[0]->find_position(v, b);
    truth = std::make_pair(0, 1);
    EXPECT_EQ(p.first, truth.first);
    EXPECT_EQ(p.second, truth.second);

    // one overlapping the end
    b = { 0, 0, 0 };
    v = { 5, 9, 9 };
    p = pg.reads[0]->find_position(v, b);
    truth = std::make_pair(9, 10);
    EXPECT_EQ(p.first, truth.first);
    EXPECT_EQ(p.second, truth.second);

    // one in reverse overlapping end
    b = { 1, 1, 1 };
    v = { 0, 9, 5 };
    p = pg.reads[0]->find_position(v, b);
    truth = std::make_pair(9, 10);
    EXPECT_EQ(p.first, truth.first);
    EXPECT_EQ(p.second, truth.second);

    // one not a match
    b = { 0, 0, 0 };
    v = { 8, 8, 8 };
    p = pg.reads[0]->find_position(v, b);
    truth = std::make_pair(
        std::numeric_limits<uint>::max(), std::numeric_limits<uint>::max());
    EXPECT_EQ(p.first, truth.first);
    EXPECT_EQ(p.second, truth.second);

    // one where orientations mean not a match
    v = { 3, 2, 7 };
    p = pg.reads[0]->find_position(v, b);
    EXPECT_EQ(p.first, truth.first);
    EXPECT_EQ(p.second, truth.second);

    // and when is whole read
    v = { 0, 1, 2 };
    p = pg.reads[1]->find_position(v, b);
    truth = std::make_pair(0, 2);
    EXPECT_EQ(p.first, truth.first);
    EXPECT_EQ(p.second, truth.second);
}

TEST(PangenomeReadTest, remove_node)
{

    std::set<MinimizerHitPtr, pComp> dummy_cluster;

    PGraphTester pg;
    std::vector<WeakNodePtr> exp_read_nodes;
    std::vector<bool> exp_read_orientations;

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");

    // read 0: 0->1->2->3
    pg.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    // read 1: -4 -> -3 -> -1
    pg.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 1, dummy_cluster);

    // read 2: 0 -> 1 -> 3 -> 4
    pg.add_hits_between_PRG_and_read(l0, 2, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l4, 2, dummy_cluster);

    // check all as expected
    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[2], pg.nodes[3] };
    exp_read_orientations = { 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4] };
    exp_read_orientations = { 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[2]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example with a node replacing an old node which only appears in one read
    pg.reads[0]->remove_all_nodes_with_this_id(pg.nodes[2]->node_id);

    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[3] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4] };
    exp_read_orientations = { 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[2]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example where old node appears in more than one read
    pg.reads[0]->remove_all_nodes_with_this_id(pg.nodes[1]->node_id);

    exp_read_nodes = { pg.nodes[0], pg.nodes[3] };
    exp_read_orientations = { 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4] };
    exp_read_orientations = { 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[2]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example where have actual hit
    Interval i(0, 5);
    std::deque<Interval> d = { Interval(7, 8), Interval(10, 14) };
    prg::Path p;
    p.initialize(d);
    Minimizer m1(0, i.start, i.get_end(), 0); // kmer, start, end, strand
    MiniRecord mr1(4, p, 0, 0);
    MinimizerHitPtr mh(std::make_shared<MinimizerHit>(4, m1, mr1));
    std::set<MinimizerHitPtr, pComp> c;
    c.insert(mh);
    pg.reads[2]->add_hits(pg.nodes[4], c);

    pg.reads[2]->remove_all_nodes_with_this_id(pg.nodes[4]->node_id);

    exp_read_nodes = { pg.nodes[0], pg.nodes[3] };
    exp_read_orientations = { 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[3] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[2]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example where node appears twice in read
    pg.add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pg.reads[2]->remove_all_nodes_with_this_id(pg.nodes[1]->node_id);

    exp_read_nodes = { pg.nodes[0], pg.nodes[3] };
    exp_read_orientations = { 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[0], pg.nodes[3] };
    exp_read_orientations = { 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[2]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);
}

TEST(PangenomeReadTest, remove_node_it)
{
    std::set<MinimizerHitPtr, pComp> dummy_cluster;

    PGraphTester pg;
    std::vector<WeakNodePtr> exp_read_nodes;
    std::vector<bool> exp_read_orientations;

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");

    // read 0: 0->1->2->3
    pg.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    // read 1: -4 -> -3 -> -1
    pg.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 1, dummy_cluster);

    // read 2: 0 -> 1 -> 3 -> 4
    pg.add_hits_between_PRG_and_read(l0, 2, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l4, 2, dummy_cluster);

    // check all as expected
    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[2], pg.nodes[3] };
    exp_read_orientations = { 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4] };
    exp_read_orientations = { 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[2]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example removing a node which only appears in one read
    auto it = pg.reads[0]->find_node_by_id(pg.nodes[2]->node_id);
    pg.reads[0]->remove_node_with_iterator(it);

    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[3] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4] };
    exp_read_orientations = { 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[2]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example where old node appears in more than one read
    pg.reads[0]->remove_all_nodes_with_this_id(pg.nodes[1]->node_id);

    exp_read_nodes = { pg.nodes[0], pg.nodes[3] };
    exp_read_orientations = { 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[3], pg.nodes[4] };
    exp_read_orientations = { 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[2]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example where have actual hit
    Interval i(0, 5);
    std::deque<Interval> d = { Interval(7, 8), Interval(10, 14) };
    prg::Path p;
    p.initialize(d);
    Minimizer m1(0, i.start, i.get_end(), 0); // kmer, start, end, strand
    MiniRecord mr1(4, p, 0, 0);
    MinimizerHitPtr mh(std::make_shared<MinimizerHit>(4, m1, mr1));
    std::set<MinimizerHitPtr, pComp> c;
    c.insert(mh);
    pg.reads[2]->add_hits(pg.nodes[4], c);

    pg.reads[2]->remove_all_nodes_with_this_id(pg.nodes[4]->node_id);

    exp_read_nodes = { pg.nodes[0], pg.nodes[3] };
    exp_read_orientations = { 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[0], pg.nodes[1], pg.nodes[3] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[2]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);

    // example where node appears twice in read
    pg.add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pg.reads[2]->remove_all_nodes_with_this_id(pg.nodes[1]->node_id);

    exp_read_nodes = { pg.nodes[0], pg.nodes[3] };
    exp_read_orientations = { 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[0]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    exp_read_orientations = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[1]->node_orientations, exp_read_orientations);
    exp_read_nodes = { pg.nodes[0], pg.nodes[3] };
    exp_read_orientations = { 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[2]->get_nodes(), exp_read_nodes, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(
        std::vector<bool>, pg.reads[2]->node_orientations, exp_read_orientations);
}

TEST(PangenomeReadTest, replace_node)
{
    std::set<MinimizerHitPtr, pComp> dummy_cluster;

    pangenome::Graph pg;

    auto l0 = std::make_shared<LocalPRG>(0, "0", "");
    auto l1 = std::make_shared<LocalPRG>(1, "1", "");
    auto l2 = std::make_shared<LocalPRG>(2, "2", "");
    auto l3 = std::make_shared<LocalPRG>(3, "3", "");
    auto l4 = std::make_shared<LocalPRG>(4, "4", "");

    // read 0: 0->1->2->3->1
    pg.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);

    // read 1: 4 -> 3 -> 1
    pg.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pg.add_hits_between_PRG_and_read(l1, 1, dummy_cluster);

    // check what we expect to start with
    EXPECT_EQ((uint)5, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[2]->node_id, (uint)2);
    EXPECT_EQ(pg.nodes[2]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[3]->node_id, (uint)3);
    EXPECT_EQ(pg.nodes[3]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[4]->node_id, (uint)4);
    EXPECT_EQ(pg.nodes[4]->covg, (uint)1);

    EXPECT_EQ((uint)2, pg.reads.size());
    std::vector<WeakNodePtr> read_exp
        = { pg.nodes[0], pg.nodes[1], pg.nodes[2], pg.nodes[3], pg.nodes[1] };
    std::vector<bool> read_o_exp = { 0, 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), read_exp, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(std::vector<bool>, read_o_exp, pg.reads[0]->node_orientations);
    read_exp = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    read_o_exp = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), read_exp, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(std::vector<bool>, read_o_exp, pg.reads[1]->node_orientations);

    // example with a node replacing an old node which only appears in one read
    NodePtr n = std::make_shared<pangenome::Node>(l2, 5);
    pg.nodes[5] = n;
    auto it = pg.reads[0]->get_nodes().begin() + 2;
    pg.reads[0]->replace_node_with_iterator(it, n);

    EXPECT_EQ((uint)6, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[2]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);

    EXPECT_EQ((uint)2, pg.reads.size());
    read_exp = { pg.nodes[0], pg.nodes[1], pg.nodes[5], pg.nodes[3], pg.nodes[1] };
    read_o_exp = { 0, 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), read_exp, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(std::vector<bool>, read_o_exp, pg.reads[0]->node_orientations);
    read_exp = { pg.nodes[4], pg.nodes[3], pg.nodes[1] };
    read_o_exp = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), read_exp, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(std::vector<bool>, read_o_exp, pg.reads[1]->node_orientations);

    // example where old node appears in more than one read
    n = std::make_shared<pangenome::Node>(l3, 6);
    pg.nodes[6] = n;
    it = pg.reads[1]->get_nodes().begin() + 1;
    pg.reads[1]->replace_node_with_iterator(it, n);

    EXPECT_EQ((uint)7, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[2]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[6]->prg_id, (uint)3);

    EXPECT_EQ((uint)2, pg.reads.size());
    read_exp = { pg.nodes[0], pg.nodes[1], pg.nodes[5], pg.nodes[3], pg.nodes[1] };
    read_o_exp = { 0, 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), read_exp, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(std::vector<bool>, read_o_exp, pg.reads[0]->node_orientations);
    read_exp = { pg.nodes[4], pg.nodes[6], pg.nodes[1] };
    read_o_exp = { 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), read_exp, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(std::vector<bool>, read_o_exp, pg.reads[1]->node_orientations);

    // example where move actual hits
    Interval i(0, 5);
    std::deque<Interval> d = { Interval(7, 8), Interval(10, 14) };
    prg::Path p;
    p.initialize(d);
    Minimizer m1(0, i.start, i.get_end(), 0); // kmer, start, end, strand
    MiniRecord mr1(4, p, 0, 0);
    MinimizerHitPtr mh(std::make_shared<MinimizerHit>(4, m1, mr1));
    std::set<MinimizerHitPtr, pComp> c;
    c.insert(mh);
    pg.reads[1]->add_hits(pg.nodes[4], c);
    EXPECT_EQ(pg.reads[1]->get_hits_as_unordered_map()[4].size(), (uint)1);
    n = std::make_shared<pangenome::Node>(l4, 7);
    pg.nodes[7] = n;
    it = pg.reads[1]->get_nodes().begin();
    pg.reads[1]->replace_node_with_iterator(it, n);

    EXPECT_EQ((uint)8, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[2]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[6]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[7]->prg_id, (uint)4);

    EXPECT_EQ((uint)2, pg.reads.size());
    read_exp = { pg.nodes[0], pg.nodes[1], pg.nodes[5], pg.nodes[3], pg.nodes[1] };
    read_o_exp = { 0, 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), read_exp, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(std::vector<bool>, read_o_exp, pg.reads[0]->node_orientations);
    read_exp = { pg.nodes[7], pg.nodes[6], pg.nodes[1], pg.nodes[4] };
    read_o_exp = { 0, 0, 0, 1 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), read_exp, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(std::vector<bool>, read_o_exp, pg.reads[1]->node_orientations);
    EXPECT_EQ(pg.reads[1]->get_hits_as_unordered_map()[7].size(), (uint)0);
    EXPECT_EQ(pg.reads[1]->get_hits_as_unordered_map()[4].size(), (uint)1);

    // example where node appears twice in read
    n = std::make_shared<pangenome::Node>(l1, 8);
    pg.nodes[8] = n;
    it = pg.reads[0]->get_nodes().begin() + 4;
    pg.reads[0]->replace_node_with_iterator(it, n);

    EXPECT_EQ((uint)9, pg.nodes.size());
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[2]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[5]->prg_id, (uint)2);
    EXPECT_EQ(pg.nodes[6]->prg_id, (uint)3);
    EXPECT_EQ(pg.nodes[7]->prg_id, (uint)4);
    EXPECT_EQ(pg.nodes[8]->prg_id, (uint)1);

    EXPECT_EQ((uint)2, pg.reads.size());
    read_exp = { pg.nodes[0], pg.nodes[1], pg.nodes[5], pg.nodes[3], pg.nodes[8] };
    read_o_exp = { 0, 0, 0, 0, 0 };
    EXPECT_TRUE(equal_containers(
        pg.reads[0]->get_nodes(), read_exp, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(std::vector<bool>, read_o_exp, pg.reads[0]->node_orientations);
    read_exp = { pg.nodes[7], pg.nodes[6], pg.nodes[1], pg.nodes[4] };
    read_o_exp = { 0, 0, 0, 1 };
    EXPECT_TRUE(equal_containers(
        pg.reads[1]->get_nodes(), read_exp, EqualComparatorWeakNodePtr()));
    EXPECT_ITERABLE_EQ(std::vector<bool>, read_o_exp, pg.reads[1]->node_orientations);
}

TEST(PangenomeReadTest, equals)
{
    Read pr1(1);
    Read pr2(2);
    EXPECT_EQ(pr1, pr1);
    EXPECT_EQ(pr2, pr2);
    EXPECT_EQ((pr1 == pr2), false);
    EXPECT_EQ((pr2 == pr1), false);
}

TEST(PangenomeReadTest, nequals)
{
    Read pr1(1);
    Read pr2(2);
    EXPECT_NE(pr1, pr2);
    EXPECT_NE(pr2, pr1);
    EXPECT_EQ((pr1 != pr1), false);
    EXPECT_EQ((pr2 != pr2), false);
}

TEST(PangenomeReadTest, less)
{
    Read pr1(1);
    Read pr2(2);
    EXPECT_EQ((pr1 < pr1), false);
    EXPECT_EQ((pr2 < pr2), false);
    EXPECT_EQ((pr1 < pr2), true);
    EXPECT_EQ((pr2 < pr1), false);
}
