#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "pangenome/ns.cpp"
#include "pangenome_graph_class.h"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "pangenome/pansample.h"
#include "minihit.h"
#include "localPRG.h"
#include <stdint.h>
#include <numeric>
#include <cassert>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>


using namespace pangenome;

TEST(PangenomeGraphGetRead, AddRead_PangenomeGraphReadsContainsReadId) {
    uint32_t read_id = 2;
    PGraphTester pg;

    auto result = pg.reads.find(read_id) == pg.reads.end();
    EXPECT_TRUE(result);
    EXPECT_EQ(pg.reads.size(), (uint) 0);

    pg.get_read(read_id);

    result = pg.reads.find(read_id) != pg.reads.end();
    EXPECT_TRUE(result);
    EXPECT_EQ(pg.reads.size(), (uint) 1);
}

TEST(PangenomeGraphGetRead, AddReadTwice_PangenomeGraphReadsContainsReadId) {
    uint32_t read_id = 2;
    PGraphTester pg;

    pg.get_read(read_id);
    pg.get_read(read_id);

    auto result = pg.reads.find(read_id) != pg.reads.end();
    EXPECT_TRUE(result);
    EXPECT_EQ(pg.reads.size(), (uint) 1);
}

TEST(PangenomeGraphAddCoverage, NodeDoesntAlreadyExist_PangenomeGraphNodesContainsNodeId) {
    PGraphTester pg;

    uint32_t read_id = 2;
    ReadPtr read_ptr = pg.get_read(read_id);

    std::set<MinimizerHitPtr, pComp> mhs;
    uint32_t node_id = 0;
    uint32_t prg_id = 1;

    EXPECT_EQ(pg.nodes.size(), (uint) 0);

    pg.add_coverage(read_ptr, node_id, prg_id, "0");

    auto result = pg.nodes.find(node_id) != pg.nodes.end();
    EXPECT_TRUE(result);

    result = pg.nodes.find(prg_id) == pg.nodes.end();
    EXPECT_TRUE(result);
    EXPECT_EQ(pg.reads.size(), (uint) 1);
}

TEST(PangenomeGraphAddCoverage, NodeDoesntAlreadyExist_PangenomeGraphNodeContainsReadPtr) {
    PGraphTester pg;

    uint32_t read_id = 2;
    ReadPtr read_ptr = pg.get_read(read_id);

    std::set<MinimizerHitPtr, pComp> mhs;
    uint32_t node_id = 0;
    uint32_t prg_id = 1;

    NodePtr node_ptr = pg.add_coverage(read_ptr, node_id, prg_id, "0");

    auto result = node_ptr->reads.find(read_ptr) != node_ptr->reads.end();
    EXPECT_TRUE(result);
}

TEST(PangenomeGraphAddCoverage, NodeAlreadyExists_PangenomeGraphNodeCoverageIncreases) {
    PGraphTester pg;

    uint32_t read_id = 2;
    ReadPtr read_ptr = pg.get_read(read_id);

    std::set<MinimizerHitPtr, pComp> mhs;
    uint32_t node_id = 0;
    uint32_t prg_id = 1;

    EXPECT_EQ(pg.nodes.size(), (uint) 0);

    NodePtr node_ptr = pg.add_coverage(read_ptr, node_id, prg_id, "0");
    uint32_t covg = node_ptr->covg;
    node_ptr = pg.add_coverage(read_ptr, node_id, prg_id, "0");
    EXPECT_EQ(node_ptr->covg - covg, (uint) 1);
}

TEST(PangenomeGraphAddCoverage, NodeAlreadyExists_PangenomeGraphNodeReadsContainsReadTwice) {
    PGraphTester pg;

    uint32_t read_id = 2;
    ReadPtr read_ptr = pg.get_read(read_id);

    std::set<MinimizerHitPtr, pComp> mhs;
    uint32_t node_id = 0;
    uint32_t prg_id = 1;

    pg.add_coverage(read_ptr, node_id, prg_id, "0");
    NodePtr node_ptr = pg.add_coverage(read_ptr, node_id, prg_id, "0");

    auto result = node_ptr->reads.count(read_ptr);
    uint expected = 2;
    EXPECT_EQ(result, expected);
}

TEST(PangenomeGraphAddNode, AddClusterWrongReadId_AssertCatches) {
    uint32_t read_id = 1;
    uint32_t not_read_id = 7;

    std::set<MinimizerHitPtr, pComp> cluster;
    uint32_t prg_id = 4;
    Interval interval(0, 5);
    std::deque<Interval> raw_path = {Interval(7, 8), Interval(10, 14)};
    prg::Path path;
    path.initialize(raw_path);
    MinimizerHitPtr minimizer_hit(std::make_shared<MinimizerHit>(not_read_id, interval, prg_id, path, 0, 0));
    cluster.insert(minimizer_hit);

    PGraphTester pg;

    EXPECT_DEATH(pg.add_node(prg_id, "", read_id, cluster), "");
}

TEST(PangenomeGraphAddNode, AddClusterWrongPrgId_AssertCatches) {
    uint32_t read_id = 1;
    Read read(read_id);

    std::set<MinimizerHitPtr, pComp> cluster;
    uint32_t prg_id = 4;
    uint32_t not_prg_id = 7;
    Interval interval(0, 5);
    std::deque<Interval> raw_path = {Interval(7, 8), Interval(10, 14)};
    prg::Path path;
    path.initialize(raw_path);
    MinimizerHitPtr minimizer_hit(std::make_shared<MinimizerHit>(read_id, interval, not_prg_id, path, 0, 0));
    cluster.insert(minimizer_hit);

    PGraphTester pg;

    EXPECT_DEATH(pg.add_node(prg_id, "", read_id, cluster), "");
}

TEST(PangenomeGraphAddNode, AddNode_PangenomeGraphNodesContainsNodeId) {
    std::set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg;
    uint32_t node_id = 0;
    uint32_t read_id = 1;
    pg.add_node(node_id, "0", read_id, mhs);

    auto result = pg.nodes.find(node_id) != pg.nodes.end();
    EXPECT_TRUE(result);
}

TEST(PangenomeGraphAddNode, AddNode_PangenomeGraphNodeHasRightProperties) {
    std::set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg;
    uint32_t node_id = 0;
    uint32_t read_id = 1;
    pg.add_node(node_id, "0", read_id, mhs);

    NodePtr pan_node = std::make_shared<pangenome::Node>(node_id, node_id, "0");
    EXPECT_EQ(*pg.nodes[0], *pan_node);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->name, "0");
    EXPECT_EQ(pg.nodes[0]->covg, (uint) 1);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint) 1);
}

TEST(PangenomeGraphAddNode, AddNode_PangenomeGraphReadHasRightProperties) {
    std::set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg;
    uint32_t node_id = 0;
    uint32_t read_id = 1;
    pg.add_node(node_id, "0", read_id, mhs);

    ReadPtr pr = std::make_shared<Read>(1);
    EXPECT_EQ(pg.reads.size(), (uint) 1);
    EXPECT_EQ(*pg.reads[1], *pr);
    EXPECT_EQ(pg.reads[1]->hits.size(), (uint) 1);
    EXPECT_EQ(pg.reads[1]->hits[0].size(), (uint) 0);
}

TEST(PangenomeGraphTest, add_node_sample) {
    // add node and check it's there
    PGraphTester pg;

    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA"));
    std::vector<KmerNodePtr> kmp;

    pg.add_node(0, "zero", "sample", 0, l0, kmp);

    EXPECT_EQ(pg.nodes.size(), (uint) 1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint) 1);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint) 0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint) 1);

    EXPECT_EQ(pg.samples.size(), (uint) 1);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint) 1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint) 1);

    EXPECT_EQ(pg.reads.size(), (uint) 0);

    // add a second time
    pg.add_node(0, "zero", "sample", 0, l0, kmp);
    EXPECT_EQ(pg.nodes.size(), (uint) 1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint) 2);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint) 0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint) 1);

    EXPECT_EQ(pg.samples.size(), (uint) 1);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint) 1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint) 2);

    EXPECT_EQ(pg.reads.size(), (uint) 0);

    // add a node with a different sample
    pg.add_node(0, "zero", "sample1", 1, l0, kmp);
    EXPECT_EQ(pg.nodes.size(), (uint) 1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint) 3);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint) 0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint) 2);

    EXPECT_EQ(pg.samples.size(), (uint) 2);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint) 1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint) 2);
    EXPECT_EQ(pg.samples["sample1"]->name, "sample1");
    EXPECT_EQ(pg.samples["sample1"]->paths.size(), (uint) 1);
    EXPECT_EQ(pg.samples["sample1"]->paths[0].size(), (uint) 1);

    EXPECT_EQ(pg.reads.size(), (uint) 0);

    // add a node with a different prg
    pg.add_node(1, "one", "sample1", 1, l0, kmp);
    EXPECT_EQ(pg.nodes.size(), (uint) 2);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint) 3);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint) 0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint) 2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint) 1);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint) 1);
    EXPECT_EQ(pg.nodes[1]->name, "one");
    EXPECT_EQ(pg.nodes[1]->covg, (uint) 1);
    EXPECT_EQ(pg.nodes[1]->reads.size(), (uint) 0);
    EXPECT_EQ(pg.nodes[1]->samples.size(), (uint) 1);

    EXPECT_EQ(pg.samples.size(), (uint) 2);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint) 1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint) 2);
    EXPECT_EQ(pg.samples["sample1"]->name, "sample1");
    EXPECT_EQ(pg.samples["sample1"]->paths.size(), (uint) 2);
    EXPECT_EQ(pg.samples["sample1"]->paths[0].size(), (uint) 1);
    EXPECT_EQ(pg.samples["sample1"]->paths[1].size(), (uint) 1);

    EXPECT_EQ(pg.reads.size(), (uint) 0);
}

TEST(PangenomeGraphTest, clear) {
    // read pg
    std::set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg;
    pg.add_node(0, "0", 1, mhs);
    EXPECT_EQ(pg.nodes.size(), (uint) 1);
    EXPECT_EQ(pg.reads.size(), (uint) 1);
    EXPECT_EQ(pg.samples.size(), (uint) 0);
    pg.clear();
    EXPECT_EQ(pg.nodes.size(), (uint) 0);
    EXPECT_EQ(pg.reads.size(), (uint) 0);
    EXPECT_EQ(pg.samples.size(), (uint) 0);

    // sample pg
    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA"));
    std::vector<KmerNodePtr> kmp;
    pg.add_node(0, "zero", "sample", 0, l0, kmp);
    EXPECT_EQ(pg.reads.size(), (uint) 0);
    EXPECT_EQ(pg.samples.size(), (uint) 1);
    pg.clear();
    EXPECT_EQ(pg.nodes.size(), (uint) 0);
    EXPECT_EQ(pg.reads.size(), (uint) 0);
    EXPECT_EQ(pg.samples.size(), (uint) 0);
}


TEST(PangenomeGraphTest, equals) {
    std::set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg1;
    pg1.add_node(0, "0", 0, mhs);
    pg1.add_node(1, "1", 2, mhs);
    pg1.add_node(1, "1", 0, mhs);
    pg1.add_node(2, "2", 2, mhs);

    PGraphTester pg2;
    pg2.add_node(1, "1", 2, mhs);
    pg2.add_node(0, "0", 0, mhs);
    pg2.add_node(2, "2", 2, mhs);
    pg2.add_node(1, "1", 0, mhs);

    // adding nodes in different order should make no difference
    EXPECT_EQ(pg1, pg1);
    EXPECT_EQ(pg2, pg2);
    EXPECT_EQ(pg1, pg2);
    EXPECT_EQ(pg2, pg1);

    // should not matter if node_id is different provided prg_id is same
    pg2.nodes[7] = std::make_shared<pangenome::Node>(2, 7, "2");
    pg2.nodes.erase(2);
    EXPECT_EQ(pg2, pg2);
    EXPECT_EQ(pg1, pg2);
    EXPECT_EQ(pg2, pg1);

    // or one extra node
    pg2.add_node(3, "3", 0, mhs);
    EXPECT_EQ((pg1 == pg2), false);
    EXPECT_EQ((pg2 == pg1), false);

    // should not break when have a cycle in pangraph
    pg1.add_node(0, "0", 0, mhs);
    EXPECT_EQ(pg1, pg1);
}

TEST(PangenomeGraphTest, not_equals) {
    std::set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg1;
    pg1.add_node(0, "0", 0, mhs);
    pg1.add_node(1, "1", 2, mhs);
    pg1.add_node(1, "1", 0, mhs);
    pg1.add_node(2, "2", 2, mhs);

    PGraphTester pg2;
    pg2.add_node(1, "1", 2, mhs);
    pg2.add_node(0, "0", 0, mhs);
    pg2.add_node(2, "2", 2, mhs);
    pg2.add_node(1, "1", 0, mhs);

    // adding nodes in different order should make no difference
    EXPECT_EQ((pg1 != pg1), false);
    EXPECT_EQ((pg2 != pg2), false);
    EXPECT_EQ((pg1 != pg2), false);
    EXPECT_EQ((pg2 != pg1), false);

    // or one extra node
    pg2.add_node(3, "3", 0, mhs);
    EXPECT_EQ((pg1 != pg2), true);
    EXPECT_EQ((pg2 != pg1), true);

    // should not break when have a cycle in pangraph
    pg1.add_node(0, "0", 0, mhs);
    EXPECT_EQ((pg1 != pg1), false);
}

TEST(PangenomeGraphTest, remove_node) {
    std::set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(0, "0", 0, mhs);
    pg1.add_node(1, "1", 0, mhs);
    pg1.add_node(2, "2", 0, mhs);
    pg1.add_node(3, "3", 0, mhs);

    // read 0: 0->1->3
    pg2.add_node(0, "0", 0, mhs);
    pg2.add_node(1, "1", 0, mhs);
    pg2.add_node(3, "3", 0, mhs);

    pg1.remove_node(pg1.nodes[2]);
    EXPECT_EQ(pg1, pg2);
}

TEST(PangenomeGraphTest, remove_read) {
    std::set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg1, pg2, pg3;
    // read 0: 0->1->2->3
    pg1.add_node(0, "0", 0, mhs);
    pg1.add_node(1, "1", 0, mhs);
    pg1.add_node(2, "2", 0, mhs);
    pg1.add_node(3, "3", 0, mhs);

    // read 1: 4->5->0->5
    pg1.add_node(4, "4", 1, mhs);
    pg1.add_node(5, "5", 1, mhs);
    pg1.add_node(0, "0", 1, mhs);
    pg1.add_node(5, "5", 1, mhs);

    // read 1: 4->5->0->5
    pg2.add_node(4, "0", 1, mhs);
    pg2.add_node(5, "5", 1, mhs);
    pg2.add_node(0, "0", 1, mhs);
    pg2.add_node(5, "5", 1, mhs);

    pg1.remove_read(0);
    EXPECT_EQ(pg1, pg2);

    EXPECT_EQ(pg1.nodes[4]->covg, pg2.nodes[4]->covg);
    EXPECT_EQ(pg1.nodes[5]->covg, pg2.nodes[5]->covg);
    EXPECT_EQ(pg1.nodes[0]->covg, pg2.nodes[0]->covg);
    EXPECT_EQ(pg1.nodes[4]->reads.size(), pg2.nodes[4]->reads.size());
    EXPECT_EQ(pg1.nodes[5]->reads.size(), pg2.nodes[5]->reads.size());
    EXPECT_EQ(pg1.nodes[0]->reads.size(), pg2.nodes[0]->reads.size());

    pg1.remove_read(1);
    EXPECT_EQ(pg1, pg3);
}

TEST(PangenomeGraphTest, remove_low_covg_nodes) {
    std::set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg1, pg2, pg3;
    // read 0: 0->1->2->3
    pg1.add_node(0, "0", 0, mhs);
    pg1.add_node(1, "1", 0, mhs);
    pg1.add_node(2, "2", 0, mhs);
    pg1.add_node(3, "3", 0, mhs);
    // read 1: -4 -> -3 -> -1
    pg1.add_node(1, "1", 1, mhs);
    pg1.add_node(3, "3", 1, mhs);
    pg1.add_node(4, "4", 1, mhs);
    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(0, "0", 2, mhs);
    pg1.add_node(1, "1", 2, mhs);
    pg1.add_node(3, "3", 2, mhs);
    pg1.add_node(4, "4", 2, mhs);
    // read 3: 0 -> 5
    pg1.add_node(0, "0", 3, mhs);
    pg1.add_node(5, "5", 3, mhs);
    // read 4: 5 -> 1
    pg1.add_node(5, "5", 4, mhs);
    pg1.add_node(1, "1", 4, mhs);

    // read 0: 0->1->3
    pg2.add_node(0, "0", 0, mhs);
    pg2.add_node(1, "1", 0, mhs);
    pg2.add_node(3, "3", 0, mhs);
    // read 1: -4 -> -3 -> -1
    pg2.add_node(1, "1", 1, mhs);
    pg2.add_node(3, "3", 1, mhs);
    pg2.add_node(4, "4", 1, mhs);
    // read 2: 0 -> 1 -> 3 -> 4
    pg2.add_node(0, "0", 2, mhs);
    pg2.add_node(1, "1", 2, mhs);
    pg2.add_node(3, "3", 2, mhs);
    pg2.add_node(4, "4", 2, mhs);
    // read 3: 0 -> 5
    pg2.add_node(0, "0", 3, mhs);
    pg2.add_node(5, "5", 3, mhs);
    // read 4: 5 -> 1
    pg2.add_node(5, "5", 4, mhs);
    pg2.add_node(1, "1", 4, mhs);

    pg1.remove_low_covg_nodes(1);
    EXPECT_EQ(pg1, pg2);

    // read 0: 0->1->3
    pg3.add_node(0, "0", 0, mhs);
    pg3.add_node(1, "1", 0, mhs);
    pg3.add_node(3, "3", 0, mhs);
    // read 1: -4 -> -3 -> -1
    pg3.add_node(1, "1", 1, mhs);
    pg3.add_node(3, "3", 1, mhs);
    // read 2: 0 -> 1 -> 3 -> 4
    pg3.add_node(0, "0", 2, mhs);
    pg3.add_node(1, "1", 2, mhs);
    pg3.add_node(3, "3", 2, mhs);
    // read 3: 0 -> 5
    pg3.add_node(0, "0", 3, mhs);
    // read 4: 5 -> 1
    pg3.add_node(1, "1", 4, mhs);

    pg1.remove_low_covg_nodes(2);
    EXPECT_EQ(pg1, pg3);
}

TEST(PangenomeGraphTest, split_node_by_reads) {
    std::set<MinimizerHitPtr, pComp> mhs;

    PGraphTester pg1, pg2, pg3;
    // read 0: 0->1->2->3
    pg1.add_node(0, "0", 0, mhs);
    pg1.add_node(1, "1", 0, mhs);
    pg1.add_node(2, "2", 0, mhs);
    pg1.add_node(3, "3", 0, mhs);

    // read 1: 4->5->0->5
    pg1.add_node(4, "4", 1, mhs);
    pg1.add_node(5, "5", 1, mhs);
    pg1.add_node(0, "0", 1, mhs);
    pg1.add_node(5, "5", 1, mhs);

    EXPECT_EQ((uint) 6, pg1.nodes.size());
    EXPECT_EQ(pg1.nodes[0]->prg_id, (uint) 0);
    EXPECT_EQ(pg1.nodes[0]->covg, (uint) 2);
    EXPECT_EQ(pg1.nodes[1]->prg_id, (uint) 1);
    EXPECT_EQ(pg1.nodes[1]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[2]->prg_id, (uint) 2);
    EXPECT_EQ(pg1.nodes[2]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[3]->prg_id, (uint) 3);
    EXPECT_EQ(pg1.nodes[3]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[4]->prg_id, (uint) 4);
    EXPECT_EQ(pg1.nodes[4]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[5]->prg_id, (uint) 5);
    EXPECT_EQ(pg1.nodes[5]->covg, (uint) 2);

    // read 0: 0->1->2->3
    pg2.add_node(0, "0", 0, mhs);
    pg2.add_node(1, "1", 0, mhs);
    NodePtr n = std::make_shared<pangenome::Node>(2, 7, "2");
    pg2.nodes[7] = n;
    pg2.add_node(3, "3", 0, mhs);

    // read 1: 4->5->0->5
    pg2.add_node(4, "4", 1, mhs);
    pg2.add_node(5, "5", 1, mhs);
    pg2.add_node(0, "0", 1, mhs);
    pg2.add_node(5, "5", 1, mhs);

    std::unordered_set<ReadPtr> reads = {pg1.reads[0]};
    std::vector<uint_least32_t> node_ids = {1, 2, 3};
    std::vector<uint_least32_t> node_ids_exp = {1, 6, 3};
    std::vector<bool> node_orients = {0, 0, 0};
    pg1.split_node_by_reads(reads, node_ids, node_orients, 2);
    EXPECT_EQ(pg1, pg2);
    EXPECT_ITERABLE_EQ(std::vector<uint_least32_t>, node_ids_exp, node_ids);

    EXPECT_EQ((uint) 6, pg1.nodes.size());
    EXPECT_EQ(pg1.nodes[0]->prg_id, (uint) 0);
    EXPECT_EQ(pg1.nodes[0]->covg, (uint) 2);
    EXPECT_EQ(pg1.nodes[1]->prg_id, (uint) 1);
    EXPECT_EQ(pg1.nodes[1]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[6]->prg_id, (uint) 2);
    EXPECT_EQ(pg1.nodes[6]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[3]->prg_id, (uint) 3);
    EXPECT_EQ(pg1.nodes[3]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[4]->prg_id, (uint) 4);
    EXPECT_EQ(pg1.nodes[4]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[5]->prg_id, (uint) 5);
    EXPECT_EQ(pg1.nodes[5]->covg, (uint) 2);

    // read 0: 0->1->2->3
    pg3.add_node(0, "0", 0, mhs);
    pg3.add_node(1, "1", 0, mhs);
    n = std::make_shared<pangenome::Node>(2, 7, "2");
    pg3.nodes[7] = n;
    pg3.add_node(3, "3", 0, mhs);

    // read 1: 4->5->0->5
    pg3.add_node(4, "4", 1, mhs);
    n = std::make_shared<pangenome::Node>(5, 8, "5");
    pg3.nodes[8] = n;
    pg3.add_node(0, "0", 1, mhs);
    pg3.add_node(5, "5", 1, mhs);

    reads = {pg1.reads[1]};
    node_ids = {5, 0, 5};
    node_ids_exp = {7, 0, 5};
    pg1.split_node_by_reads(reads, node_ids, node_orients, 5);
    EXPECT_EQ(pg1, pg3);
    EXPECT_ITERABLE_EQ(std::vector<uint_least32_t>, node_ids_exp, node_ids);

    EXPECT_EQ((uint) 7, pg1.nodes.size());
    EXPECT_EQ(pg1.nodes[0]->prg_id, (uint) 0);
    EXPECT_EQ(pg1.nodes[0]->covg, (uint) 2);
    EXPECT_EQ(pg1.nodes[1]->prg_id, (uint) 1);
    EXPECT_EQ(pg1.nodes[1]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[6]->prg_id, (uint) 2);
    EXPECT_EQ(pg1.nodes[6]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[3]->prg_id, (uint) 3);
    EXPECT_EQ(pg1.nodes[3]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[4]->prg_id, (uint) 4);
    EXPECT_EQ(pg1.nodes[4]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[5]->prg_id, (uint) 5);
    EXPECT_EQ(pg1.nodes[5]->covg, (uint) 1);
    EXPECT_EQ(pg1.nodes[7]->prg_id, (uint) 5);
    EXPECT_EQ(pg1.nodes[7]->covg, (uint) 1);

}

TEST(PangenomeGraphTest, add_hits_to_kmergraph) {
}

TEST(PangenomeGraphTest, save_matrix) {
    // add node and check it's there
    PGraphTester pg;

    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA"));
    std::vector<KmerNodePtr> kmp;

    pg.add_node(0, "zero", "sample1", 0, l0, kmp);
    pg.add_node(0, "zero", "sample1", 0, l0, kmp);
    pg.add_node(0, "zero", "sample2", 0, l0, kmp);
    pg.add_node(1, "one", "sample1", 0, l0, kmp);
    pg.add_node(2, "two", "sample3", 0, l0, kmp);

    std::vector<std::string> names = {"sample1", "sample2", "sample3", "sample4"};

    pg.save_matrix("../../test/test_cases/pangraph_test_save.matrix", names);
}

TEST(PangenomeGraphTest, save_mapped_read_strings) {
    PGraphTester pg;
    pangenome::ReadPtr pr;
    MinimizerHits mhits;

    Minimizer m;
    std::deque<Interval> d;
    prg::Path p;
    MiniRecord *mr;

    // read1
    m = Minimizer(0, 1, 6, 0); // kmer, start, end, strand
    d = {Interval(7, 8), Interval(10, 14)};
    p.initialize(d);
    mr = new MiniRecord(0, p, 0, 0);
    mhits.add_hit(1, m, mr); // read 1

    m = Minimizer(0, 0, 5, 0);
    d = {Interval(6, 10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    mhits.add_hit(1, m, mr);

    d = {Interval(6, 10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    mhits.add_hit(1, m, mr);

    mhits.sort();
    pg.add_node(0, "zero", 1, mhits.hits);
    mhits.clear();

    //read 2
    m = Minimizer(0, 2, 7, 1);
    d = {Interval(6, 10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    mhits.add_hit(2, m, mr);

    m = Minimizer(0, 5, 10, 1);
    d = {Interval(6, 10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    mhits.add_hit(2, m, mr);

    mhits.sort();
    delete mr;
    pg.add_node(0, "zero", 2, mhits.hits);

    std::string expected1 = ">read1 pandora: 1 0:6 + \nshould\n>read2 pandora: 2 2:10 - \nis time \n";
    std::string expected2 = ">read2 pandora: 2 2:10 - \nis time \n>read1 pandora: 1 0:6 + \nshould\n";

    pg.save_mapped_read_strings("../../test/test_cases/reads.fa", "save_mapped_read_strings");
    std::ifstream ifs("save_mapped_read_strings/zero/zero.reads.fa");
    std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    EXPECT_TRUE((content == expected1) or (content == expected2));

    pg.save_mapped_read_strings("../../test/test_cases/reads.fa", ".");
    std::ifstream ifs2("zero/zero.reads.fa");
    std::string content2((std::istreambuf_iterator<char>(ifs2)), (std::istreambuf_iterator<char>()));
    EXPECT_TRUE((content2 == expected1) or (content2 == expected2));
}

TEST(PangenomeGraphTest, get_node_closest_vcf_reference_no_paths) {
    uint32_t prg_id = 3, w=1, k=3;
    std::string prg_name = "nested varsite";
    LocalPRG l3(prg_id, prg_name , "A 5 G 7 C 8 T 7  6 G 5 T");
    auto index = std::make_shared<Index>();
    l3.minimizer_sketch(index, w, k);
    auto prg_ptr = std::make_shared<LocalPRG>(l3);

    pangenome::Graph pangraph;
    std::string sample_name = "null_test_sample";
    std::vector<KmerNodePtr> sample_kmer_path = {};

    pangraph.add_node(prg_id, prg_name, sample_name, 0, prg_ptr, sample_kmer_path);
    auto path = pangraph.get_node_closest_vcf_reference(*pangraph.nodes[prg_id], w, l3);
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, path, l3.prg.top_path());
}

TEST(PangenomeGraphTest, get_node_closest_vcf_reference_one_path) {
    uint32_t prg_id = 3, w=1, k=3;
    std::string prg_name = "nested varsite";
    LocalPRG l3(prg_id, prg_name , "A 5 G 7 C 8 T 7  6 G 5 T");
    auto index = std::make_shared<Index>();
    l3.minimizer_sketch(index, w, k);
    auto prg_ptr = std::make_shared<LocalPRG>(l3);

    pangenome::Graph pangraph;
    std::string sample_name = "single_test_sample";

    auto &kg = l3.kmer_prg;
    std::vector<KmerNodePtr> sample_kmer_path = {kg.nodes[0], kg.nodes[2], kg.nodes[5], kg.nodes[6]};

    pangraph.add_node(prg_id, prg_name, sample_name, 0, prg_ptr, sample_kmer_path);
    auto &node = *pangraph.nodes[prg_id];

    auto path = pangraph.get_node_closest_vcf_reference(node, w, l3);
    std::vector<LocalNodePtr> exp_path = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};

    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, path, exp_path);
}

TEST(PangenomeGraphTest, get_node_closest_vcf_reference_three_paths) {
    uint32_t prg_id = 3, w=1, k=3;
    std::string prg_name = "nested varsite";
    LocalPRG l3(prg_id, prg_name , "A 5 G 7 C 8 T 7  6 G 5 T");
    auto index = std::make_shared<Index>();
    l3.minimizer_sketch(index, w, k);
    auto prg_ptr = std::make_shared<LocalPRG>(l3);

    pangenome::Graph pangraph;
    auto &kg = l3.kmer_prg;

    std::string sample_name = "test_sample1";
    std::vector<KmerNodePtr> sample_kmer_path = {kg.nodes[0], kg.nodes[2], kg.nodes[5], kg.nodes[6]};
    pangraph.add_node(prg_id, prg_name, sample_name, 0, prg_ptr, sample_kmer_path);

    sample_name = "test_sample1_again";
    sample_kmer_path = {kg.nodes[0], kg.nodes[2], kg.nodes[5], kg.nodes[6]};
    pangraph.add_node(prg_id, prg_name, sample_name, 1, prg_ptr, sample_kmer_path);

    sample_name = "test_sample2";
    sample_kmer_path = {kg.nodes[0], kg.nodes[4], kg.nodes[6]};
    pangraph.add_node(prg_id, prg_name, sample_name, 2, prg_ptr, sample_kmer_path);

    auto &node = *pangraph.nodes[prg_id];
    auto path = pangraph.get_node_closest_vcf_reference(node, w, l3);
    std::vector<LocalNodePtr> exp_path = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[3], l3.prg.nodes[4], l3.prg.nodes[6]};

    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, path, exp_path);
}

TEST(PangenomeGraphTest, copy_coverages_to_kmergraphs){
    uint32_t prg_id = 3, w=1, k=3;
    std::string prg_name = "nested varsite", sample_name = "sample";
    LocalPRG l3(prg_id, prg_name , "A 5 G 7 C 8 T 7  6 G 5 T");
    auto index = std::make_shared<Index>();
    l3.minimizer_sketch(index, w, k);
    auto prg_ptr = std::make_shared<LocalPRG>(l3);

    pangenome::Graph ref_pangraph;
    auto sample_id = 0;
    std::vector<KmerNodePtr> empty = {};
    ref_pangraph.add_node(prg_id, prg_name, sample_name, sample_id, prg_ptr, empty);

    EXPECT_TRUE(ref_pangraph.nodes.find(prg_id) != ref_pangraph.nodes.end());
    ref_pangraph.nodes[prg_id]->kmer_prg = l3.kmer_prg;
    auto &kg = ref_pangraph.nodes[prg_id]->kmer_prg;
    EXPECT_EQ(kg.nodes.size(), (uint)7);
    kg.nodes[2]->set_covg(5, 1, sample_id);
    kg.nodes[4]->set_covg(8, 0, sample_id);
    kg.nodes[5]->set_covg(2, 1, sample_id);
    kg.nodes[6]->set_covg(5, 0, sample_id);

    pangenome::Graph pangraph;
    sample_id = 3;
    pangraph.add_node(prg_id, prg_name, sample_name, sample_id, prg_ptr, empty);

    LocalPRG dummy(0,"null","");
    auto dummy_prg_ptr = std::make_shared<LocalPRG>(dummy);
    pangraph.setup_kmergraphs({dummy_prg_ptr, dummy_prg_ptr, dummy_prg_ptr, prg_ptr}, 4);

    pangraph.copy_coverages_to_kmergraphs(ref_pangraph, sample_id);

    for (uint32_t id = 0; id < 3; ++id){
        for (const auto &node : pangraph.nodes[prg_id]->kmer_prg.nodes){
            EXPECT_EQ(node->get_covg(0,id), (uint)0);
            EXPECT_EQ(node->get_covg(1,id), (uint)0);
        }
    }
    auto id = sample_id;
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[0]->get_covg(0,id), (uint)0);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[0]->get_covg(1,id), (uint)0);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[1]->get_covg(0,id), (uint)0);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[1]->get_covg(1,id), (uint)0);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[2]->get_covg(0,id), (uint)0);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[2]->get_covg(1,id), (uint)5);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[3]->get_covg(0,id), (uint)0);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[3]->get_covg(1,id), (uint)0);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[4]->get_covg(0,id), (uint)8);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[4]->get_covg(1,id), (uint)0);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[5]->get_covg(0,id), (uint)0);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[5]->get_covg(1,id), (uint)2);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[6]->get_covg(0,id), (uint)5);
    EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg.nodes[6]->get_covg(1,id), (uint)0);
}

TEST(PangenomeGraphTest, infer_node_vcf_reference_path_no_file_strings)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    std::vector<std::string> prg_strings;
    prg_strings.push_back("ATGCCGGTAATTAAAGTACGTGAAAAGAAACTGGCTC 5 A 6 G 5 CGAAAACGCACGCCGCACTCGTCTGTAC");
    prg_strings.push_back("A 5 G 7 C 8 T 7  6 G 5 T");
    prg_strings.push_back("TC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AG");
    prg_strings.push_back("A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");

    pangenome::Graph pangraph;
    std::vector<KmerNodePtr> empty;
    auto index = std::make_shared<Index>();
    uint32_t prg_id = 0, sample_id = 0, w = 1, k = 3;
    std::unordered_map<std::string, std::string> vcf_refs;
    std::vector<std::vector<LocalNodePtr>> vcf_ref_paths;
    for (const auto &prg_string : prg_strings){
        std::string prg_name = "prg" + std::to_string(prg_id), sample_name = "sample";
        prgs.emplace_back(std::make_shared<LocalPRG>(LocalPRG(prg_id, prg_name, prg_string)));
        prgs.back()->minimizer_sketch(index, w, k);
        pangraph.add_node(prg_id, prg_name, sample_name, sample_id, prgs.back(), empty);
        vcf_ref_paths.emplace_back(pangraph.infer_node_vcf_reference_path(*pangraph.nodes[prg_id], prgs.back(), w, vcf_refs));
        prg_id++;
    }

    EXPECT_EQ(vcf_ref_paths.size(), (uint)4);
    for (uint j=0; j<vcf_ref_paths.size(); ++j)
        EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, vcf_ref_paths[j], prgs[j]->prg.top_path());
}

TEST(PangenomeGraphTest, infer_node_vcf_reference_path_with_file_strings)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    std::vector<std::string> prg_strings;
    prg_strings.push_back("ATGCCGGTAATTAAAGTACGTGAAAAGAAACTGGCTC 5 A 6 G 5 CGAAAACGCACGCCGCACTCGTCTGTAC");
    prg_strings.push_back("A 5 G 7 C 8 T 7  6 G 5 T");
    prg_strings.push_back("TC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AG");
    prg_strings.push_back("AATTTTTTTGGGGTTGGTTTTAAA 5 GGGGG 7 CCCCCC 8 TTTTTT 7 TTTTTT 9 CCGCCGCCGCCG 10 CGGCCGCCG 9  6 GGGGG 5 TATAAAAATTTTTT");
    std::unordered_map<std::string, std::string> vcf_ref_strings;
    vcf_ref_strings["prg0"] = "ATGCCGGTAATTAAAGTACGTGAAAAGAAACTGGCTCGCGAAAACGCACGCCGCACTCGTCTGTAC"; // valid
    vcf_ref_strings["prg1"] = "AGT"; //invalid, too short
    vcf_ref_strings["prg2"] = "ATGCCGGTAATTAAAGTACGTGAAAAGAAACTGGCTCGCGAAAACGCACGCCGCACTCGTCTGTAC"; //invalid, is not a path through prg
    vcf_ref_strings["prg3"] = "AATTTTTTTGGGGTTGGTTTTAAAGGGGGTTTTTTTTTTTTCCGCCGCCGCCGTATAAAAATTTTTT"; //valid

    pangenome::Graph pangraph;
    std::vector<KmerNodePtr> empty;
    auto index = std::make_shared<Index>();
    uint32_t prg_id = 0, sample_id = 0, w = 1, k = 3;
    std::vector<std::vector<LocalNodePtr>> vcf_ref_paths;
    for (const auto &prg_string : prg_strings){
        std::string prg_name = "prg" + std::to_string(prg_id), sample_name = "sample";
        prgs.emplace_back(std::make_shared<LocalPRG>(LocalPRG(prg_id, prg_name, prg_string)));
        prgs.back()->minimizer_sketch(index, w, k);
        pangraph.add_node(prg_id, prg_name, sample_name, sample_id, prgs.back(), empty);
        vcf_ref_paths.emplace_back(pangraph.infer_node_vcf_reference_path(*pangraph.nodes[prg_id], prgs.back(),
                                                                          w, vcf_ref_strings));
        prg_id++;
    }

    std::vector<LocalNodePtr> exp_path0 = {prgs[0]->prg.nodes[0], prgs[0]->prg.nodes[2], prgs[0]->prg.nodes[3]};
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, vcf_ref_paths[0], exp_path0);
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, vcf_ref_paths[1], prgs[1]->prg.top_path());
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, vcf_ref_paths[2], prgs[2]->prg.top_path());
    std::vector<LocalNodePtr> exp_path3 = {prgs[3]->prg.nodes[0], prgs[3]->prg.nodes[1], prgs[3]->prg.nodes[3],
                                           prgs[3]->prg.nodes[4], prgs[3]->prg.nodes[5], prgs[3]->prg.nodes[7],
                                           prgs[3]->prg.nodes[9]};
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, vcf_ref_paths[3], exp_path3);
}