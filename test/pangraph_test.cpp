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

    set<MinimizerHitPtr, pComp> mhs;
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

    set<MinimizerHitPtr, pComp> mhs;
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

    set<MinimizerHitPtr, pComp> mhs;
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

    set<MinimizerHitPtr, pComp> mhs;
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

    set<MinimizerHitPtr, pComp> cluster;
    uint32_t prg_id = 4;
    Interval interval(0, 5);
    deque<Interval> raw_path = {Interval(7, 8), Interval(10, 14)};
    Path path;
    path.initialize(raw_path);
    MinimizerHitPtr minimizer_hit(make_shared<MinimizerHit>(not_read_id, interval, prg_id, path, 0, 0));
    cluster.insert(minimizer_hit);

    PGraphTester pg;

    EXPECT_DEATH(pg.add_node(prg_id, "", read_id, cluster), "");
}

TEST(PangenomeGraphAddNode, AddClusterWrongPrgId_AssertCatches) {
    uint32_t read_id = 1;
    Read read(read_id);

    set<MinimizerHitPtr, pComp> cluster;
    uint32_t prg_id = 4;
    uint32_t not_prg_id = 7;
    Interval interval(0, 5);
    deque<Interval> raw_path = {Interval(7, 8), Interval(10, 14)};
    Path path;
    path.initialize(raw_path);
    MinimizerHitPtr minimizer_hit(make_shared<MinimizerHit>(read_id, interval, not_prg_id, path, 0, 0));
    cluster.insert(minimizer_hit);

    PGraphTester pg;

    EXPECT_DEATH(pg.add_node(prg_id, "", read_id, cluster), "");
}

TEST(PangenomeGraphAddNode, AddNode_PangenomeGraphNodesContainsNodeId) {
    set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg;
    uint32_t node_id = 0;
    uint32_t read_id = 1;
    pg.add_node(node_id, "0", read_id, mhs);

    auto result = pg.nodes.find(node_id) != pg.nodes.end();
    EXPECT_TRUE(result);
}

TEST(PangenomeGraphAddNode, AddNode_PangenomeGraphNodeHasRightProperties) {
    set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg;
    uint32_t node_id = 0;
    uint32_t read_id = 1;
    pg.add_node(node_id, "0", read_id, mhs);

    NodePtr pn = make_shared<Node>(node_id, node_id, "0");
    EXPECT_EQ(*pg.nodes[0], *pn);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint) 0);
    EXPECT_EQ(pg.nodes[0]->name, "0");
    EXPECT_EQ(pg.nodes[0]->covg, (uint) 1);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint) 1);
}

TEST(PangenomeGraphAddNode, AddNode_PangenomeGraphReadHasRightProperties) {
    set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg;
    uint32_t node_id = 0;
    uint32_t read_id = 1;
    pg.add_node(node_id, "0", read_id, mhs);

    ReadPtr pr = make_shared<Read>(1);
    EXPECT_EQ(pg.reads.size(), (uint) 1);
    EXPECT_EQ(*pg.reads[1], *pr);
    EXPECT_EQ(pg.reads[1]->hits.size(), (uint) 1);
    EXPECT_EQ(pg.reads[1]->hits[0].size(), (uint) 0);
}

TEST(PangenomeGraphTest, add_node_sample) {
    // add node and check it's there
    PGraphTester pg;

    LocalPRG *l0;
    l0 = new LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA");
    vector<KmerNodePtr> kmp;

    pg.add_node(0, "zero", "sample", kmp, l0);

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
    pg.add_node(0, "zero", "sample", kmp, l0);
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
    pg.add_node(0, "zero", "sample1", kmp, l0);
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
    pg.add_node(1, "one", "sample1", kmp, l0);
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

    delete l0;
}

TEST(PangenomeGraphTest, clear) {
    // read pg
    set<MinimizerHitPtr, pComp> mhs;

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
    LocalPRG *l0;
    l0 = new LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA");
    vector<KmerNodePtr> kmp;
    pg.add_node(0, "zero", "sample", kmp, l0);
    EXPECT_EQ(pg.reads.size(), (uint) 0);
    EXPECT_EQ(pg.samples.size(), (uint) 1);
    pg.clear();
    EXPECT_EQ(pg.nodes.size(), (uint) 0);
    EXPECT_EQ(pg.reads.size(), (uint) 0);
    EXPECT_EQ(pg.samples.size(), (uint) 0);
    delete l0;
}


TEST(PangenomeGraphTest, equals) {
    set<MinimizerHitPtr, pComp> mhs;
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
    pg2.nodes[7] = make_shared<Node>(2, 7, "2");
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
    set<MinimizerHitPtr, pComp> mhs;
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
    set<MinimizerHitPtr, pComp> mhs;

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
    set<MinimizerHitPtr, pComp> mhs;

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
    set<MinimizerHitPtr, pComp> mhs;

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
    set<MinimizerHitPtr, pComp> mhs;

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
    NodePtr n = make_shared<Node>(2, 7, "2");
    pg2.nodes[7] = n;
    pg2.add_node(3, "3", 0, mhs);

    // read 1: 4->5->0->5
    pg2.add_node(4, "4", 1, mhs);
    pg2.add_node(5, "5", 1, mhs);
    pg2.add_node(0, "0", 1, mhs);
    pg2.add_node(5, "5", 1, mhs);

    unordered_set<ReadPtr> reads = {pg1.reads[0]};
    vector<uint_least32_t> node_ids = {1, 2, 3};
    vector<uint_least32_t> node_ids_exp = {1, 6, 3};
    vector<bool> node_orients = {0, 0, 0};
    pg1.split_node_by_reads(reads, node_ids, node_orients, 2);
    EXPECT_EQ(pg1, pg2);
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, node_ids_exp, node_ids);

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
    n = make_shared<Node>(2, 7, "2");
    pg3.nodes[7] = n;
    pg3.add_node(3, "3", 0, mhs);

    // read 1: 4->5->0->5
    pg3.add_node(4, "4", 1, mhs);
    n = make_shared<Node>(5, 8, "5");
    pg3.nodes[8] = n;
    pg3.add_node(0, "0", 1, mhs);
    pg3.add_node(5, "5", 1, mhs);

    reads = {pg1.reads[1]};
    node_ids = {5, 0, 5};
    node_ids_exp = {7, 0, 5};
    pg1.split_node_by_reads(reads, node_ids, node_orients, 5);
    EXPECT_EQ(pg1, pg3);
    EXPECT_ITERABLE_EQ(vector<uint_least32_t>, node_ids_exp, node_ids);

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

    LocalPRG *l0;
    l0 = new LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA");
    vector<KmerNodePtr> kmp;

    pg.add_node(0, "zero", "sample1", kmp, l0);
    pg.add_node(0, "zero", "sample1", kmp, l0);
    pg.add_node(0, "zero", "sample2", kmp, l0);
    pg.add_node(1, "one", "sample1", kmp, l0);
    pg.add_node(2, "two", "sample3", kmp, l0);

    pg.save_matrix("../../test/test_cases/pangraph_test_save.matrix");
}

TEST(PangenomeGraphTest, save_mapped_read_strings) {
    PGraphTester pg;
    pangenome::ReadPtr pr;
    MinimizerHits mhits;

    Minimizer m;
    deque<Interval> d;
    Path p;
    MiniRecord* mr;

    // read1
    m = Minimizer(0,1,6,0); // kmer, start, end, strand
    d = {Interval(7,8), Interval(10, 14)};
    p.initialize(d);
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr); // read 1

    m = Minimizer(0,0,5,0);
    d = {Interval(6,10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr);

    d = {Interval(6,10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr);

    mhits.sort();
    pg.add_node(0,"zero", 1, mhits.hits);
    mhits.clear();

    //read 2
    m = Minimizer(0,2,7,1);
    d = {Interval(6,10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(2, m, mr);

    m = Minimizer(0,5,10,1);
    d = {Interval(6,10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(2, m, mr);

    mhits.sort();
    delete mr;
    pg.add_node(0,"zero", 2, mhits.hits);

    string expected1 = ">read1 pandora: 1 0:6 + \nshould\n>read2 pandora: 2 2:10 - \nis time \n";
    string expected2 = ">read2 pandora: 2 2:10 - \nis time \n>read1 pandora: 1 0:6 + \nshould\n";

    pg.save_mapped_read_strings("../../test/test_cases/reads.fa", "save_mapped_read_strings");
    ifstream ifs("save_mapped_read_strings/zero/zero.reads.fa");
    string content( (std::istreambuf_iterator<char>(ifs) ),(std::istreambuf_iterator<char>()) );
    EXPECT_TRUE((content == expected1) or (content == expected2));

    pg.save_mapped_read_strings("../../test/test_cases/reads.fa", ".");
    ifstream ifs2("zero/zero.reads.fa");
    string content2( (std::istreambuf_iterator<char>(ifs2) ),(std::istreambuf_iterator<char>()) );
    EXPECT_TRUE((content2 == expected1) or (content2 == expected2));
}

TEST(PangenomeGraphTest, save_kmergraph_coverages) {
    boost::filesystem::remove_all("coverages");

    PGraphTester pg;
    MinimizerHits mhits;
    pg.add_node(0,"zero", 2, mhits.hits); //node 0
    pg.add_node(1,"one", 2, mhits.hits);  //node 1
    EXPECT_EQ((uint)2, pg.nodes.size());

    KmerGraph kg;
    deque<Interval> d = {Interval(0, 0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    EXPECT_EQ((uint)7, kg.nodes.size());

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[1], kg.nodes[2]);
    kg.add_edge(kg.nodes[0], kg.nodes[3]);
    kg.add_edge(kg.nodes[3], kg.nodes[4]);
    kg.add_edge(kg.nodes[0], kg.nodes[5]);
    kg.add_edge(kg.nodes[2], kg.nodes[6]);
    kg.add_edge(kg.nodes[4], kg.nodes[6]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);

    pg.nodes[0]->kmer_prg = kg;
    pg.nodes[0]->kmer_prg.nodes[1]->covg[0] += 4;
    pg.nodes[0]->kmer_prg.nodes[2]->covg[0] += 3;
    pg.nodes[1]->kmer_prg = kg;
    pg.nodes[1]->kmer_prg.nodes[5]->covg[1] += 5;

    pg.save_kmergraph_coverages(".", "test_sample");

    string expected1 = "sample\t0\t1\t2\t3\t4\t5\t6\ntest_sample\t0,0\t4,0\t3,0\t0,0\t0,0\t0,0\t0,0\n";
    string expected2 = "sample\t0\t1\t2\t3\t4\t5\t6\ntest_sample\t0,0\t0,0\t0,0\t0,0\t0,0\t0,5\t0,0\n";

    ifstream ifs("coverages/zero.csv");
    string content( (std::istreambuf_iterator<char>(ifs) ),(std::istreambuf_iterator<char>()) );
    EXPECT_EQ(content,expected1);

    ifstream ifs2("coverages/one.csv");
    string content2( (std::istreambuf_iterator<char>(ifs2) ),(std::istreambuf_iterator<char>()) );
    EXPECT_EQ(content2,expected2);
}
