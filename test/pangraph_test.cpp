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

const std::string TEST_CASE_DIR = "../../test/test_cases/";

TEST(PangenomeGraphConstructor, constructors_and_get_sample)
{
    {
        // testing default constructor
        PGraphTester pg;
        std::vector<std::string> expectedSamples = { "sample_1" };
        EXPECT_EQ(pg.samples.size(), 1);
        pg.get_sample("sample_1")->name == "sample_1";
        EXPECT_EQ(pg.next_id, 0);
        EXPECT_EQ(pg.reads.size(), 0);
        EXPECT_EQ(pg.nodes.size(), 0);
    }
    {
        // testing 1-argument constructor
        PGraphTester pg({ "test_sample" });
        std::vector<std::string> expectedSamples = { "test_sample" };
        EXPECT_EQ(pg.samples.size(), 1);
        pg.get_sample("test_sample")->name == "test_sample";
        EXPECT_EQ(pg.next_id, 0);
        EXPECT_EQ(pg.reads.size(), 0);
        EXPECT_EQ(pg.nodes.size(), 0);
    }
    {
        // testing 3-argument constructor
        std::vector<std::string> expectedSamples = { "test1", "test2", "test3" };
        PGraphTester pg(expectedSamples);
        EXPECT_EQ(pg.samples.size(), expectedSamples.size());
        pg.get_sample("test1")->name == "test1";
        pg.get_sample("test2")->name == "test2";
        pg.get_sample("test3")->name == "test3";
        EXPECT_EQ(pg.next_id, 0);
        EXPECT_EQ(pg.reads.size(), 0);
        EXPECT_EQ(pg.nodes.size(), 0);
    }
}

TEST(PangenomeGraphRead, add_read_and_get_read)
{
    PGraphTester pg;
    uint32_t read_id = 2;

    // first trivial tests - read is not in the graph
    EXPECT_EQ(pg.reads.size(), (uint)0);
    try {
        pg.get_read(read_id);
        FAIL() << "pg.get_read(read_id) should throw std::out_of_range\n";
    } catch (const std::out_of_range&) {
    } catch (...) {
        FAIL() << "pg.get_read(read_id) should throw std::out_of_range\n";
    }

    // add read
    pg.add_read(read_id);
    // and test
    EXPECT_EQ(pg.reads.size(), (uint)1);
    EXPECT_EQ(pg.get_read(read_id)->id, read_id);

    // add read again
    pg.add_read(read_id); // this should not change anything
    // and test
    EXPECT_EQ(pg.reads.size(), (uint)1);
    EXPECT_EQ(pg.get_read(read_id)->id, read_id);

    // add a different read
    uint32_t read_id_2 = 10;
    pg.add_read(read_id_2);
    // and test
    EXPECT_EQ(pg.reads.size(), (uint)2);
    EXPECT_EQ(pg.get_read(read_id)->id, read_id);
    EXPECT_EQ(pg.reads.size(), (uint)2);
    EXPECT_EQ(pg.get_read(read_id_2)->id, read_id_2);
}

TEST(PangenomeGraphNode, add_node_and_get_node)
{
    PGraphTester pg;
    uint32_t prg_id = 0;
    auto prg_pointer = std::make_shared<LocalPRG>(prg_id, "", "");

    // trivial test - nodes are empty
    EXPECT_EQ(pg.nodes.size(), (uint)0);
    try {
        pg.get_node(prg_id);
        FAIL() << "pg.get_node() should throw std::out_of_range\n";
    } catch (const std::out_of_range&) {
    } catch (...) {
        FAIL() << "pg.get_node() should throw std::out_of_range\n";
    }
    try {
        pg.get_node(prg_pointer);
        FAIL() << "pg.get_node() should throw std::out_of_range\n";
    } catch (const std::out_of_range&) {
    } catch (...) {
        FAIL() << "pg.get_node() should throw std::out_of_range\n";
    }

    // add a node
    pg.add_node(prg_pointer);
    // test if the node was added
    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.get_node(prg_id)->prg_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_id)->node_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_pointer)->prg_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_pointer)->node_id, prg_id);

    // add the node again - nothing should change
    pg.add_node(prg_pointer);
    // test if the node was added
    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.get_node(prg_id)->prg_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_id)->node_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_pointer)->prg_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_pointer)->node_id, prg_id);

    // add a second node
    uint32_t prg_id_2 = 10;
    auto prg_pointer_2 = std::make_shared<LocalPRG>(prg_id_2, "", "");
    pg.add_node(prg_pointer_2);
    // test if the node was added
    EXPECT_EQ(pg.nodes.size(), (uint)2);
    EXPECT_EQ(pg.get_node(prg_id)->prg_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_id)->node_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_pointer)->prg_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_pointer)->node_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_id_2)->prg_id, prg_id_2);
    EXPECT_EQ(pg.get_node(prg_id_2)->node_id, prg_id_2);
    EXPECT_EQ(pg.get_node(prg_pointer_2)->prg_id, prg_id_2);
    EXPECT_EQ(pg.get_node(prg_pointer_2)->node_id, prg_id_2);

    // add a third node, with a prg that was already added, but with a different node_id
    uint32_t other_node_id_for_prg_pointer = prg_id + 1;
    pg.add_node(prg_pointer, other_node_id_for_prg_pointer);
    // test if it was added
    EXPECT_EQ(pg.nodes.size(), (uint)3);
    EXPECT_EQ(pg.get_node(prg_id)->prg_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_id)->node_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_pointer)->prg_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_pointer)->node_id, prg_id);
    EXPECT_EQ(pg.get_node(prg_id_2)->prg_id, prg_id_2);
    EXPECT_EQ(pg.get_node(prg_id_2)->node_id, prg_id_2);
    EXPECT_EQ(pg.get_node(prg_pointer_2)->prg_id, prg_id_2);
    EXPECT_EQ(pg.get_node(prg_pointer_2)->node_id, prg_id_2);
    EXPECT_EQ(pg.get_node(other_node_id_for_prg_pointer)->prg_id, prg_id);
    EXPECT_EQ(pg.get_node(other_node_id_for_prg_pointer)->node_id,
        other_node_id_for_prg_pointer);
}

// function that will make things easier by setting up MiniRecord, MinimizerHit, cluster
// and prg
void setup_minimizerhit_cluster_prg_function(uint32_t prg_id, uint32_t read_id,
    MiniRecord* mr1, MinimizerHitPtr* minimizer_hit,
    std::shared_ptr<std::set<MinimizerHitPtr, pComp>>* cluster_pointer,
    std::shared_ptr<LocalPRG>* prg_pointer)
{
    std::deque<Interval> raw_path = { Interval(7, 8), Interval(10, 14) };
    prg::Path path;
    path.initialize(raw_path);

    Interval interval(0, 5);
    Minimizer m1(0, interval.start, interval.get_end(), 0); // kmer, start, end, strand

    (*mr1) = MiniRecord(prg_id, path, prg_id, 0);
    (*minimizer_hit) = std::make_shared<MinimizerHit>(read_id, m1, *mr1);

    (*cluster_pointer) = std::make_shared<std::set<MinimizerHitPtr, pComp>>();
    (*cluster_pointer)->insert(*minimizer_hit);

    (*prg_pointer) = std::make_shared<LocalPRG>(prg_id, "", "");
}

TEST(PangenomeGraph_add_hits_between_PRG_and_read, AddClusters)
{
    PGraphTester pg;

    // test 1: add a cluster
    uint32_t prg_id_1 = 1;
    uint32_t read_id_1 = 100;
    MiniRecord mr1;
    MinimizerHitPtr minimizer_hit_1;
    std::shared_ptr<std::set<MinimizerHitPtr, pComp>> cluster_pointer_1;
    std::shared_ptr<LocalPRG> prg_pointer_1;
    {
        setup_minimizerhit_cluster_prg_function(prg_id_1, read_id_1, &mr1,
            &minimizer_hit_1, &cluster_pointer_1, &prg_pointer_1);
        pg.add_hits_between_PRG_and_read(prg_pointer_1, read_id_1, *cluster_pointer_1);

        // test if everything is fine
        EXPECT_EQ(pg.nodes.size(), 1); // is the node really inserted?
        EXPECT_EQ(
            pg.get_node(prg_id_1)->node_id, prg_id_1); // is the node really inserted?
        EXPECT_EQ(pg.get_node(prg_id_1)->covg, 1); // is the coverage of the node 1?
        EXPECT_EQ(pg.get_node(prg_id_1)->reads.size(),
            1); // is the read really inserted in the node?
        EXPECT_EQ(pg.get_node(prg_id_1)->reads.begin()->get()->id,
            read_id_1); // is the read really inserted in the node?
        EXPECT_EQ(pg.reads.size(), 1);
        EXPECT_EQ(pg.get_read(read_id_1)->id, read_id_1); // is the read really
                                                          // inserted?
        EXPECT_EQ(pg.get_read(read_id_1)->get_hits_as_unordered_map().size(),
            1); // is the hit really inserted?
        EXPECT_EQ(*pg.get_read(read_id_1)->get_hits_as_unordered_map()[prg_id_1][0],
            *minimizer_hit_1); // is the hit really inserted?
        EXPECT_EQ(pg.get_read(read_id_1)
                      ->get_hits_as_unordered_map()[prg_id_1][0]
                      ->get_kmer_node_id(),
            prg_id_1); // is the node really inserted in the read?
        EXPECT_EQ(pg.get_read(read_id_1)->node_orientations.size(),
            1); // is the node_orientation was inserted in the read?
        EXPECT_EQ(pg.get_read(read_id_1)->node_orientations[0],
            true); // is the node_orientation was inserted in the read?
    }

    // test 2: add a cluster with the same read_id, but different prg_id
    uint32_t prg_id_2 = 2;
    MiniRecord mr2;
    MinimizerHitPtr minimizer_hit_2;
    std::shared_ptr<std::set<MinimizerHitPtr, pComp>> cluster_pointer_2;
    std::shared_ptr<LocalPRG> prg_pointer_2;
    {
        setup_minimizerhit_cluster_prg_function(prg_id_2, read_id_1, &mr2,
            &minimizer_hit_2, &cluster_pointer_2, &prg_pointer_2);
        pg.add_hits_between_PRG_and_read(prg_pointer_2, read_id_1, *cluster_pointer_2);

        // test if everything is fine
        EXPECT_EQ(pg.nodes.size(), 2); // is the node really inserted?
        EXPECT_EQ(pg.get_node(prg_id_1)->node_id, prg_id_1);
        EXPECT_EQ(
            pg.get_node(prg_id_2)->node_id, prg_id_2); // is the node really inserted?
        EXPECT_EQ(pg.get_node(prg_id_1)->covg, 1); // is the coverage of the node 1?
        EXPECT_EQ(pg.get_node(prg_id_2)->covg, 1); // is the coverage of the node 1?
        EXPECT_EQ(pg.get_node(prg_id_1)->reads.size(),
            1); // is the read really inserted in the node?
        EXPECT_EQ(pg.get_node(prg_id_2)->reads.size(),
            1); // is the read really inserted in the node?
        EXPECT_EQ(pg.get_node(prg_id_1)->reads.begin()->get()->id,
            read_id_1); // is the read really inserted in the node?
        EXPECT_EQ(pg.get_node(prg_id_2)->reads.begin()->get()->id,
            read_id_1); // is the read really inserted in the node?
        EXPECT_EQ(pg.reads.size(), 1);
        EXPECT_EQ(pg.get_read(read_id_1)->id, read_id_1); // is the read really
                                                          // inserted?
        EXPECT_EQ(pg.get_read(read_id_1)->get_hits_as_unordered_map().size(),
            2); // we should have another hit here
        EXPECT_EQ(*pg.get_read(read_id_1)->get_hits_as_unordered_map()[prg_id_1][0],
            *minimizer_hit_1); // is the hit really inserted?
        EXPECT_EQ(*pg.get_read(read_id_1)->get_hits_as_unordered_map()[prg_id_2][0],
            *minimizer_hit_2); // is the hit really inserted?
        EXPECT_EQ(pg.get_read(read_id_1)->node_orientations.size(),
            2); // is the node_orientation was inserted in the read?
        EXPECT_EQ(pg.get_read(read_id_1)->node_orientations[0],
            true); // is the node_orientation was inserted in the read?
        EXPECT_EQ(pg.get_read(read_id_1)->node_orientations[1],
            true); // is the node_orientation was inserted in the read?
    }

    // test 3: add a cluster with the same prg_id, but different read_id
    uint32_t read_id_3 = 101;
    MiniRecord mr3;
    MinimizerHitPtr minimizer_hit_3;
    std::shared_ptr<std::set<MinimizerHitPtr, pComp>> cluster_pointer_3;
    std::shared_ptr<LocalPRG> prg_pointer_3;
    {
        setup_minimizerhit_cluster_prg_function(prg_id_1, read_id_3, &mr3,
            &minimizer_hit_3, &cluster_pointer_3, &prg_pointer_3);
        pg.add_hits_between_PRG_and_read(prg_pointer_3, read_id_3, *cluster_pointer_3);

        // test if everything is fine
        EXPECT_EQ(pg.nodes.size(), 2); // nothing should change from the previous test
        EXPECT_EQ(pg.get_node(prg_id_1)->node_id,
            prg_id_1); // nothing should change from the previous test
        EXPECT_EQ(pg.get_node(prg_id_2)->node_id,
            prg_id_2); // nothing should change from the previous test
        EXPECT_EQ(pg.get_node(prg_id_1)->covg, 2); // coverage should change to 2
        EXPECT_EQ(pg.get_node(prg_id_2)->covg,
            1); // nothing should change from the previous test
        EXPECT_EQ(pg.get_node(prg_id_1)->reads.size(),
            2); // is the read really inserted in the node?
        EXPECT_EQ(pg.get_node(prg_id_2)->reads.size(),
            1); // is the read really inserted in the node?

        std::vector<ReadPtr> reads_in_node_as_vector(
            pg.get_node(prg_id_1)->reads.begin(), pg.get_node(prg_id_1)->reads.end());
        EXPECT_TRUE((reads_in_node_as_vector[0]->id == read_id_1
                        && reads_in_node_as_vector[1]->id == read_id_3)
            || (reads_in_node_as_vector[0]->id == read_id_3
                && reads_in_node_as_vector[1]->id == read_id_1));

        EXPECT_EQ(pg.reads.size(), 2);
        EXPECT_EQ(pg.get_read(read_id_3)->id, read_id_3);
        EXPECT_EQ(pg.get_read(read_id_3)->get_hits_as_unordered_map().size(), 1);
        EXPECT_EQ(*pg.get_read(read_id_3)->get_hits_as_unordered_map()[prg_id_1][0],
            *minimizer_hit_3);
        EXPECT_EQ(pg.get_read(read_id_3)->get_nodes().size(), 1);
        EXPECT_EQ(pg.get_read(read_id_3)->get_nodes()[0].lock()->node_id, prg_id_1);
        EXPECT_EQ(pg.get_read(read_id_3)->node_orientations.size(), 1);
        EXPECT_EQ(pg.get_read(read_id_3)->node_orientations[0], true);
    }

    // TODO: improve this? - add the cluster again, add another cluster, other
    // situations, etc...
}

TEST(PangenomeGraphAddNode, NodeDoesntAlreadyExist_PangenomeGraphNodesContainsNodeId)
{
    PGraphTester pg;

    EXPECT_EQ(pg.nodes.size(), (uint)0);

    uint32_t node_id = 5;
    uint32_t prg_id = 10;
    auto prg_pointer = std::make_shared<LocalPRG>(node_id, "prg_name", "");
    pg.add_node(prg_pointer, node_id);

    NodePtr pan_node = std::make_shared<pangenome::Node>(prg_pointer);
    EXPECT_EQ(*pg.nodes[node_id], *pan_node);
    EXPECT_EQ(pg.nodes[node_id]->node_id, (uint)node_id);
    EXPECT_EQ(pg.nodes[node_id]->prg_id, (uint)node_id);
    EXPECT_EQ(pg.nodes[node_id]->name, "prg_name");
    EXPECT_EQ(pg.nodes[node_id]->covg, (uint)0);
    EXPECT_EQ(pg.nodes[node_id]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[node_id]->samples.size(), (uint)0);

    auto result = pg.nodes.find(node_id) != pg.nodes.end();
    EXPECT_TRUE(result);

    result = pg.nodes.find(prg_id) == pg.nodes.end();
    EXPECT_TRUE(result);
}

/* this test is now comprised on TEST(PangenomeGraph_add_hits_between_PRG_and_read,
AddClusters) TEST(PangenomeGraphAddCoverage,
NodeDoesntAlreadyExist_PangenomeGraphNodeContainsReadPtr) { PGraphTester pg;

    uint32_t read_id = 2;
    pg.add_read(read_id);
    ReadPtr read_ptr = pg.get_read(read_id);

    std::set<MinimizerHitPtr, pComp> mhs;
    uint32_t node_id = 0;
    uint32_t prg_id = 1;

    NodePtr node_ptr = pg.add_coverage(read_ptr, node_id, prg_id, "0");

    auto result = node_ptr->reads.find(read_ptr) != node_ptr->reads.end();
    EXPECT_TRUE(result);
}
 */

TEST(PangenomeGraph_add_hits_between_PRG_and_read, AddTheSameClusterTwice)
{
    PGraphTester pg;

    // test 1: add a cluster
    uint32_t prg_id_1 = 1;
    uint32_t read_id_1 = 100;
    MiniRecord mr1;
    MinimizerHitPtr minimizer_hit_1;
    std::shared_ptr<std::set<MinimizerHitPtr, pComp>> cluster_pointer_1;
    std::shared_ptr<LocalPRG> prg_pointer_1;
    {
        setup_minimizerhit_cluster_prg_function(prg_id_1, read_id_1, &mr1,
            &minimizer_hit_1, &cluster_pointer_1, &prg_pointer_1);
        pg.add_hits_between_PRG_and_read(prg_pointer_1, read_id_1, *cluster_pointer_1);

        // test if everything is fine
        EXPECT_EQ(pg.nodes.size(), 1); // is the node really inserted?
        EXPECT_EQ(
            pg.get_node(prg_id_1)->node_id, prg_id_1); // is the node really inserted?
        EXPECT_EQ(pg.get_node(prg_id_1)->covg, 1); // is the coverage of the node 1?
        EXPECT_EQ(pg.get_node(prg_id_1)->reads.size(),
            1); // is the read really inserted in the node?
        EXPECT_EQ(pg.get_node(prg_id_1)->reads.begin()->get()->id,
            read_id_1); // is the read really inserted in the node?
        EXPECT_EQ(pg.reads.size(), 1);
        EXPECT_EQ(pg.get_read(read_id_1)->id, read_id_1); // is the read really
                                                          // inserted?
        EXPECT_EQ(pg.get_read(read_id_1)->get_hits_as_unordered_map().size(),
            1); // is the hit really inserted?
        EXPECT_EQ(*pg.get_read(read_id_1)->get_hits_as_unordered_map()[prg_id_1][0],
            *minimizer_hit_1); // is the hit really inserted?
        EXPECT_EQ(pg.get_read(read_id_1)
                      ->get_hits_as_unordered_map()[prg_id_1][0]
                      ->get_kmer_node_id(),
            prg_id_1); // is the node really inserted in the read?
        EXPECT_EQ(pg.get_read(read_id_1)->node_orientations.size(),
            1); // is the node_orientation was inserted in the read?
        EXPECT_EQ(pg.get_read(read_id_1)->node_orientations[0],
            true); // is the node_orientation was inserted in the read?

        // add the cluster again
        EXPECT_DEATH(pg.add_hits_between_PRG_and_read(
                         prg_pointer_1, read_id_1, *cluster_pointer_1),
            "");

        /*
        EXPECT_EQ(pg.nodes.size(), 1); //should not change
        EXPECT_EQ(pg.get_node(prg_id_1)->node_id, prg_id_1); //should not change
        EXPECT_EQ(pg.get_node(prg_id_1)->covg, 2); //coverage should increase (TODO:
        reflect if this is really what should happen, since the read id is exactly the
        same) EXPECT_EQ(pg.get_node(prg_id_1)->reads.size(), 2); //we should now have 2
        reads EXPECT_EQ(pg.get_node(prg_id_1)->reads.begin()->get()->id, read_id_1);
        //should not change
        EXPECT_EQ((++pg.get_node(prg_id_1)->reads.begin())->get()->id, read_id_1); //the
        new read id should be also read_id_1 EXPECT_EQ(pg.reads.size(), 1); //should not
        change EXPECT_EQ(pg.get_read(read_id_1)->id, read_id_1); //should not change
        EXPECT_EQ(pg.get_read(read_id_1)->get_hits_as_unordered_map().size(), 1);
        //should not change
        EXPECT_EQ(*pg.get_read(read_id_1)->get_hits_as_unordered_map()[prg_id_1][0],
        *minimizer_hit_1); //should not change
        EXPECT_EQ(pg.get_read(read_id_1)->get_hits_as_unordered_map()[prg_id_1][0]->get_kmer_node_id(),
        prg_id_1); //should not change
        EXPECT_EQ(pg.get_read(read_id_1)->node_orientations.size(), 1); //should not
        change EXPECT_EQ(pg.get_read(read_id_1)->node_orientations[0], true); //should
        not change
        */
    }
}

/* this test is now comprised on TEST(PangenomeGraph_add_hits_between_PRG_and_read,
AddTheSameClusterTwice) TEST(PangenomeGraphAddCoverage,
NodeAlreadyExists_PangenomeGraphNodeReadsContainsReadTwice) { PGraphTester pg;

    uint32_t read_id = 2;
    pg.add_read(read_id);
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
*/

TEST(PangenomeGraphAddNode, AddClusterWrongReadId_AssertCatches)
{
    uint32_t read_id = 1;
    uint32_t not_read_id = 7;

    std::deque<Interval> raw_path = { Interval(7, 8), Interval(10, 14) };
    prg::Path path;
    path.initialize(raw_path);

    Interval interval(0, 5);
    Minimizer m1(0, interval.start, interval.get_end(), 0); // kmer, start, end, strand

    uint32_t prg_id = 4;
    MiniRecord mr1(prg_id, path, 0, 0);
    MinimizerHitPtr minimizer_hit(std::make_shared<MinimizerHit>(not_read_id, m1, mr1));

    std::set<MinimizerHitPtr, pComp> cluster;
    cluster.insert(minimizer_hit);

    PGraphTester pg;
    auto prg_pointer = std::make_shared<LocalPRG>(prg_id, "", "");

    EXPECT_DEATH(pg.add_hits_between_PRG_and_read(prg_pointer, read_id, cluster), "");
}

TEST(PangenomeGraphAddNode, AddClusterWrongPrgId_AssertCatches)
{
    uint32_t read_id = 1;
    Read read(read_id);

    std::set<MinimizerHitPtr, pComp> cluster;
    uint32_t prg_id = 4;
    uint32_t not_prg_id = 7;
    Interval interval(0, 5);
    std::deque<Interval> raw_path = { Interval(7, 8), Interval(10, 14) };
    prg::Path path;
    path.initialize(raw_path);
    Minimizer m1(0, interval.start, interval.get_end(), 0); // kmer, start, end, strand
    MiniRecord mr1(not_prg_id, path, 0, 0);
    MinimizerHitPtr minimizer_hit(std::make_shared<MinimizerHit>(read_id, m1, mr1));
    cluster.insert(minimizer_hit);

    PGraphTester pg;
    auto prg_pointer = std::make_shared<LocalPRG>(prg_id, "", "");

    EXPECT_DEATH(pg.add_hits_between_PRG_and_read(prg_pointer, read_id, cluster), "");
}

/* this test is now comprised in TEST(PangenomeGraphNode, add_node_and_get_node)
TEST(PangenomeGraphAddNode, AddNode_PangenomeGraphNodesContainsNodeId) {
    std::set<MinimizerHitPtr, pComp> mhs;
    PGraphTester pg;
    uint32_t node_id = 0;
    uint32_t read_id = 1;
    pg.add_node(node_id, "0", read_id, mhs);

    auto result = pg.nodes.find(node_id) != pg.nodes.end();
    EXPECT_TRUE(result);
}
*/

TEST(PangenomeGraphAddNode, AddNode_PangenomeGraphNodeHasRightProperties)
{
    PGraphTester pg;
    uint32_t node_id = 15;
    auto prg_pointer = std::make_shared<LocalPRG>(node_id, "prg_name", "");
    pg.add_node(prg_pointer);

    NodePtr pan_node = std::make_shared<pangenome::Node>(prg_pointer);
    EXPECT_EQ(*pg.nodes[node_id], *pan_node);
    EXPECT_EQ(pg.nodes[node_id]->node_id, (uint)node_id);
    EXPECT_EQ(pg.nodes[node_id]->prg_id, (uint)node_id);
    EXPECT_EQ(pg.nodes[node_id]->name, "prg_name");
    EXPECT_EQ(pg.nodes[node_id]->covg, (uint)0);
    EXPECT_EQ(pg.nodes[node_id]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[node_id]->samples.size(), (uint)0);

    // TODO: fix this
    // EXPECT_EQ(pg.nodes[node_id]->kmer_prg_with_coverage,
    // KmerGraphWithCoverage(&prg_pointer->kmer_prg));
}

TEST(PangenomeGraphAddNode, AddNode_PangenomeGraphReadHasRightProperties)
{
    std::set<MinimizerHitPtr, pComp> dummy_cluster;
    PGraphTester pg;
    uint32_t node_id = 0;
    uint32_t read_id = 1;
    auto prg_pointer = std::make_shared<LocalPRG>(node_id, "zero", "");
    pg.add_hits_between_PRG_and_read(prg_pointer, read_id, dummy_cluster);

    ReadPtr pr = std::make_shared<Read>(read_id);
    EXPECT_EQ(pg.reads.size(), (uint)1);
    EXPECT_EQ(*pg.reads[1], *pr);
    EXPECT_EQ(pg.reads[1]->get_hits_as_unordered_map().size(), (uint)1);
    EXPECT_EQ(pg.reads[1]->get_hits_as_unordered_map()[0].size(), (uint)0);
}

TEST(PangenomeGraphTest, add_node_sample)
{
    // add node and check it's there
    PGraphTester pg({ "sample", "sample1" });

    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA"));
    std::vector<KmerNodePtr> kmp;

    pg.add_node(l0);
    pg.add_hits_between_PRG_and_sample(0, "sample", kmp);

    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)1);

    EXPECT_EQ(pg.samples.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)1);

    EXPECT_EQ(pg.reads.size(), (uint)0);

    // add a second time
    pg.add_hits_between_PRG_and_sample(0, "sample", kmp);
    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)1);

    EXPECT_EQ(pg.samples.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)2);

    EXPECT_EQ(pg.reads.size(), (uint)0);

    // add a node with a different sample
    pg.add_hits_between_PRG_and_sample(0, "sample1", kmp);
    EXPECT_EQ(pg.nodes.size(), (uint)1);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)2);

    EXPECT_EQ(pg.samples.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)2);
    EXPECT_EQ(pg.samples["sample1"]->name, "sample1");
    EXPECT_EQ(pg.samples["sample1"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample1"]->paths[0].size(), (uint)1);

    EXPECT_EQ(pg.reads.size(), (uint)0);

    // add a node with a different prg
    auto l1 = std::make_shared<LocalPRG>(LocalPRG(1, "one", "AGCTGCTAGCTTCGGACGCACA"));
    pg.add_node(l1);
    pg.add_hits_between_PRG_and_sample(1, "sample1", kmp);
    EXPECT_EQ(pg.nodes.size(), (uint)2);
    EXPECT_EQ(pg.nodes[0]->node_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg.nodes[0]->name, "zero");
    EXPECT_EQ(pg.nodes[0]->covg, (uint)3);
    EXPECT_EQ(pg.nodes[0]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[0]->samples.size(), (uint)2);
    EXPECT_EQ(pg.nodes[1]->node_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg.nodes[1]->name, "one");
    EXPECT_EQ(pg.nodes[1]->covg, (uint)1);
    EXPECT_EQ(pg.nodes[1]->reads.size(), (uint)0);
    EXPECT_EQ(pg.nodes[1]->samples.size(), (uint)1);

    EXPECT_EQ(pg.samples.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample"]->name, "sample");
    EXPECT_EQ(pg.samples["sample"]->paths.size(), (uint)1);
    EXPECT_EQ(pg.samples["sample"]->paths[0].size(), (uint)2);
    EXPECT_EQ(pg.samples["sample1"]->name, "sample1");
    EXPECT_EQ(pg.samples["sample1"]->paths.size(), (uint)2);
    EXPECT_EQ(pg.samples["sample1"]->paths[0].size(), (uint)1);
    EXPECT_EQ(pg.samples["sample1"]->paths[1].size(), (uint)1);

    EXPECT_EQ(pg.reads.size(), (uint)0);
}

/* this test is not needed anymore as we don't have a clear function anymore
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
 */

TEST(PangenomeGraphTest, equals)
{
    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "0", ""));
    auto l1 = std::make_shared<LocalPRG>(LocalPRG(1, "1", ""));
    auto l2 = std::make_shared<LocalPRG>(LocalPRG(2, "2", ""));
    auto l3 = std::make_shared<LocalPRG>(LocalPRG(3, "3", ""));

    PGraphTester pg1;
    pg1.add_node(l0);
    pg1.add_node(l1);
    pg1.add_node(l1);
    pg1.add_node(l2);

    PGraphTester pg2;
    pg2.add_node(l1);
    pg2.add_node(l0);
    pg2.add_node(l2);
    pg2.add_node(l1);

    // adding nodes in different order should make no difference
    EXPECT_EQ(pg1, pg1);
    EXPECT_EQ(pg2, pg2);
    EXPECT_EQ(pg1, pg2);
    EXPECT_EQ(pg2, pg1);

    // should not matter if node_id is different provided prg_id is same
    pg2.nodes[7] = std::make_shared<pangenome::Node>(l2, 7);
    pg2.nodes.erase(2);
    EXPECT_EQ(pg2, pg2);
    EXPECT_EQ(pg1, pg2);
    EXPECT_EQ(pg2, pg1);

    // or one extra node
    pg2.add_node(l3);
    EXPECT_EQ((pg1 == pg2), false);
    EXPECT_EQ((pg2 == pg1), false);

    // should not break when have a cycle in pangraph
    pg1.add_node(l0);
    EXPECT_EQ(pg1, pg1);
}

TEST(PangenomeGraphTest, not_equals)
{
    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "0", ""));
    auto l1 = std::make_shared<LocalPRG>(LocalPRG(1, "1", ""));
    auto l2 = std::make_shared<LocalPRG>(LocalPRG(2, "2", ""));
    auto l3 = std::make_shared<LocalPRG>(LocalPRG(3, "3", ""));

    PGraphTester pg1;
    pg1.add_node(l0);
    pg1.add_node(l1);
    pg1.add_node(l1);
    pg1.add_node(l2);

    PGraphTester pg2;
    pg2.add_node(l1);
    pg2.add_node(l0);
    pg2.add_node(l2);
    pg2.add_node(l1);

    // adding nodes in different order should make no difference
    EXPECT_EQ((pg1 != pg1), false);
    EXPECT_EQ((pg2 != pg2), false);
    EXPECT_EQ((pg1 != pg2), false);
    EXPECT_EQ((pg2 != pg1), false);

    // or one extra node
    pg2.add_node(l3);
    EXPECT_EQ((pg1 != pg2), true);
    EXPECT_EQ((pg2 != pg1), true);

    // should not break when have a cycle in pangraph
    pg1.add_node(l0);
    EXPECT_EQ((pg1 != pg1), false);
}

TEST(PangenomeGraphTest, remove_node)
{
    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "0", ""));
    auto l1 = std::make_shared<LocalPRG>(LocalPRG(1, "1", ""));
    auto l2 = std::make_shared<LocalPRG>(LocalPRG(2, "2", ""));
    auto l3 = std::make_shared<LocalPRG>(LocalPRG(3, "3", ""));
    std::set<MinimizerHitPtr, pComp> dummy_cluster;

    PGraphTester pg1, pg2;
    // read 0: 0->1->2->3
    pg1.add_node(l0);
    pg1.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg1.add_node(l1);
    pg1.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg1.add_node(l2);
    pg1.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg1.add_node(l3);
    pg1.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    // read 0: 0->1->3
    pg2.add_node(l0);
    pg1.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg2.add_node(l1);
    pg1.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg2.add_node(l3);
    pg1.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    pg1.remove_node(pg1.nodes[2]);
    EXPECT_EQ(pg1, pg2);
}

TEST(PangenomeGraphTest, remove_read)
{
    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "0", ""));
    auto l1 = std::make_shared<LocalPRG>(LocalPRG(1, "1", ""));
    auto l2 = std::make_shared<LocalPRG>(LocalPRG(2, "2", ""));
    auto l3 = std::make_shared<LocalPRG>(LocalPRG(3, "3", ""));
    auto l4 = std::make_shared<LocalPRG>(LocalPRG(4, "4", ""));
    auto l5 = std::make_shared<LocalPRG>(LocalPRG(5, "5", ""));
    std::set<MinimizerHitPtr, pComp> dummy_cluster;

    PGraphTester pg1, pg2, pg3;
    // read 0: 0->1->2->3
    pg1.add_node(l0);
    pg1.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg1.add_node(l1);
    pg1.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg1.add_node(l2);
    pg1.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg1.add_node(l3);
    pg1.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    // read 1: 4->5->0->5
    pg1.add_node(l4);
    pg1.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg1.add_node(l5);
    pg1.add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pg1.add_node(l0);
    pg1.add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg1.add_node(l5);
    pg1.add_hits_between_PRG_and_read(l5, 1, dummy_cluster);

    // read 1: 4->5->0->5
    pg2.add_node(l4);
    pg2.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg2.add_node(l5);
    pg2.add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pg2.add_node(l0);
    pg2.add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg2.add_node(l5);
    pg2.add_hits_between_PRG_and_read(l5, 1, dummy_cluster);

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

TEST(PangenomeGraphTest, remove_low_covg_nodes)
{
    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "0", ""));
    auto l1 = std::make_shared<LocalPRG>(LocalPRG(1, "1", ""));
    auto l2 = std::make_shared<LocalPRG>(LocalPRG(2, "2", ""));
    auto l3 = std::make_shared<LocalPRG>(LocalPRG(3, "3", ""));
    auto l4 = std::make_shared<LocalPRG>(LocalPRG(4, "4", ""));
    auto l5 = std::make_shared<LocalPRG>(LocalPRG(5, "5", ""));
    std::set<MinimizerHitPtr, pComp> dummy_cluster;

    PGraphTester pg1, pg2, pg3;
    // read 0: 0->1->2->3
    pg1.add_node(l0);
    pg1.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg1.add_node(l1);
    pg1.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg1.add_node(l2);
    pg1.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg1.add_node(l3);
    pg1.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    // read 1: -4 -> -3 -> -1
    pg1.add_node(l1);
    pg1.add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pg1.add_node(l3);
    pg1.add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pg1.add_node(l4);
    pg1.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);

    // read 2: 0 -> 1 -> 3 -> 4
    pg1.add_node(l0);
    pg1.add_hits_between_PRG_and_read(l0, 2, dummy_cluster);
    pg1.add_node(l1);
    pg1.add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pg1.add_node(l3);
    pg1.add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pg1.add_node(l4);
    pg1.add_hits_between_PRG_and_read(l4, 2, dummy_cluster);

    // read 3: 0 -> 5
    pg1.add_node(l0);
    pg1.add_hits_between_PRG_and_read(l0, 3, dummy_cluster);
    pg1.add_node(l5);
    pg1.add_hits_between_PRG_and_read(l5, 3, dummy_cluster);

    // read 4: 5 -> 1
    pg1.add_node(l5);
    pg1.add_hits_between_PRG_and_read(l5, 4, dummy_cluster);
    pg1.add_node(l1);
    pg1.add_hits_between_PRG_and_read(l1, 4, dummy_cluster);

    // read 0: 0->1->3
    pg2.add_node(l0);
    pg2.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg2.add_node(l1);
    pg2.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg2.add_node(l3);
    pg2.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    // read 1: -4 -> -3 -> -1
    pg2.add_node(l1);
    pg2.add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pg2.add_node(l3);
    pg2.add_hits_between_PRG_and_read(l3, 1, dummy_cluster);
    pg2.add_node(l4);
    pg2.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);

    // read 2: 0 -> 1 -> 3 -> 4
    pg2.add_node(l0);
    pg2.add_hits_between_PRG_and_read(l0, 2, dummy_cluster);
    pg2.add_node(l1);
    pg2.add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pg2.add_node(l3);
    pg2.add_hits_between_PRG_and_read(l3, 2, dummy_cluster);
    pg2.add_node(l4);
    pg2.add_hits_between_PRG_and_read(l4, 2, dummy_cluster);

    // read 3: 0 -> 5
    pg2.add_node(l0);
    pg2.add_hits_between_PRG_and_read(l0, 3, dummy_cluster);
    pg2.add_node(l5);
    pg2.add_hits_between_PRG_and_read(l5, 3, dummy_cluster);

    // read 4: 5 -> 1
    pg2.add_node(l5);
    pg2.add_hits_between_PRG_and_read(l5, 4, dummy_cluster);
    pg2.add_node(l1);
    pg2.add_hits_between_PRG_and_read(l1, 4, dummy_cluster);

    pg1.remove_low_covg_nodes(1);
    EXPECT_EQ(pg1, pg2);

    // read 0: 0->1->3
    pg3.add_node(l0);
    pg3.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg3.add_node(l1);
    pg3.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg3.add_node(l3);
    pg3.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    // read 1: 1 -> 3
    pg3.add_node(l1);
    pg3.add_hits_between_PRG_and_read(l1, 1, dummy_cluster);
    pg3.add_node(l3);
    pg3.add_hits_between_PRG_and_read(l3, 1, dummy_cluster);

    // read 2: 0 -> 1 -> 3
    pg3.add_node(l0);
    pg3.add_hits_between_PRG_and_read(l0, 2, dummy_cluster);
    pg3.add_node(l1);
    pg3.add_hits_between_PRG_and_read(l1, 2, dummy_cluster);
    pg3.add_node(l3);
    pg3.add_hits_between_PRG_and_read(l3, 2, dummy_cluster);

    // read 3: 0
    pg3.add_node(l0);
    pg3.add_hits_between_PRG_and_read(l0, 3, dummy_cluster);

    // read 4: 1
    pg3.add_node(l1);
    pg3.add_hits_between_PRG_and_read(l1, 4, dummy_cluster);

    pg1.remove_low_covg_nodes(2);
    EXPECT_EQ(pg1, pg3);
}

TEST(PangenomeGraphTest, split_node_by_reads)
{
    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "0", ""));
    auto l1 = std::make_shared<LocalPRG>(LocalPRG(1, "1", ""));
    auto l2 = std::make_shared<LocalPRG>(LocalPRG(2, "2", ""));
    auto l3 = std::make_shared<LocalPRG>(LocalPRG(3, "3", ""));
    auto l4 = std::make_shared<LocalPRG>(LocalPRG(4, "4", ""));
    auto l5 = std::make_shared<LocalPRG>(LocalPRG(5, "5", ""));
    std::set<MinimizerHitPtr, pComp> dummy_cluster;

    PGraphTester pg1, pg2, pg3;
    // read 0: 0->1->2->3
    pg1.add_node(l0);
    pg1.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg1.add_node(l1);
    pg1.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    pg1.add_node(l2);
    pg1.add_hits_between_PRG_and_read(l2, 0, dummy_cluster);
    pg1.add_node(l3);
    pg1.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    // read 1: 4->5->0->5
    pg1.add_node(l4);
    pg1.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg1.add_node(l5);
    pg1.add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pg1.add_node(l0);
    pg1.add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg1.add_node(l5);
    pg1.add_hits_between_PRG_and_read(l5, 1, dummy_cluster);

    EXPECT_EQ((uint)6, pg1.nodes.size());
    EXPECT_EQ(pg1.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg1.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg1.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg1.nodes[1]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[2]->prg_id, (uint)2);
    EXPECT_EQ(pg1.nodes[2]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg1.nodes[3]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg1.nodes[4]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[5]->prg_id, (uint)5);
    EXPECT_EQ(pg1.nodes[5]->covg, (uint)2);

    // read 0: 0->1->2->3
    pg2.add_node(l0);
    pg2.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg2.add_node(l1);
    pg2.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    NodePtr n = std::make_shared<pangenome::Node>(l2, 7);
    pg2.nodes[7] = n;
    pg2.add_node(l3);
    pg2.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    // read 1: 4->5->0->5
    pg2.add_node(l4);
    pg2.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    pg2.add_node(l5);
    pg2.add_hits_between_PRG_and_read(l5, 1, dummy_cluster);
    pg2.add_node(l0);
    pg2.add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg2.add_node(l5);
    pg2.add_hits_between_PRG_and_read(l5, 1, dummy_cluster);

    std::unordered_set<ReadPtr> reads = { pg1.reads[0] };
    std::vector<uint_least32_t> node_ids = { 1, 2, 3 };
    std::vector<uint_least32_t> node_ids_exp = { 1, 6, 3 };
    std::vector<bool> node_orients = { 0, 0, 0 };
    pg1.split_node_by_reads(reads, node_ids, node_orients, 2);
    EXPECT_EQ(pg1, pg2);
    EXPECT_ITERABLE_EQ(std::vector<uint_least32_t>, node_ids_exp, node_ids);

    EXPECT_EQ((uint)6, pg1.nodes.size());
    EXPECT_EQ(pg1.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg1.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg1.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg1.nodes[1]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[6]->prg_id, (uint)2);
    EXPECT_EQ(pg1.nodes[6]->covg, (uint)0);
    EXPECT_EQ(pg1.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg1.nodes[3]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg1.nodes[4]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[5]->prg_id, (uint)5);
    EXPECT_EQ(pg1.nodes[5]->covg, (uint)2);

    // read 0: 0->1->2->3
    pg3.add_node(l0);
    pg3.add_hits_between_PRG_and_read(l0, 0, dummy_cluster);
    pg3.add_node(l1);
    pg3.add_hits_between_PRG_and_read(l1, 0, dummy_cluster);
    n = std::make_shared<pangenome::Node>(l2, 7);
    pg3.nodes[7] = n;
    pg3.add_node(l3);
    pg3.add_hits_between_PRG_and_read(l3, 0, dummy_cluster);

    // read 1: 4->5->0->5
    pg3.add_node(l4);
    pg3.add_hits_between_PRG_and_read(l4, 1, dummy_cluster);
    n = std::make_shared<pangenome::Node>(l5, 8);
    pg3.nodes[8] = n;
    pg3.add_node(l0);
    pg3.add_hits_between_PRG_and_read(l0, 1, dummy_cluster);
    pg3.add_node(l5);
    pg3.add_hits_between_PRG_and_read(l5, 1, dummy_cluster);

    reads = { pg1.reads[1] };
    node_ids = { 5, 0, 5 };
    node_ids_exp = { 7, 0, 5 };
    pg1.split_node_by_reads(reads, node_ids, node_orients, 5);
    EXPECT_EQ(pg1, pg3);
    EXPECT_ITERABLE_EQ(std::vector<uint_least32_t>, node_ids_exp, node_ids);

    EXPECT_EQ((uint)7, pg1.nodes.size());
    EXPECT_EQ(pg1.nodes[0]->prg_id, (uint)0);
    EXPECT_EQ(pg1.nodes[0]->covg, (uint)2);
    EXPECT_EQ(pg1.nodes[1]->prg_id, (uint)1);
    EXPECT_EQ(pg1.nodes[1]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[6]->prg_id, (uint)2);
    EXPECT_EQ(pg1.nodes[6]->covg, (uint)0);
    EXPECT_EQ(pg1.nodes[3]->prg_id, (uint)3);
    EXPECT_EQ(pg1.nodes[3]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[4]->prg_id, (uint)4);
    EXPECT_EQ(pg1.nodes[4]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[5]->prg_id, (uint)5);
    EXPECT_EQ(pg1.nodes[5]->covg, (uint)1);
    EXPECT_EQ(pg1.nodes[7]->prg_id, (uint)5);
    EXPECT_EQ(pg1.nodes[7]->covg, (uint)0);
}

TEST(PangenomeGraphTest, add_hits_to_kmergraph) { }

TEST(PangenomeGraphTest, save_matrix)
{
    // add node and check it's there
    std::vector<std::string> names = { "sample1", "sample2", "sample3", "sample4" };
    PGraphTester pg(names);

    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "zero", "AGCTGCTAGCTTCGGACGCACA"));
    auto l1 = std::make_shared<LocalPRG>(LocalPRG(1, "one", ""));
    auto l2 = std::make_shared<LocalPRG>(LocalPRG(2, "two", ""));
    std::vector<KmerNodePtr> kmp;

    pg.add_node(l0);
    pg.add_node(l1);
    pg.add_node(l2);

    pg.add_hits_between_PRG_and_sample(0, "sample1", kmp);
    pg.add_hits_between_PRG_and_sample(0, "sample1", kmp);
    pg.add_hits_between_PRG_and_sample(0, "sample2", kmp);
    pg.add_hits_between_PRG_and_sample(1, "sample1", kmp);
    pg.add_hits_between_PRG_and_sample(2, "sample3", kmp);

    pg.save_matrix(TEST_CASE_DIR + "pangraph_test_save.matrix", names);
}

TEST(PangenomeGraphTest, save_mapped_read_strings)
{
    PGraphTester pg;
    pangenome::ReadPtr pr;
    MinimizerHits mhits;
    std::deque<Interval> d;
    prg::Path p;

    // read1
    Minimizer m1(0, 1, 6, 0); // kmer, start, end, strand
    d = { Interval(7, 8), Interval(10, 14) };
    p.initialize(d);
    MiniRecord mr1(0, p, 0, 0);
    mhits.add_hit(1, m1, mr1); // read 1

    Minimizer m2(0, 0, 5, 0);
    d = { Interval(6, 10), Interval(11, 12) };
    p.initialize(d);
    MiniRecord mr2(0, p, 0, 0);
    mhits.add_hit(1, m2, mr2);

    Minimizer m3(0, 0, 5, 0);
    d = { Interval(6, 10), Interval(12, 13) };
    p.initialize(d);
    MiniRecord mr3(0, p, 0, 0);
    mhits.add_hit(1, m3, mr3);

    auto l0 = std::make_shared<LocalPRG>(LocalPRG(0, "zero", ""));
    pg.add_node(l0);
    pg.add_hits_between_PRG_and_read(l0, 1, mhits.hits);
    mhits.clear();

    // read 2
    Minimizer m4(0, 2, 7, 1);
    d = { Interval(6, 10), Interval(11, 12) };
    p.initialize(d);
    MiniRecord mr4(0, p, 0, 0);
    mhits.add_hit(2, m4, mr4);

    Minimizer m5(0, 5, 10, 1);
    d = { Interval(6, 10), Interval(12, 13) };
    p.initialize(d);
    MiniRecord mr5(0, p, 0, 0);
    mhits.add_hit(2, m5, mr5);

    pg.add_hits_between_PRG_and_read(l0, 2, mhits.hits);

    std::string expected1
        = ">read1 pandora: 1 0:6 + \nshould\n>read2 pandora: 2 2:10 - \nis time \n";
    std::string expected2
        = ">read2 pandora: 2 2:10 - \nis time \n>read1 pandora: 1 0:6 + \nshould\n";

    pg.save_mapped_read_strings(TEST_CASE_DIR + "reads.fa", "save_mapped_read_strings");
    std::ifstream ifs("save_mapped_read_strings/zero/zero.reads.fa");
    std::string content(
        (std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    EXPECT_TRUE((content == expected1) or (content == expected2));

    pg.save_mapped_read_strings(TEST_CASE_DIR + "reads.fa", ".");
    std::ifstream ifs2("zero/zero.reads.fa");
    std::string content2(
        (std::istreambuf_iterator<char>(ifs2)), (std::istreambuf_iterator<char>()));
    EXPECT_TRUE((content2 == expected1) or (content2 == expected2));
}

TEST(PangenomeGraphTest, get_node_closest_vcf_reference_no_paths)
{
    uint32_t prg_id = 3, w = 1, k = 3, max_num_kmers_to_average = 100;
    std::string prg_name = "nested varsite";
    auto l3 = std::make_shared<LocalPRG>(prg_id, prg_name, "A 5 G 7 C 8 T 7  6 G 5 T");
    auto index = std::make_shared<Index>();
    l3->minimizer_sketch(index, w, k);

    std::string sample_name = "null_test_sample";
    pangenome::Graph pangraph({ sample_name });
    std::vector<KmerNodePtr> sample_kmer_path = {};

    pangraph.add_node(l3);
    pangraph.add_hits_between_PRG_and_sample(prg_id, sample_name, sample_kmer_path);
    auto path = pangraph.get_node_closest_vcf_reference(
        *pangraph.nodes[prg_id], w, *l3, max_num_kmers_to_average);
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, path, l3->prg.top_path());
}

TEST(PangenomeGraphTest, get_node_closest_vcf_reference_one_path)
{
    uint32_t prg_id = 3, w = 1, k = 3, max_num_kmers_to_average = 100;
    std::string prg_name = "nested varsite";
    auto l3 = std::make_shared<LocalPRG>(prg_id, prg_name, "A 5 G 7 C 8 T 7  6 G 5 T");
    auto index = std::make_shared<Index>();
    l3->minimizer_sketch(index, w, k);

    std::string sample_name = "single_test_sample";
    pangenome::Graph pangraph({ sample_name });

    auto& kg = l3->kmer_prg;
    std::vector<KmerNodePtr> sample_kmer_path
        = { kg.nodes[0], kg.nodes[2], kg.nodes[5], kg.nodes[6] };

    pangraph.add_node(l3);
    pangraph.add_hits_between_PRG_and_sample(prg_id, sample_name, sample_kmer_path);
    auto& node = *pangraph.nodes[prg_id];

    auto path = pangraph.get_node_closest_vcf_reference(
        node, w, *l3, max_num_kmers_to_average);
    std::vector<LocalNodePtr> exp_path = { l3->prg.nodes[0], l3->prg.nodes[1],
        l3->prg.nodes[3], l3->prg.nodes[4], l3->prg.nodes[6] };

    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, path, exp_path);
}

TEST(PangenomeGraphTest, get_node_closest_vcf_reference_three_paths)
{
    uint32_t prg_id = 3, w = 1, k = 3, max_num_kmers_to_average = 100;
    std::string prg_name = "nested varsite";
    auto l3 = std::make_shared<LocalPRG>(prg_id, prg_name, "A 5 G 7 C 8 T 7  6 G 5 T");
    auto index = std::make_shared<Index>();
    l3->minimizer_sketch(index, w, k);

    pangenome::Graph pangraph({ "test_sample1", "test_sample1_again", "test_sample2" });
    auto& kg = l3->kmer_prg;

    std::string sample_name = "test_sample1";
    std::vector<KmerNodePtr> sample_kmer_path
        = { kg.nodes[0], kg.nodes[2], kg.nodes[5], kg.nodes[6] };
    pangraph.add_node(l3);
    pangraph.add_hits_between_PRG_and_sample(prg_id, sample_name, sample_kmer_path);

    sample_name = "test_sample1_again";
    sample_kmer_path = { kg.nodes[0], kg.nodes[2], kg.nodes[5], kg.nodes[6] };
    pangraph.add_hits_between_PRG_and_sample(prg_id, sample_name, sample_kmer_path);

    sample_name = "test_sample2";
    sample_kmer_path = { kg.nodes[0], kg.nodes[4], kg.nodes[6] };
    pangraph.add_hits_between_PRG_and_sample(prg_id, sample_name, sample_kmer_path);

    auto& node = *pangraph.nodes[prg_id];
    auto path = pangraph.get_node_closest_vcf_reference(
        node, w, *l3, max_num_kmers_to_average);
    std::vector<LocalNodePtr> exp_path = { l3->prg.nodes[0], l3->prg.nodes[1],
        l3->prg.nodes[3], l3->prg.nodes[4], l3->prg.nodes[6] };

    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, path, exp_path);
}

TEST(PangenomeGraphTest, copy_coverages_to_kmergraphs)
{
    uint32_t prg_id = 3, w = 1, k = 3;
    std::string prg_name = "nested varsite", sample_name = "sample";
    auto l3 = std::make_shared<LocalPRG>(prg_id, prg_name, "A 5 G 7 C 8 T 7  6 G 5 T");
    auto index = std::make_shared<Index>();
    l3->minimizer_sketch(index, w, k);

    pangenome::Graph ref_pangraph({ sample_name });
    auto sample_id = 0;
    std::vector<KmerNodePtr> empty = {};
    ref_pangraph.add_node(l3);
    ref_pangraph.add_hits_between_PRG_and_sample(prg_id, sample_name, empty);

    EXPECT_TRUE(ref_pangraph.nodes.find(prg_id) != ref_pangraph.nodes.end());
    auto& kg = *(ref_pangraph.nodes[prg_id]->kmer_prg_with_coverage.kmer_prg);
    EXPECT_EQ(kg.nodes.size(), (uint)7);
    auto& kgWithCoverage = ref_pangraph.nodes[prg_id]->kmer_prg_with_coverage;
    kgWithCoverage.set_covg(2, 5, 1, sample_id);
    kgWithCoverage.set_covg(4, 8, 0, sample_id);
    kgWithCoverage.set_covg(5, 2, 1, sample_id);
    kgWithCoverage.set_covg(6, 5, 0, sample_id);

    pangenome::Graph pangraph({ "sample_0", "sample_1", "sample_2", "sample_3" });
    sample_id = 3;
    pangraph.add_node(l3);
    pangraph.add_hits_between_PRG_and_sample(prg_id, "sample_3", empty);

    pangraph.copy_coverages_to_kmergraphs(ref_pangraph, sample_id);

    for (uint32_t id = 0; id < 3; ++id) {
        for (const auto& node :
            pangraph.nodes[prg_id]->kmer_prg_with_coverage.kmer_prg->nodes) {
            EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(
                          node->id, 0, id),
                (uint)0);
            EXPECT_EQ(pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(
                          node->id, 1, id),
                (uint)0);
        }
    }
    auto id = sample_id;
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(0, 0, id), (uint)0);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(0, 1, id), (uint)0);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(1, 0, id), (uint)0);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(1, 1, id), (uint)0);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(2, 0, id), (uint)0);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(2, 1, id), (uint)5);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(3, 0, id), (uint)0);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(3, 1, id), (uint)0);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(4, 0, id), (uint)8);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(4, 1, id), (uint)0);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(5, 0, id), (uint)0);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(5, 1, id), (uint)2);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(6, 0, id), (uint)5);
    EXPECT_EQ(
        pangraph.nodes[prg_id]->kmer_prg_with_coverage.get_covg(6, 1, id), (uint)0);
}

TEST(PangenomeGraphTest, infer_node_vcf_reference_path_no_file_strings)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    std::vector<std::string> prg_strings;
    prg_strings.push_back(
        "ATGCCGGTAATTAAAGTACGTGAAAAGAAACTGGCTC 5 A 6 G 5 CGAAAACGCACGCCGCACTCGTCTGTAC");
    prg_strings.push_back("A 5 G 7 C 8 T 7  6 G 5 T");
    prg_strings.push_back("TC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AG");
    prg_strings.push_back("A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");

    std::string sample_name = "sample";
    pangenome::Graph pangraph({ sample_name });
    std::vector<KmerNodePtr> empty;
    auto index = std::make_shared<Index>();
    uint32_t prg_id = 0, sample_id = 0, w = 1, k = 3, max_num_kmers_to_average = 100;
    std::unordered_map<std::string, std::string> vcf_refs;
    std::vector<std::vector<LocalNodePtr>> vcf_ref_paths;
    for (const auto& prg_string : prg_strings) {
        std::string prg_name = "prg" + std::to_string(prg_id);
        prgs.emplace_back(
            std::make_shared<LocalPRG>(LocalPRG(prg_id, prg_name, prg_string)));
        prgs.back()->minimizer_sketch(index, w, k);
        pangraph.add_node(prgs.back());
        pangraph.add_hits_between_PRG_and_sample(prg_id, sample_name, empty);
        vcf_ref_paths.emplace_back(
            pangraph.infer_node_vcf_reference_path(*pangraph.nodes[prg_id], prgs.back(),
                w, vcf_refs, max_num_kmers_to_average));
        prg_id++;
    }

    EXPECT_EQ(vcf_ref_paths.size(), (uint)4);
    for (uint j = 0; j < vcf_ref_paths.size(); ++j)
        EXPECT_ITERABLE_EQ(
            std::vector<LocalNodePtr>, vcf_ref_paths[j], prgs[j]->prg.top_path());
}

TEST(PangenomeGraphTest, infer_node_vcf_reference_path_with_file_strings)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    std::vector<std::string> prg_strings;
    prg_strings.push_back(
        "ATGCCGGTAATTAAAGTACGTGAAAAGAAACTGGCTC 5 A 6 G 5 CGAAAACGCACGCCGCACTCGTCTGTAC");
    prg_strings.push_back("A 5 G 7 C 8 T 7  6 G 5 T");
    prg_strings.push_back("TC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AG");
    prg_strings.push_back("AATTTTTTTGGGGTTGGTTTTAAA 5 GGGGG 7 CCCCCC 8 TTTTTT 7 TTTTTT "
                          "9 CCGCCGCCGCCG 10 CGGCCGCCG 9  6 GGGGG 5 TATAAAAATTTTTT");
    std::unordered_map<std::string, std::string> vcf_ref_strings;
    vcf_ref_strings["prg0"]
        = "ATGCCGGTAATTAAAGTACGTGAAAAGAAACTGGCTCGCGAAAACGCACGCCGCACTCGTCTGTAC"; // valid
    vcf_ref_strings["prg1"] = "AGT"; // invalid, too short
    vcf_ref_strings["prg2"] = "ATGCCGGTAATTAAAGTACGTGAAAAGAAACTGGCTCGCGAAAACGCACGCCGCAC"
                              "TCGTCTGTAC"; // invalid, is not a path through prg
    vcf_ref_strings["prg3"] = "AATTTTTTTGGGGTTGGTTTTAAAGGGGGTTTTTTTTTTTTCCGCCGCCGCCGTAT"
                              "AAAAATTTTTT"; // valid

    std::string sample_name = "sample";
    pangenome::Graph pangraph({ sample_name });
    std::vector<KmerNodePtr> empty;
    auto index = std::make_shared<Index>();
    uint32_t prg_id = 0, sample_id = 0, w = 1, k = 3, max_num_kmers_to_average = 100;
    std::vector<std::vector<LocalNodePtr>> vcf_ref_paths;
    for (const auto& prg_string : prg_strings) {
        std::string prg_name = "prg" + std::to_string(prg_id);
        prgs.emplace_back(
            std::make_shared<LocalPRG>(LocalPRG(prg_id, prg_name, prg_string)));
        prgs.back()->minimizer_sketch(index, w, k);
        pangraph.add_node(prgs.back());
        pangraph.add_hits_between_PRG_and_sample(prg_id, sample_name, empty);
        vcf_ref_paths.emplace_back(
            pangraph.infer_node_vcf_reference_path(*pangraph.nodes[prg_id], prgs.back(),
                w, vcf_ref_strings, max_num_kmers_to_average));
        prg_id++;
    }

    std::vector<LocalNodePtr> exp_path0
        = { prgs[0]->prg.nodes[0], prgs[0]->prg.nodes[2], prgs[0]->prg.nodes[3] };
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, vcf_ref_paths[0], exp_path0);
    EXPECT_ITERABLE_EQ(
        std::vector<LocalNodePtr>, vcf_ref_paths[1], prgs[1]->prg.top_path());
    EXPECT_ITERABLE_EQ(
        std::vector<LocalNodePtr>, vcf_ref_paths[2], prgs[2]->prg.top_path());
    std::vector<LocalNodePtr> exp_path3 = { prgs[3]->prg.nodes[0],
        prgs[3]->prg.nodes[1], prgs[3]->prg.nodes[3], prgs[3]->prg.nodes[4],
        prgs[3]->prg.nodes[5], prgs[3]->prg.nodes[7], prgs[3]->prg.nodes[9] };
    EXPECT_ITERABLE_EQ(std::vector<LocalNodePtr>, vcf_ref_paths[3], exp_path3);
}