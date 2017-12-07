#include "gtest/gtest.h"
#include "test_macro.cpp"
#include <cmath>
#include "index.h"
#include "pangraph.h"
#include "minihits.h"
#include "utils.h"
#include "kmergraph.h"
#include "pathAbundanceEstimator.h"

using namespace std;

class SimulatedMixtureTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    /*string prgfile = "../test/test_cases/mixtures/GC00000005_4.prg.fa";
    uint32_t w = 14, k = 15, min_cluster_size = 10, genome_size = 5000000; // default parameters
    int max_diff = 500;
    float e_rate = 0.11;
    Index *idx;
    idx = new Index();
    idx->load(prgfile, w, k);
    vector<LocalPRG *> prgs;
    read_prg_file(prgs, prgfile);
    load_PRG_kmergraphs(prgs, w, k, prgfile);
    MinimizerHits *mhs;
    mhs = new MinimizerHits(100 * idx->minhash.size());
    PanGraph *pangraph;
    pangraph = new PanGraph();*/
  }
  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
    /*idx->clear();
    delete idx;
    delete mhs;
    delete pangraph;*/
  }
};

// These tests run a simplified version of pandora_map

TEST_F(SimulatedMixtureTest, gene1gene2_5050) {
/*    string prgfile = "../test/test_cases/mixtures/GC00000005_4.prg.fa";
    uint32_t w = 14, k = 15, min_cluster_size = 10, genome_size = 5000000; // default parameters
    int max_diff = 500;
    float e_rate = 0.11;
    Index *idx;
    idx = new Index();
    idx->load(prgfile, w, k);
    vector<LocalPRG *> prgs;
    read_prg_file(prgs, prgfile);
    load_PRG_kmergraphs(prgs, w, k, prgfile);
    MinimizerHits *mhs;
    mhs = new MinimizerHits(100 * idx->minhash.size());
    PanGraph *pangraph;
    pangraph = new PanGraph();

    string readfile = "../test/test_cases/mixtures/read_mixtures/gene1gene2_50.50_300x.fa"; 
    pangraph_from_read_file(readfile, mhs, pangraph, idx, prgs, w, k, max_diff, e_rate, min_cluster_size, genome_size, false);
    update_localPRGs_with_hits(pangraph, prgs);  

    ASSERT_GE((uint)1,pangraph->nodes.size());
    EXPECT_TRUE((pangraph->nodes[0]->kmer_prg.covgs.size() >= 200) && (pangraph->nodes[0]->kmer_prg.covgs.size() <= 300)); // we have the gene in all 300 reads, 
															// so expect to find it ~200-300 times

    vector<deque<KmerNodePtr>> paths;
    vector<vector<pair<uint16_t,uint16_t>>> hit_pairs;
    pangraph->nodes[0]->kmer_prg.find_all_compatible_paths(paths, hit_pairs);

    EXPECT_EQ((uint)2, paths.size());

    double eps = 10;
    PathAbundanceEstimator pae(hit_pairs, paths, 1e-8, 1000);
    std::vector<double> pathCnts = pae.runEM();

    EXPECT_EQ(pathCnts.size(), static_cast<size_t>(2));
    EXPECT_LE(std::abs(pathCnts[0]-static_cast<double>(pangraph->nodes[0]->kmer_prg.covgs.size()/2)), eps);
    EXPECT_LE(std::abs(pathCnts[1]-static_cast<double>(pangraph->nodes[0]->kmer_prg.covgs.size()/2)), eps);

    idx->clear();
    delete idx;
    delete mhs;
    delete pangraph;*/
}

