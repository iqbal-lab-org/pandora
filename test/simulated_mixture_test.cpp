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
    string prgfile = "../test/test_cases/mixtures/GC00000005_4.prg.fa";
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

    // read in the coverage information from the readfile
    //string readfile = "../test/test_cases/mixtures/read_mixtures/gene1gene5_50.50_300x.fa";
    string readfile = "../test/test_cases/mixtures/read_mixtures/gene1gene2_50.50_300x.fa";
    //string readfile = "../test/test_cases/mixtures/read_mixtures/gene2gene3_50.50_300x.fa";

    pangraph_from_read_file(readfile, mhs, pangraph, idx, prgs, w, k, max_diff, e_rate, min_cluster_size, genome_size, false);
    update_localPRGs_with_hits(pangraph, prgs);  

    // sanity check
    ASSERT_GE((uint)1,pangraph->nodes.size());
    EXPECT_TRUE((pangraph->nodes[0]->kmer_prg.covgs.size() >= 200) && (pangraph->nodes[0]->kmer_prg.covgs.size() <= 300)); // we have the gene in all 300 reads, 
															// so expect to find it ~200-300 times
    pangraph->nodes[0]->kmer_prg.save("../test/test_cases/simulated_mixtures.kg.gfa");
    prgs[0]->prg.write_gfa("../test/test_cases/simulated_mixtures.lg.gfa");

    // choose a sensible parameter for min_covg
    vector<KmerNodePtr> kmp;
    kmp.reserve(800);
    uint8_t min_covg = 30;
    pangraph->nodes[0]->kmer_prg.set_p(0.11);
    pangraph->nodes[0]->kmer_prg.num_reads = pangraph->nodes[0]->kmer_prg.covgs.size();
    pangraph->nodes[0]->kmer_prg.find_max_path(kmp);
    cout << "Update min covg from default " << +min_covg;
    for (auto n : kmp)
    {
        min_covg = (uint8_t)min((uint)min_covg, (n->covg[0]+n->covg[1]));
    }
    min_covg = min_covg/2;
    cout << " to " << +min_covg << endl;

    // find the compatible paths to use as input
    vector<deque<KmerNodePtr>> paths;
    vector<vector<pair<uint16_t,uint16_t>>> hit_pairs;
    pangraph->nodes[0]->kmer_prg.remove_shortcut_edges();
    pangraph->nodes[0]->kmer_prg.find_all_compatible_paths(paths, hit_pairs, min_covg);

    EXPECT_GE(paths.size(), (uint)2);

    // compare the sequences chosen with the truth
    // NB these are the first 5 sequences in aligned gene file
    // truth1 is very different from truth2 and truth3
    // truth1 is very similar but not identical to truth4
    // truth1 and truth5 are similar, but further apart than truth 1 and 4
    string truth1 = "ATGACTCAGAAAAATTTCGTTGAACTGCGCAACGTCACTAAACGATTTGGCAGTAATACGGTAATCGACAATATCAACCTCACCATCCCGCAGGGGCAAATGGTGACGCTGCTCGGCCCGTCCGGCTGCGGCAAAACCACTATTTTGCGCCTGGTTGCCGGGCTGGAAAAACCGAGCGAAGGGCAAATTTTCATTGATGGCGAAGACGTCACCCATCGCTCTATTCAGCAGCGCGATATCTGTATGGTGTTTCAGTCCTATGCCCTGTTCCCGCATATGTCGCTGGGAGAGAATGTCGGTTATGGCCTGAAAATGCTCGGCGTACCGCGCGCAGAGCTGAAAGCCCGCGTCAAAGAGGCGTTGGCGATGGTGGATCTGGAAGGATTCGAAGACCGCTTTGTCGATCAGATCTCCGGCGGGCAGCAGCAGCAGCGCGTGGCGCTGGCCCGCGCGCTGATCCTCAAGCCGAAAGTGCTGCTGTTTGATGAGCCGTTGAGTAACCTCGACGCCAACCTGCGTCGCAGCATGCGCGACAAGATCCGCGAGTTGCAAAAGCAGTTTGATATCACCTCGCTGTACGTCACCCACGATCAGAGCGAAGCCTTTGCGGTTTCTGATACTGTGCTGGTGATGAACAAGGGACACATCATGCAGATCGGCTCACCGCAGGATCTTTACCGCCAGCCCGCCTCCCGCTTTATGGCGAGCTTTATGGGCGATGCCAACCTGTTCCCGGCAACCTTCAGCGACGGATACGTTGATATCTACGGCTATCATCTGCCGCGCCCGCTGCACTTTGGTACACAGGGTGAAGGGATGGTCGGTGTGCGCCCGGAAGCGATCACGCTCAGCGATCGCGGCGAAGAGAGCCAGCGCTGCGTGATCCGCCATGTCGCCTATATGGGGCCGCAGTATGAAGTGACGGTGGAATGGCACGGGCAGGAGATATTATTGCAGGTCAACGCTACGCGTCTGCAACCGGACGTCGGCGAGCAGTATTATCTTGAAATCCATCCGTACGGCATGTTTGTTCTGGCGGATGCGGCA";
    string truth2 = "ATGAGTCAGAAAAATTTTGTTGAACTGCGCAACGTCACTAAACGATTCGGCAGTAATACGGTTATCGACAATATCAATCTCACCATCCCACAAGGGCAAATGGTGACGCTGCTCGGTCCTTCCGGCTGTGGCAAAACCACCATTTTGCGTCTGGTTGCCGGGCTGGAAAAACCGAGCGAAGGGCAAATTTTTATTGATGGCGAAGATGTCACGCATCGTTCCATTCAACAGCGCGATATCTGCATGGTGTTTCAGTCATACGCTCTGTTCCCGCATATGTCGCTGGGTGAAAACGTTGGCTACGGGCTGAAGATGCTTGGCGTGTCGCGTAGCGAAGTGAAACAGCGGGTGAAGGAGGCGCTGGCAATGGTGGATCTGGAAGGCTTCGAGGACCGCTATGTCGATCAGATTTCCGGTGGTCAGCAACAGCGTGTGGCACTGGCCCGCGCGTTGATCCTCAAGCCAAAGGTGCTGCTGTTTGATGAGCCGTTAAGTAACCTCGATGCCAACCTGCGCCGCAGTATGCGCGATAAGATCCGCGAGCTGCAAAAGCAGTTTAATATCACGTCGCTCTACGTCACTCACGATCAAAGTGAAGCTTTCGCGGTGTCCGACACTGTGCTGGTAATGAATAAAGGTCACATCATGCAGATTGGCTCACCGCAGGATCTCTATCGTCAGCCAGCCTCCCGATTTATGGCAAGTTTTATGGGCGACGCCAACCTGTTCCCGGCGAACTTTAGCGAAGAGTATGTCGATATCTACGGTTATCGCCTGCCGCGCGCGGCGCATTTCCCGGCGCAAGGTAGCGGCACCGTCGGCGTTCGCCCGGAAGCCATCACCTTAAGCAATCACGGCGAAGAGAGTCAGCGTTGTGTTATTCGCCATGTCGCCTACATGGGGCCGCAGTACGAAGTGACCGTAGAGTGGCATGGACAGGAGATTTTATTACAAGTAAACGCCACCCGTTTACAGCCCGATATTGGTGAGCACTATTACCTCGAAATCCATCCTTACGGGATGTTTGTTTTAGCGGATGCGGCA";
    string truth3 = "ATGAGTCAGAAAAATTTTGTGGAACTGCGCAACGTCACTAAACGATTCGGCAGTAATACGGTTATCGACAATATCAATCTCACCATCCCACAAGGGCAAATGGTGACGCTGCTTGGCCCTTCCGGCTGCGGCAAAACCACCATTTTGCGTCTGGTTGCCGGGCTGGAAAAACCGAGTGAAGGGCAAATCTTTATTGATGGCGAAGATGTCACGCATCGTTCCATTCAGCAGCGCGATATCTGCATGGTGTTTCAGTCATACGCTCTGTTCCCGCATATGTCGCTGGGCGAAAACGTCGGCTACGGGTTAAAGATGCTCGGCGTGTCGCGTAGCGAAGTAAAGCAGCGGGTGAAAGAGGCGCTGGCAATGGTGGATCTGGAAGGGTTCGAGGACCGCTATGTCGACCAGATTTCCGGTGGTCAGCAACAGCGCGTGGCGCTGGCCCGCGCGTTGATCCTCAAACCAAAGGTGCTGCTGTTTGATGAGCCGTTAAGTAACCTCGATGCCAACCTGCGCCGCAGCATGCGCGATAAGATCCGCGAGCTGCAAAAGCAGTTTAATATCACGTCGCTCTACGTGACTCACGATCAAAGTGAGGCCTTCGCGGTTTCCGACACTGTGCTGGTGATGAATAAAGGTCACATCATGCAGATTGGCTCACCGCAGGATCTCTATCGTCAGCCAGCCTCCCGCTTTATGGCAAGTTTTATGGGCGACGCCAACCTGTTCCCGGCGAACTTTAGCGAAGAGTATGTCGATATCTACGGTTATCGCCTGCCGCGCGCGGCGCATTTCCCGGTGCAAGGTAGCGGCACCGTCGGCGTTCGCCCGGAAGCCATCACCTTAAGCAATCACGGCGAAGAGAGTCAGCGTTGTGTTATTCGCCATGTCGCCTACATGGGGCCGCAGTACGAAGTGACCGTAGAGTGGCATGGACAGGAGATTTTATTACAAGTAAACGCCACCCGTTTACAGCCCGATATTGGTGAGCACTATTACCTCGAAATCCATCCTTACGGGATGTTTGTTTTAGCGGATGCGGCA";
    string truth4 = "ATGACTCAGAAAAATTTCGTTGAACTGCGCAACGTCACTAAACGATTTGGCAGTAATACGGTAATCGACAATATCAACCTCACCATCCCGCAGGGGCAAATGGTGACGCTGCTCGGCCCGTCCGGCTGCGGCAAAACCACTATTTTGCGCCTGGTTGCCGGGCTGGAAAAACCGAGCGAAGGGCAAATTTTCATTGATGGCGAAGACGTCACCCATCGCTCTATTCAGCAGCGCGATATCTGTATGGTGTTTCAGTCCTATGCCCTGTTCCCGCATATGTCGCTGGGAGAGAATGTCGGTTATGGCCTGAAAATGCTCGGCGTACCGCGCGCAGAGCTGAAAGCCCGCGTCAAAGAGGCGTTGGCGATGGTGGATCTGGAAGGATTCGAAGACCGCTTTGTCGATCAGATCTCCGGCGGGCAGCAGCAGCGCGTGGCGCTGGCCCGCGCGCTGATCCTCAAGCCGAAAGTGCTGCTGTTTGATGAGCCGTTGAGTAACCTCGACGCCAACCTGCGTCGCAGCATGCGCGACAAGATCCGCGAGTTGCAAAAGCAGTTTGATATCACCTCGCTGTACGTCACCCACGATCAGAGCGAAGCCTTTGCGGTTTCTGATACTGTGCTGGTGATGAACAAGGGACACATCATGCAGATCGGCTCACCGCAGGATCTTTACCGCCAGCCCGCCTCCCGCTTTATGGCGAGCTTTATGGGCGATGCCAACCTGTTCCCGGCAACCTTCAGCGACGGATACGTTGATATCTACGGCTATCATCTGCCGCGCCCGCTGCACTTTGGTACACAGGGTGAAGGGATGGTCGGTGTGCGCCCGGAAGCGATCACGCTCAGCGATCGCGGCGAAGAGAGCCAGCGCTGCGTGATCCGCCATGTCGCCTATATGGGGCCGCAGTATGAAGTGACGGTGGAATGGCACGGGCAGGAGATATTATTGCAGGTCAACGCTACGCGTCTGCAACCGGACGTCGGCGAGCAGTATTATCTTGAAATCCATCCGTACGGCATGTTTGTTCTGGCGGATGCGGCA";
    string truth5 = "ATGACTCAGAAAAATTTCGTTGAACTGCGCAACGTCACTAAACGATTTGGCAGTAATACGGTAATCGACAATATCAACCTCACCATCCCGCAGGGGCAAATGGTGACGCTGTTCGGCCCGTCCGGCTGCGGCAAAACCACTATTTTGCGCCTGGTTGCCGGGCTGGAAAAACCGAGCGAAGGGCAAATTTTCATTGATGGCGAAGACGTCACCCATCGCTCTATTCAGCAGCGCGATATCTGTATGGTGTTTCAGTCCTATGCCCTGTTCCCGCATATGTCGCTGGGAGAGAATGTCGGTTATGGCCTGAAAATGCTCGGCGTACCGCGCGCAGAGCTGAAAGCCCGCGTCAAAGAGGCGTTGGCGATGGTGGATCTGGAAGGATTCGAAGACCGCTTTGTCGATCAGATCTCCGGCGGGCAGCAGCAGCGCGTGGCGCTGGCCCGCGCGCTGATCCTCAAGCCGAAAGTGCTGCTGTTTGATGAGCCGTTGAGTAACCTCGACGCCAACCTGCGTCGCAGCATGCGCGACAAGATCCGCGAGTTGCAAAAGCAGTTTGATATCACCTCGCTGTACGTCACCCACGATCAGAGCGAAGCCTTTGCGGTTTCTGATACTGTGCTGGTGATGAACAAGGGACACATCATGCAGATCGGCTCACCGCAGGATCTTTACCGCCAGCCCGCCTCCCGCTTTATGGCGAGCTTTATGGGCGATGCCAACCTGTTCCCGGCAACCTTCAGCGACGGATACGTTGATATCTACGGCTATCATCTGCCGCGCCCGCTGCACTTTGGTACACAGGGTGAAGGGATGGTCGGTGTGCGCCCGGAAGCGATCACGCTCAGCGATCGCGGCGAAGAGAGCCAGCGCTGCGTGATCCGCCATGTCGCCTATATGGGGCCGCAGTATGAAGTGACGGTGGAATGGCACGGGCAGGAGATATTATTGCAGGTCAACGCTACGCGTCTGCAACCGGACGTCGGCGAGCAGTATTATCTTGAAATCCATCCGTACGGCATGTTTGTTCTGGCGGATGCGGCA";
    truth1.erase(truth1.length()-15, truth1.length()-15);
    truth1.erase(0, 15);
    //cout << truth1 << endl;
    truth2.erase(truth2.length()-15, truth2.length()-15);
    truth2.erase(0, 15);
    //cout << truth2 << endl;
    truth3.erase(truth3.length()-15, truth3.length()-15);
    truth3.erase(0, 15);
    //cout << truth3 << endl;
    truth4.erase(truth4.length()-15, truth4.length()-15);
    truth4.erase(0, 15);
    //cout << truth4 << endl;
    truth5.erase(truth5.length()-15, truth5.length()-15);
    truth5.erase(0, 15);
    //cout << truth5 << endl;

    vector<LocalNodePtr> lmp;
    lmp.reserve(100);
    string result, result1, result2, result3, result4, result5;
    uint found = 0;

    // check the true paths are in the input compatible paths
    uint count = 0;
    for (auto p : paths)
    {
        cout << ".";
	    kmp = vector<KmerNodePtr>(p.begin(), p.end());
        for (auto n : kmp)
        {
            cout << n->path << " ";
        }
        cout << endl;
    	lmp = prgs[0]->localnode_path_from_kmernode_path(kmp, w);
	    result = "";
	    for (auto n : lmp)
    	{
            cout << *n ";
            result += n->seq;
    	}
        cout << endl;
        result1 = result.substr(15,truth1.length());
        result2 = result.substr(15,truth2.length());
        result3 = result.substr(15,truth3.length());
        result4 = result.substr(15,truth4.length());
        result5 = result.substr(15,truth5.length());
	    if ((result1 == truth1)
            or (result2 == truth2)
            or (result3 == truth3)
            or (result4 == truth4)
            or (result5 == truth5))
	    {
            cout << endl << (result1 == truth1) << (result2 == truth2) << (result3 == truth3) << (result4 == truth4) << (result5 == truth5) << endl;
            cout << result << endl;
            for (auto l : kmp)
            {
                cout << l->id << " ";
            }
            /*for (auto l : lmp)
            {
                cout << *l << " ";
            }*/
            cout << endl;

            found += 1;
	    }
        count += 1
	    if (count == 20)
	    {
	        break;
	    }
    }
    cout << "done" << endl;
    EXPECT_LE((uint)2, found);

    // now run EM
    double eps = 50;
    PathAbundanceEstimator pae(hit_pairs, paths, 1e-8, 1000);
    std::vector<double> pathCnts = pae.runEM();
    EXPECT_EQ(pathCnts.size(), paths.size());

    // check the best 2 paths at end of EM
    cout << "check top path";
    uint path_i = distance(pathCnts.begin(), max_element(pathCnts.begin(),pathCnts.end()));
    EXPECT_LE(std::abs(pathCnts[path_i]-static_cast<double>(pangraph->nodes[0]->kmer_prg.covgs.size()/2)), eps);

    kmp = vector<KmerNodePtr>(paths[path_i].begin(), paths[path_i].end());
    lmp = prgs[0]->localnode_path_from_kmernode_path(kmp, w); 
    result = "";
    for (auto n : lmp)
    {
        result += n->seq;
    }
    result1 = result.substr(15,truth1.length());
    result2 = result.substr(15,truth2.length());
    result3 = result.substr(15,truth3.length());
    result4 = result.substr(15,truth4.length());
    result5 = result.substr(15,truth5.length());
    cout << " gives " << endl;
    EXPECT_EQ((result1==truth1) or (result2==truth2) or (result3==truth3) or (result4==truth4) or (result5==truth5), true);
    cout << endl << (result1 == truth1) << (result2 == truth2) << (result3 == truth3) << (result4 == truth4) << (result5 == truth5) << endl;
    cout << result << endl;

    cout << "check next path";
    uint next_i = max(distance(pathCnts.begin(), max_element(pathCnts.begin(),pathCnts.begin()+path_i)),
                  distance(pathCnts.begin(), max_element(pathCnts.begin()+path_i+1,pathCnts.end())));
    EXPECT_LE(std::abs(pathCnts[next_i]-static_cast<double>(pangraph->nodes[0]->kmer_prg.covgs.size()/2)), eps);

    kmp = vector<KmerNodePtr>(paths[next_i].begin(), paths[next_i].end());
    lmp = prgs[0]->localnode_path_from_kmernode_path(kmp, w);
    result = "";
    for (auto n : lmp)
    {   
        result += n->seq;
    }
    result1 = result.substr(15,truth1.length());
    result2 = result.substr(15,truth2.length());
    result3 = result.substr(15,truth3.length());
    result4 = result.substr(15,truth4.length());
    result5 = result.substr(15,truth5.length());
    cout << " gives " << endl;
    EXPECT_EQ((result1==truth1) or (result2==truth2) or (result3==truth3) or (result4==truth4) or (result5==truth5), true);
    cout << endl << (result1 == truth1) << (result2 == truth2) << (result3 == truth3) << (result4 == truth4) << (result5 == truth5) << endl;
    cout << result << endl;

    // clear up
    idx->clear();
    delete idx;
    delete mhs;
    delete pangraph;
}

