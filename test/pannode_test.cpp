#include <cstdint>
#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "pangenome/ns.cpp"
#include "pangenome/pannode.h"
#include "pangenome/pansample.h"
#include "pangenome/pangraph.h"
#include "pangenome/panread.h"
#include "minihit.h"
#include "localPRG.h"


using namespace pangenome;

TEST(PangenomeNodeTest, create) {

    pangenome::Node pan_node(4, 3, "3");
    uint32_t j = 3;
    EXPECT_EQ(j, pan_node.node_id);
    EXPECT_EQ((uint) 4, pan_node.prg_id);
    EXPECT_EQ("3", pan_node.name);
    EXPECT_EQ((uint) 1, pan_node.covg);
    EXPECT_EQ((uint) 0, pan_node.reads.size());
    EXPECT_EQ((uint) 0, pan_node.samples.size());
}

TEST(PangenomeNodeTest, get_name) {
    pangenome::Node pn1(3, 3, "3");
    pangenome::Node pn2(2, 2, "2");
    pangenome::Node pn3(2, 4, "2");

    EXPECT_EQ(pn1.get_name(), "3");
    EXPECT_EQ(pn2.get_name(), "2");
    EXPECT_EQ(pn3.get_name(), "2.4");
}

TEST(PangenomeNodeTest, add_path) {
    pangenome::Node pn1(3, 3, "3");
    std::vector<KmerNodePtr> kmp;
    pn1.add_path(kmp, 0);

    KmerGraph kg;
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
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
    EXPECT_EQ((uint) 7, kg.nodes.size());

    pn1.kmer_prg = kg;
    EXPECT_EQ((uint) 7, pn1.kmer_prg.nodes.size());
    pn1.kmer_prg.sort_topologically();
    EXPECT_EQ((uint) 7, pn1.kmer_prg.sorted_nodes.size());
    kmp = {pn1.kmer_prg.sorted_nodes[0], pn1.kmer_prg.sorted_nodes[3], pn1.kmer_prg.sorted_nodes[4],
           pn1.kmer_prg.sorted_nodes[6]};
    pn1.add_path(kmp, 0);
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[0]->get_covg(0, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[1]->get_covg(0, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[2]->get_covg(0, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[3]->get_covg(0, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[4]->get_covg(0, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[5]->get_covg(0, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[6]->get_covg(0, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[0]->get_covg(1, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[1]->get_covg(1, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[2]->get_covg(1, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[3]->get_covg(1, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[4]->get_covg(1, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[5]->get_covg(1, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[6]->get_covg(1, 0));
}

TEST(PangenomeNodeTest, get_read_overlap_coordinates) {
    pangenome::Node pan_node(3, 3, "3");
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
    pr = std::make_shared<pangenome::Read>(1);
    pr->add_hits(3, mhits.hits);
    pan_node.reads.insert(pr);
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
    pr = std::make_shared<pangenome::Read>(2);
    pr->add_hits(3, mhits.hits);
    pan_node.reads.insert(pr);
    mhits.clear();

    delete mr;

    std::vector<std::vector<uint32_t>> read_overlap_coordinates;
    pan_node.get_read_overlap_coordinates(read_overlap_coordinates);
    std::vector<std::vector<uint32_t>> expected_read_overlap_coordinates = {{1, 0, 6,  1},
                                                                            {2, 2, 10, 0}};
    for (const auto &coord : read_overlap_coordinates) {
        if (coord[0] == 1) {
            EXPECT_ITERABLE_EQ(std::vector<uint32_t>, expected_read_overlap_coordinates[0], coord);
        } else {
            EXPECT_ITERABLE_EQ(std::vector<uint32_t>, expected_read_overlap_coordinates[1], coord);
        }
    }
}

TEST(PangenomeNodeTest,construct_multisample_vcf_single_prg)
{
    uint32_t prg_id = 3, w=1, k=3, min_kmer_covg=0;
    std::string prg_name = "nested varsite";
    LocalPRG local_prg(prg_id, prg_name , "A 5 G 7 C 8 T 8 CT 7  6 G 5 T");

    auto index = std::make_shared<Index>();
    local_prg.minimizer_sketch(index, w, k);
    auto prg_ptr = std::make_shared<LocalPRG>(local_prg);
    auto &kg = local_prg.kmer_prg;

    pangenome::Graph pangraph;

    //sample1
    std::string sample_name = "sample1";
    std::vector<KmerNodePtr> sample_kmer_path = {kg.nodes[0], kg.nodes[2], kg.nodes[6], kg.nodes[9]};
    pangraph.add_node(prg_id, prg_name, sample_name, 0, prg_ptr, sample_kmer_path);

    //sample2 identical to sample1
    sample_name = "sample2";
    sample_kmer_path = {kg.nodes[0], kg.nodes[2], kg.nodes[6], kg.nodes[9]};
    pangraph.add_node(prg_id, prg_name, sample_name, 1, prg_ptr, sample_kmer_path);

    //sample3 with top path
    sample_name = "sample3";
    sample_kmer_path = {kg.nodes[0], kg.nodes[1], kg.nodes[5], kg.nodes[9]};
    pangraph.add_node(prg_id, prg_name, sample_name, 2, prg_ptr, sample_kmer_path);

    //sample4 with bottom path
    sample_name = "sample4";
    sample_kmer_path = {kg.nodes[0], kg.nodes[4], kg.nodes[9]};
    pangraph.add_node(prg_id, prg_name, sample_name, 3, prg_ptr, sample_kmer_path);

    LocalPRG dummy(0,"null","");
    auto dummy_prg_ptr = std::make_shared<LocalPRG>(dummy);
    pangraph.setup_kmergraphs({dummy_prg_ptr, dummy_prg_ptr, dummy_prg_ptr, prg_ptr}, 4);

    VCF master_vcf;
    std::vector<LocalNodePtr> vcf_reference_path = {local_prg.prg.nodes[0], local_prg.prg.nodes[1], local_prg.prg.nodes[3], local_prg.prg.nodes[5], local_prg.prg.nodes[7]};
    auto &pannode = *pangraph.nodes[prg_id];
    pannode.construct_multisample_vcf(master_vcf, vcf_reference_path, prg_ptr, w, min_kmer_covg);

    EXPECT_EQ((uint)2, master_vcf.records.size());
    EXPECT_EQ((uint)4, master_vcf.samples.size());

    //NB samples order changes to get index of each sample so can compare
    //samples 1 and 2 are ref, sample 3 is top path and sample 4 is bottom path
    auto iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample1");
    auto sample1_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample2");
    auto sample2_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample3");
    auto sample3_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample4");
    auto sample4_index = std::distance(master_vcf.samples.begin(), iter);
    std::vector<uint16_t> alt_gt = {1};
    std::vector<uint16_t> ref_gt = {0};
    std::vector<uint16_t> no_covg2 = {0,0};
    std::vector<uint16_t> no_covg3 = {0,0,0};

    EXPECT_EQ((uint)1, master_vcf.records[0].pos);
    EXPECT_EQ("GT", master_vcf.records[0].ref);
    EXPECT_EQ((uint)1, master_vcf.records[0].alt.size());
    EXPECT_EQ("G", master_vcf.records[0].alt[0]);
    EXPECT_EQ((uint)4, master_vcf.records[0].samples.size());
    EXPECT_FALSE(master_vcf.records[0].samples[sample4_index].find("GT") == master_vcf.records[0].samples[sample4_index].end()) ;
    EXPECT_TRUE(master_vcf.records[0].samples[sample3_index].find("GT") == master_vcf.records[0].samples[sample3_index].end()) ;
    EXPECT_FALSE(master_vcf.records[0].samples[sample2_index].find("GT") == master_vcf.records[0].samples[sample2_index].end()) ;
    EXPECT_FALSE(master_vcf.records[0].samples[sample1_index].find("GT") == master_vcf.records[0].samples[sample1_index].end()) ;
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample4_index]["GT"], alt_gt);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample2_index]["GT"], ref_gt);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample1_index]["GT"], ref_gt);
    std::vector<std::string> formats = {"MEAN_FWD_COVG", "MEAN_REV_COVG", "MED_FWD_COVG", "MED_REV_COVG",
                              "SUM_FWD_COVG", "SUM_REV_COVG"};
    for (const auto format : formats){
        EXPECT_FALSE(master_vcf.records[0].samples[sample1_index].find(format) == master_vcf.records[0].samples[sample1_index].end());
        EXPECT_FALSE(master_vcf.records[0].samples[sample2_index].find(format) == master_vcf.records[0].samples[sample2_index].end());
        EXPECT_FALSE(master_vcf.records[0].samples[sample3_index].find(format) == master_vcf.records[0].samples[sample3_index].end());
        EXPECT_FALSE(master_vcf.records[0].samples[sample4_index].find(format) == master_vcf.records[0].samples[sample4_index].end());
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample1_index][format], no_covg2);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample2_index][format], no_covg2);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample3_index][format], no_covg2);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample4_index][format], no_covg2);
    }


    EXPECT_EQ((uint)2, master_vcf.records[1].pos);
    EXPECT_EQ("T", master_vcf.records[1].ref);
    EXPECT_EQ((uint)2, master_vcf.records[1].alt.size());
    EXPECT_EQ("C", master_vcf.records[1].alt[0]);
    EXPECT_EQ("CT", master_vcf.records[1].alt[1]);
    EXPECT_EQ((uint)4, master_vcf.records[0].samples.size());
    EXPECT_TRUE(master_vcf.records[1].samples[sample4_index].find("GT") == master_vcf.records[1].samples[sample4_index].end()) ;
    EXPECT_FALSE(master_vcf.records[1].samples[sample3_index].find("GT") == master_vcf.records[1].samples[sample3_index].end()) ;
    EXPECT_FALSE(master_vcf.records[1].samples[sample2_index].find("GT") == master_vcf.records[1].samples[sample2_index].end()) ;
    EXPECT_FALSE(master_vcf.records[1].samples[sample1_index].find("GT") == master_vcf.records[1].samples[sample1_index].end()) ;
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample3_index]["GT"], alt_gt);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample2_index]["GT"], ref_gt);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample1_index]["GT"], ref_gt);

    for (const auto format : formats){
        EXPECT_FALSE(master_vcf.records[1].samples[sample1_index].find(format) == master_vcf.records[1].samples[sample1_index].end());
        EXPECT_FALSE(master_vcf.records[1].samples[sample2_index].find(format) == master_vcf.records[1].samples[sample2_index].end());
        EXPECT_FALSE(master_vcf.records[1].samples[sample3_index].find(format) == master_vcf.records[1].samples[sample3_index].end());
        EXPECT_FALSE(master_vcf.records[1].samples[sample4_index].find(format) == master_vcf.records[1].samples[sample4_index].end());
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample1_index][format], no_covg3);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample2_index][format], no_covg3);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample3_index][format], no_covg3);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample4_index][format], no_covg3);
    }
}

TEST(PangenomeNodeTest,construct_multisample_vcf_two_prg)
{
    uint32_t prg_id1 = 3, w=1, k=3, min_kmer_covg=0;
    std::string prg_name1 = "nested varsite";
    LocalPRG local_prg1(prg_id1, prg_name1 , "A 5 G 7 C 8 T 8 CT 7  6 G 5 T");
    uint32_t prg_id2 = 5;
    std::string prg_name2 = "modified";
    LocalPRG local_prg2(prg_id2, prg_name2 , "A 5 G 7 G 8 A 8 GA 7  6 G 5 T");

    auto index = std::make_shared<Index>();
    local_prg1.minimizer_sketch(index, w, k);
    auto prg_ptr1 = std::make_shared<LocalPRG>(local_prg1);
    auto &kg1 = local_prg1.kmer_prg;
    local_prg2.minimizer_sketch(index, w, k);
    auto prg_ptr2 = std::make_shared<LocalPRG>(local_prg2);
    auto &kg2 = local_prg2.kmer_prg;

    pangenome::Graph pangraph;

    //sample1
    std::string sample_name = "sample1";
    std::vector<KmerNodePtr> sample_kmer_path = {kg1.nodes[0], kg1.nodes[2], kg1.nodes[6], kg1.nodes[9]};
    pangraph.add_node(prg_id1, prg_name1, sample_name, 0, prg_ptr1, sample_kmer_path);
    sample_kmer_path = {kg2.nodes[0], kg2.nodes[1], kg2.nodes[5], kg2.nodes[9]};
    pangraph.add_node(prg_id2, prg_name2, sample_name, 0, prg_ptr2, sample_kmer_path);

    //sample2 identical to sample1 in prg1, no prg2
    sample_name = "sample2";
    sample_kmer_path = {kg1.nodes[0], kg1.nodes[2], kg1.nodes[6], kg1.nodes[9]};
    pangraph.add_node(prg_id1, prg_name1, sample_name, 1, prg_ptr1, sample_kmer_path);

    //sample3 with top path
    sample_name = "sample3";
    sample_kmer_path = {kg1.nodes[0], kg1.nodes[1], kg1.nodes[5], kg1.nodes[9]};
    pangraph.add_node(prg_id1, prg_name1, sample_name, 2, prg_ptr1, sample_kmer_path);
    sample_kmer_path = {kg2.nodes[0], kg2.nodes[4], kg2.nodes[9]};
    pangraph.add_node(prg_id2, prg_name2, sample_name, 2, prg_ptr2, sample_kmer_path);

    //sample4 with bottom path
    sample_name = "sample4";
    sample_kmer_path = {kg1.nodes[0], kg1.nodes[4], kg1.nodes[9]};
    pangraph.add_node(prg_id1, prg_name1, sample_name, 3, prg_ptr1, sample_kmer_path);
    sample_kmer_path = {kg2.nodes[0], kg2.nodes[3], kg2.nodes[7], kg2.nodes[8], kg2.nodes[9]};
    pangraph.add_node(prg_id2, prg_name2, sample_name, 3, prg_ptr2, sample_kmer_path);

    LocalPRG dummy(0,"null","");
    auto dummy_prg_ptr = std::make_shared<LocalPRG>(dummy);
    pangraph.setup_kmergraphs({dummy_prg_ptr, dummy_prg_ptr, dummy_prg_ptr, prg_ptr1, dummy_prg_ptr, prg_ptr2}, 4);

    VCF master_vcf;
    std::vector<LocalNodePtr> vcf_reference_path1 = {local_prg1.prg.nodes[0], local_prg1.prg.nodes[1], local_prg1.prg.nodes[3], local_prg1.prg.nodes[5], local_prg1.prg.nodes[7]};
    std::vector<LocalNodePtr> vcf_reference_path2 = {local_prg2.prg.nodes[0], local_prg2.prg.nodes[1], local_prg2.prg.nodes[3], local_prg2.prg.nodes[5], local_prg2.prg.nodes[7]};

    auto &pannode1 = *pangraph.nodes[prg_id1];
    pannode1.construct_multisample_vcf(master_vcf, vcf_reference_path1, prg_ptr1, w, min_kmer_covg);
    auto &pannode2 = *pangraph.nodes[prg_id2];
    pannode2.construct_multisample_vcf(master_vcf, vcf_reference_path2, prg_ptr2, w, min_kmer_covg);

    EXPECT_EQ((uint)4, master_vcf.records.size());
    EXPECT_EQ((uint)4, master_vcf.samples.size());

    //NB samples order changes to get index of each sample so can compare
    //samples 1 and 2 are ref, sample 3 is top path and sample 4 is bottom path
    auto iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample1");
    auto sample1_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample2");
    auto sample2_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample3");
    auto sample3_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample4");
    auto sample4_index = std::distance(master_vcf.samples.begin(), iter);
    std::vector<uint16_t> alt_gt = {1};
    std::vector<uint16_t> ref_gt = {0};
    std::vector<uint16_t> alt2_gt = {2};
    std::vector<uint16_t> no_covg2 = {0,0};
    std::vector<uint16_t> no_covg3 = {0,0,0};

    EXPECT_EQ((uint)1, master_vcf.records[0].pos);
    EXPECT_EQ("GT", master_vcf.records[0].ref);
    EXPECT_EQ((uint)1, master_vcf.records[0].alt.size());
    EXPECT_EQ("G", master_vcf.records[0].alt[0]);
    EXPECT_EQ((uint)4, master_vcf.records[0].samples.size());
    EXPECT_FALSE(master_vcf.records[0].samples[sample4_index].find("GT") == master_vcf.records[0].samples[sample4_index].end()) ;
    EXPECT_TRUE(master_vcf.records[0].samples[sample3_index].find("GT") == master_vcf.records[0].samples[sample3_index].end()) ;
    EXPECT_FALSE(master_vcf.records[0].samples[sample2_index].find("GT") == master_vcf.records[0].samples[sample2_index].end()) ;
    EXPECT_FALSE(master_vcf.records[0].samples[sample1_index].find("GT") == master_vcf.records[0].samples[sample1_index].end()) ;
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample4_index]["GT"], alt_gt);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample2_index]["GT"], ref_gt);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample1_index]["GT"], ref_gt);
    std::vector<std::string> formats = {"MEAN_FWD_COVG", "MEAN_REV_COVG", "MED_FWD_COVG", "MED_REV_COVG",
                                        "SUM_FWD_COVG", "SUM_REV_COVG"};
    for (const auto format : formats){
        EXPECT_FALSE(master_vcf.records[0].samples[sample1_index].find(format) == master_vcf.records[0].samples[sample1_index].end());
        EXPECT_FALSE(master_vcf.records[0].samples[sample2_index].find(format) == master_vcf.records[0].samples[sample2_index].end());
        EXPECT_FALSE(master_vcf.records[0].samples[sample3_index].find(format) == master_vcf.records[0].samples[sample3_index].end());
        EXPECT_FALSE(master_vcf.records[0].samples[sample4_index].find(format) == master_vcf.records[0].samples[sample4_index].end());
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample1_index][format], no_covg2);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample2_index][format], no_covg2);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample3_index][format], no_covg2);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample4_index][format], no_covg2);
    }


    EXPECT_EQ((uint)2, master_vcf.records[1].pos);
    EXPECT_EQ("T", master_vcf.records[1].ref);
    EXPECT_EQ((uint)2, master_vcf.records[1].alt.size());
    EXPECT_EQ("C", master_vcf.records[1].alt[0]);
    EXPECT_EQ("CT", master_vcf.records[1].alt[1]);
    EXPECT_EQ((uint)4, master_vcf.records[1].samples.size());
    EXPECT_TRUE(master_vcf.records[1].samples[sample4_index].find("GT") == master_vcf.records[1].samples[sample4_index].end()) ;
    EXPECT_FALSE(master_vcf.records[1].samples[sample3_index].find("GT") == master_vcf.records[1].samples[sample3_index].end()) ;
    EXPECT_FALSE(master_vcf.records[1].samples[sample2_index].find("GT") == master_vcf.records[1].samples[sample2_index].end()) ;
    EXPECT_FALSE(master_vcf.records[1].samples[sample1_index].find("GT") == master_vcf.records[1].samples[sample1_index].end()) ;
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample3_index]["GT"], alt_gt);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample2_index]["GT"], ref_gt);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample1_index]["GT"], ref_gt);

    for (const auto format : formats){
        EXPECT_FALSE(master_vcf.records[1].samples[sample1_index].find(format) == master_vcf.records[1].samples[sample1_index].end());
        EXPECT_FALSE(master_vcf.records[1].samples[sample2_index].find(format) == master_vcf.records[1].samples[sample2_index].end());
        EXPECT_FALSE(master_vcf.records[1].samples[sample3_index].find(format) == master_vcf.records[1].samples[sample3_index].end());
        EXPECT_FALSE(master_vcf.records[1].samples[sample4_index].find(format) == master_vcf.records[1].samples[sample4_index].end());
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample1_index][format], no_covg3);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample2_index][format], no_covg3);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample3_index][format], no_covg3);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample4_index][format], no_covg3);
    }

    EXPECT_EQ((uint)1, master_vcf.records[2].pos);
    EXPECT_EQ("GA", master_vcf.records[2].ref);
    EXPECT_EQ((uint)1, master_vcf.records[2].alt.size());
    EXPECT_EQ("G", master_vcf.records[2].alt[0]);
    EXPECT_EQ((uint)4, master_vcf.records[2].samples.size());
    EXPECT_TRUE(master_vcf.records[2].samples[sample4_index].find("GT") == master_vcf.records[2].samples[sample4_index].end()) ;
    EXPECT_FALSE(master_vcf.records[2].samples[sample3_index].find("GT") == master_vcf.records[2].samples[sample3_index].end()) ;
    EXPECT_TRUE(master_vcf.records[2].samples[sample2_index].find("GT") == master_vcf.records[2].samples[sample2_index].end()) ;
    EXPECT_TRUE(master_vcf.records[2].samples[sample1_index].find("GT") == master_vcf.records[2].samples[sample1_index].end()) ;
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[2].samples[sample3_index]["GT"], alt_gt);

    for (const auto format : formats){
        EXPECT_FALSE(master_vcf.records[2].samples[sample1_index].find(format) == master_vcf.records[2].samples[sample1_index].end());
        EXPECT_TRUE(master_vcf.records[2].samples[sample2_index].find(format) == master_vcf.records[2].samples[sample2_index].end());
        EXPECT_FALSE(master_vcf.records[2].samples[sample3_index].find(format) == master_vcf.records[2].samples[sample3_index].end());
        EXPECT_FALSE(master_vcf.records[2].samples[sample4_index].find(format) == master_vcf.records[2].samples[sample4_index].end());
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[2].samples[sample1_index][format], no_covg2);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[2].samples[sample3_index][format], no_covg2);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[2].samples[sample4_index][format], no_covg2);
    }


    EXPECT_EQ((uint)2, master_vcf.records[3].pos);
    EXPECT_EQ("A", master_vcf.records[3].ref);
    EXPECT_EQ((uint)2, master_vcf.records[3].alt.size());
    EXPECT_EQ("G", master_vcf.records[3].alt[0]);
    EXPECT_EQ("GA", master_vcf.records[3].alt[1]);
    EXPECT_EQ((uint)4, master_vcf.records[3].samples.size());
    EXPECT_FALSE(master_vcf.records[3].samples[sample4_index].find("GT") == master_vcf.records[3].samples[sample4_index].end()) ;
    EXPECT_TRUE(master_vcf.records[3].samples[sample3_index].find("GT") == master_vcf.records[3].samples[sample3_index].end()) ;
    EXPECT_TRUE(master_vcf.records[3].samples[sample2_index].find("GT") == master_vcf.records[3].samples[sample2_index].end()) ;
    EXPECT_FALSE(master_vcf.records[3].samples[sample1_index].find("GT") == master_vcf.records[3].samples[sample1_index].end()) ;
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[3].samples[sample1_index]["GT"], alt_gt);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[3].samples[sample4_index]["GT"], alt2_gt);

    for (const auto format : formats){
        EXPECT_FALSE(master_vcf.records[3].samples[sample1_index].find(format) == master_vcf.records[3].samples[sample1_index].end());
        EXPECT_TRUE(master_vcf.records[3].samples[sample2_index].find(format) == master_vcf.records[3].samples[sample2_index].end());
        EXPECT_FALSE(master_vcf.records[3].samples[sample3_index].find(format) == master_vcf.records[3].samples[sample3_index].end());
        EXPECT_FALSE(master_vcf.records[3].samples[sample4_index].find(format) == master_vcf.records[3].samples[sample4_index].end());
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[3].samples[sample1_index][format], no_covg3);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[3].samples[sample3_index][format], no_covg3);
        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[3].samples[sample4_index][format], no_covg3);
    }
}

TEST(PangenomeNodeTest,construct_multisample_vcf_two_prg_with_covgs)
{
    uint32_t prg_id1 = 3, w=1, k=3, min_kmer_covg=0;
    std::string prg_name1 = "nested varsite";
    LocalPRG local_prg1(prg_id1, prg_name1 , "A 5 G 7 C 8 T 8 CT 7  6 G 5 T");
    uint32_t prg_id2 = 5;
    std::string prg_name2 = "modified";
    LocalPRG local_prg2(prg_id2, prg_name2 , "A 5 G 7 G 8 A 8 GA 7  6 G 5 T");

    auto index = std::make_shared<Index>();
    local_prg1.minimizer_sketch(index, w, k);
    auto prg_ptr1 = std::make_shared<LocalPRG>(local_prg1);
    auto &kg1 = local_prg1.kmer_prg;
    local_prg2.minimizer_sketch(index, w, k);
    auto prg_ptr2 = std::make_shared<LocalPRG>(local_prg2);
    auto &kg2 = local_prg2.kmer_prg;

    pangenome::Graph pangraph;

    //sample1
    std::string sample_name = "sample1";
    std::vector<KmerNodePtr> sample_kmer_path = {kg1.nodes[0], kg1.nodes[2], kg1.nodes[6], kg1.nodes[9]};
    pangraph.add_node(prg_id1, prg_name1, sample_name, 0, prg_ptr1, sample_kmer_path);
    sample_kmer_path = {kg2.nodes[0], kg2.nodes[1], kg2.nodes[5], kg2.nodes[9]};
    pangraph.add_node(prg_id2, prg_name2, sample_name, 0, prg_ptr2, sample_kmer_path);

    LocalPRG dummy(0,"null","");
    auto dummy_prg_ptr = std::make_shared<LocalPRG>(dummy);
    pangraph.setup_kmergraphs({dummy_prg_ptr, dummy_prg_ptr, dummy_prg_ptr, prg_ptr1, dummy_prg_ptr, prg_ptr2}, 4);

    auto &pannode1 = *pangraph.nodes[prg_id1];
    auto &pannode2 = *pangraph.nodes[prg_id2];

    pannode1.kmer_prg.nodes[0]->set_covg(4, 0, 0);
    pannode1.kmer_prg.nodes[2]->set_covg(4, 0, 0);
    pannode1.kmer_prg.nodes[6]->set_covg(4, 0, 0);
    pannode1.kmer_prg.nodes[9]->set_covg(4, 0, 0);
    pannode2.kmer_prg.nodes[0]->set_covg(4, 0, 0);
    pannode2.kmer_prg.nodes[1]->set_covg(4, 0, 0);
    pannode2.kmer_prg.nodes[5]->set_covg(4, 0, 0);
    pannode2.kmer_prg.nodes[9]->set_covg(4, 0, 0);

    //sample2 identical to sample1 in prg1, no prg2
    sample_name = "sample2";
    sample_kmer_path = {kg1.nodes[0], kg1.nodes[2], kg1.nodes[6], kg1.nodes[9]};
    pangraph.add_node(prg_id1, prg_name1, sample_name, 1, prg_ptr1, sample_kmer_path);
    pannode1.kmer_prg.nodes[0]->set_covg(10, 0, 1);
    pannode1.kmer_prg.nodes[2]->set_covg(10, 0, 1);
    pannode1.kmer_prg.nodes[6]->set_covg(10, 0, 1);
    pannode1.kmer_prg.nodes[9]->set_covg(10, 0, 1);

    //sample3 with top path
    sample_name = "sample3";
    sample_kmer_path = {kg1.nodes[0], kg1.nodes[1], kg1.nodes[5], kg1.nodes[9]};
    pangraph.add_node(prg_id1, prg_name1, sample_name, 2, prg_ptr1, sample_kmer_path);
    sample_kmer_path = {kg2.nodes[0], kg2.nodes[4], kg2.nodes[9]};
    pangraph.add_node(prg_id2, prg_name2, sample_name, 2, prg_ptr2, sample_kmer_path);
    pannode1.kmer_prg.nodes[0]->set_covg(2, 0, 2);
    pannode1.kmer_prg.nodes[1]->set_covg(2, 0, 2);
    pannode1.kmer_prg.nodes[5]->set_covg(2, 0, 2);
    pannode1.kmer_prg.nodes[9]->set_covg(2, 0, 2);
    pannode2.kmer_prg.nodes[0]->set_covg(2, 0, 2);
    pannode2.kmer_prg.nodes[4]->set_covg(2, 0, 2);
    pannode2.kmer_prg.nodes[9]->set_covg(2, 0, 2);

    //sample4 with bottom path
    sample_name = "sample4";
    sample_kmer_path = {kg1.nodes[0], kg1.nodes[4], kg1.nodes[9]};
    pangraph.add_node(prg_id1, prg_name1, sample_name, 3, prg_ptr1, sample_kmer_path);
    sample_kmer_path = {kg2.nodes[0], kg2.nodes[3], kg2.nodes[7], kg2.nodes[8], kg2.nodes[9]};
    pangraph.add_node(prg_id2, prg_name2, sample_name, 3, prg_ptr2, sample_kmer_path);
    pannode1.kmer_prg.nodes[0]->set_covg(5, 0, 3);
    pannode1.kmer_prg.nodes[4]->set_covg(5, 0, 3);
    pannode1.kmer_prg.nodes[9]->set_covg(5, 0, 3);
    pannode2.kmer_prg.nodes[0]->set_covg(5, 0, 3);
    pannode2.kmer_prg.nodes[3]->set_covg(5, 0, 3);
    pannode2.kmer_prg.nodes[7]->set_covg(5, 0, 3);
    pannode2.kmer_prg.nodes[8]->set_covg(5, 0, 3);
    pannode2.kmer_prg.nodes[9]->set_covg(5, 0, 3);

    VCF master_vcf;
    std::vector<LocalNodePtr> vcf_reference_path1 = {local_prg1.prg.nodes[0], local_prg1.prg.nodes[1], local_prg1.prg.nodes[3], local_prg1.prg.nodes[5], local_prg1.prg.nodes[7]};
    std::vector<LocalNodePtr> vcf_reference_path2 = {local_prg2.prg.nodes[0], local_prg2.prg.nodes[1], local_prg2.prg.nodes[3], local_prg2.prg.nodes[5], local_prg2.prg.nodes[7]};

    pannode1.construct_multisample_vcf(master_vcf, vcf_reference_path1, prg_ptr1, w, min_kmer_covg);
    pannode2.construct_multisample_vcf(master_vcf, vcf_reference_path2, prg_ptr2, w, min_kmer_covg);

    EXPECT_EQ((uint)4, master_vcf.records.size());
    EXPECT_EQ((uint)4, master_vcf.samples.size());

    //NB samples order changes to get index of each sample so can compare
    //samples 1 and 2 are ref, sample 3 is top path and sample 4 is bottom path
    auto sample1_index = master_vcf.get_sample_index("sample1");
    auto sample2_index = master_vcf.get_sample_index("sample2");
    auto sample3_index = master_vcf.get_sample_index("sample3");
    auto sample4_index = master_vcf.get_sample_index("sample4");
    std::vector<uint16_t> covgs_40 = {4,0};
    std::vector<uint16_t> covgs_100 = {10,0};
    std::vector<uint16_t> covgs_02 = {0,2};
    std::vector<uint16_t> covgs_05 = {0,5};
    std::vector<uint16_t> covgs_00 = {0,0};
    std::vector<uint16_t> covgs_400 = {4,0,0};
    std::vector<uint16_t> covgs_040 = {0,4,0};
    std::vector<uint16_t> covgs_1000 = {10,0,0};
    std::vector<uint16_t> covgs_020 = {0,2,0};
    std::vector<uint16_t> covgs_005 = {0,0,5};
    std::vector<uint16_t> covgs_000 = {0,0,0};

    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample4_index]["MEAN_FWD_COVG"], covgs_05);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample2_index]["MEAN_FWD_COVG"], covgs_100);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample1_index]["MEAN_FWD_COVG"], covgs_40);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample4_index]["MEAN_REV_COVG"], covgs_00);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample2_index]["MEAN_REV_COVG"], covgs_00);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[0].samples[sample1_index]["MEAN_REV_COVG"], covgs_00);

    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample3_index]["MEAN_FWD_COVG"], covgs_020);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample2_index]["MEAN_FWD_COVG"], covgs_1000);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample1_index]["MEAN_FWD_COVG"], covgs_400);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample4_index]["MEAN_REV_COVG"], covgs_000);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample2_index]["MEAN_REV_COVG"], covgs_000);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[1].samples[sample1_index]["MEAN_REV_COVG"], covgs_000);

    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[2].samples[sample3_index]["MEAN_FWD_COVG"], covgs_02);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[2].samples[sample3_index]["MEAN_REV_COVG"], covgs_00);

    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[3].samples[sample1_index]["MEAN_FWD_COVG"], covgs_040);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[3].samples[sample4_index]["MEAN_FWD_COVG"], covgs_005);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[3].samples[sample1_index]["MEAN_REV_COVG"], covgs_000);
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.records[3].samples[sample4_index]["MEAN_REV_COVG"], covgs_000);
}

TEST(PangenomeNodeTest, equals) {
    pangenome::Node pn1(3, 3, "3");
    pangenome::Node pn2(2, 2, "2");
    pangenome::Node pn3(2, 2, "2");

    EXPECT_EQ(pn1, pn1);
    EXPECT_EQ(pn2, pn2);
    EXPECT_EQ(pn3, pn3);
    EXPECT_EQ(pn2, pn3);
    EXPECT_EQ(pn3, pn2);
    EXPECT_EQ((pn1 == pn2), false);
    EXPECT_EQ((pn1 == pn3), false);
}

TEST(PangenomeNodeTest, nequals) {
    pangenome::Node pn1(3, 3, "3");
    pangenome::Node pn2(2, 2, "2");
    pangenome::Node pn3(2, 2, "2");

    EXPECT_EQ((pn1 != pn2), true);
    EXPECT_EQ((pn2 != pn1), true);
    EXPECT_EQ((pn1 != pn1), false);
    EXPECT_EQ((pn2 != pn2), false);
    EXPECT_EQ((pn3 != pn3), false);
    EXPECT_EQ((pn2 != pn3), false);
}

TEST(PangenomeNodeTest, less) {
    pangenome::Node pn1(3, 3, "3");
    pangenome::Node pn2(2, 2, "2");
    pangenome::Node pn3(2, 2, "2");

    EXPECT_EQ((pn1 < pn1), false);
    EXPECT_EQ((pn2 < pn2), false);
    EXPECT_EQ((pn3 < pn3), false);
    EXPECT_EQ((pn1 < pn3), false);
    EXPECT_EQ((pn1 < pn2), false);
    EXPECT_EQ((pn2 < pn1), true);
    EXPECT_EQ((pn3 < pn1), true);

}

TEST(ExtractReadsTest, get_read_overlap_coordinates) {
    //
    //  Read 0 has prg 3 sequence in interval (2,12] only
    //  Read 1 has prg 3 sequence in interval (6,16] as well as noise
    //  Read 2 has prg 3 sequence in interval (4,20] stretched out
    //  Read 3 has prg 3 sequence in interval (4,14] but is missing bits
    //  Read 4 doesn't have prg 3 sequence - on all hits are noise
    //
    uint32_t pan_node_id = 3, prg_id = 3, read_id = 0, knode_id = 0;
    string pan_node_name = "three";
    bool orientation(true);
    deque<Interval> d;
    prg::Path prg_path;
    MinimizerHitPtr mh;

    PanNodePtr pan_node = make_shared<pangenome::Node>(pan_node_id, prg_id, pan_node_name);
    PanReadPtr pr = make_shared<pangenome::Read>(read_id);

    set<MinimizerHitPtr, pComp> hits;


    // READ 0
    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(2, 5), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(3, 6), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(16, 17), Interval(27, 29) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(27, 30) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // READ 1
    read_id = 1;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(12, 15), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(11, 14), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(16, 17), Interval(27, 29) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(27, 30) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(31, 33) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(78, 81) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // READ 2
    read_id = 2;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(17, 20), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(15, 18), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(16, 17), Interval(27, 29) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(27, 30) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(31, 33) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(78, 81) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // READ 3
    read_id = 3;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(16, 17), Interval(27, 29) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);


    // noise
    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // READ 4
    read_id = 4;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(17, 20), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(31, 33) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(78, 81) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // RUN GET_READ_OVERLAPS
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    const std::vector<LocalNodePtr> lmp {//l3.prg.nodes[0],
            l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6], l3.prg.nodes[7]//, l3.prg.nodes[9]
    };
    // A G C T CGG  TAT
    const std::set<ReadCoordinate> expected_overlaps {{ 0, 3, 9,  1 },
                                                      { 1, 7, 13, 1 },
                                                      { 2, 5, 13, 1 },
                                                      { 3, 6, 10, 1 }};

    prg::Path local_path;
    for (const auto &node : lmp) {
        local_path.add_end_interval(node->pos);
    }
    const auto overlaps { pan_node->get_read_overlap_coordinates(local_path) };

    EXPECT_ITERABLE_EQ(std::set<ReadCoordinate>, expected_overlaps, overlaps);
}


TEST(ExtractReadsTest, get_read_overlap_coordinates_no_duplicates) {
    //
    //  Read 0 has prg 3 sequence in interval (2,12] only
    //  Read 1 has prg 3 sequence in interval (6,16] as well as noise
    //  Read 2 has prg 3 sequence in interval (4,20] stretched out
    //  Read 3 has prg 3 sequence in interval (4,14] but is missing bits
    //  Read 4 doesn't have prg 3 sequence - on all hits are noise
    //  Read 5 is a duplicate of Read 0
    //
    uint32_t pan_node_id = 3, prg_id = 3, read_id = 0, knode_id = 0;
    string pan_node_name = "three";
    bool orientation(true);
    deque<Interval> d;
    prg::Path prg_path;
    MinimizerHitPtr mh;

    PanNodePtr pan_node = make_shared<pangenome::Node>(pan_node_id, prg_id, pan_node_name);
    PanReadPtr pr = make_shared<pangenome::Read>(read_id);

    set<MinimizerHitPtr, pComp> hits;


    // READ 0
    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(2, 5), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(3, 6), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(16, 17), Interval(27, 29) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(27, 30) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // READ 1
    read_id = 1;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(12, 15), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(11, 14), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(16, 17), Interval(27, 29) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(27, 30) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(31, 33) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(78, 81) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // READ 2
    read_id = 2;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(17, 20), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(15, 18), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(16, 17), Interval(27, 29) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(27, 30) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(31, 33) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(78, 81) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // READ 3
    read_id = 3;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10, 13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(16, 17), Interval(27, 29) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);


    // noise
    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // READ 4
    read_id = 4;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(17, 20), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1, 4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(31, 33) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9, 12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(78, 81) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13, 16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // READ 5
    read_id = 0;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(2, 5), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8, 11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7, 10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(3, 6), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4, 7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(16, 17), Interval(27, 29) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5, 8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = { Interval(27, 30) };
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6, 9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id, hits);
    pan_node->reads.insert(pr);
    hits.clear();

    // RUN GET_READ_OVERLAPS
    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    const std::vector<LocalNodePtr> lmp {//l3.prg.nodes[0],
            l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6], l3.prg.nodes[7]//, l3.prg.nodes[9]
    };
    // A G C T CGG  TAT
    const std::set<ReadCoordinate> expected_overlaps {{ 0, 3, 9,  1 },
                                                      { 1, 7, 13, 1 },
                                                      { 2, 5, 13, 1 },
                                                      { 3, 6, 10, 1 }};

    prg::Path local_path;
    for (const auto &node : lmp) {
        local_path.add_end_interval(node->pos);
    }
    const auto overlaps { pan_node->get_read_overlap_coordinates(local_path) };

    EXPECT_ITERABLE_EQ(std::set<ReadCoordinate>, expected_overlaps, overlaps);
}