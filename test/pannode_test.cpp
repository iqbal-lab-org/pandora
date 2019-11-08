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
#include "test_helpers.h"


using namespace pangenome;

TEST(PangenomeNodeTest, create) {
    auto local_graph_ptr { std::make_shared<LocalPRG>(4, "3", "") };
    pangenome::Node pan_node(local_graph_ptr, 3);
    uint32_t j = 3;
    EXPECT_EQ(j, pan_node.node_id);
    EXPECT_EQ((uint) 4, pan_node.prg_id);
    EXPECT_EQ("3", pan_node.name);
    EXPECT_EQ((uint) 0, pan_node.covg);
    EXPECT_EQ((uint) 0, pan_node.reads.size());
    EXPECT_EQ((uint) 0, pan_node.samples.size());
}

TEST(PangenomeNodeTest, get_name) {
    auto l1 { std::make_shared<LocalPRG>(3, "3", "") };
    auto l2 { std::make_shared<LocalPRG>(2, "2", "") };
    pangenome::Node pn1(l1);
    pangenome::Node pn2(l2);
    pangenome::Node pn3(l2, 4);

    EXPECT_EQ(pn1.get_name(), "3");
    EXPECT_EQ(pn2.get_name(), "2");
    EXPECT_EQ(pn3.get_name(), "2.4");
}

TEST(PangenomeNodeTest, add_path) {
    //setup the KmerGraph
    auto local_prg_ptr { std::make_shared<LocalPRG>(3, "3", "") };
    std::deque<Interval> d = {Interval(0, 0)};
    prg::Path p;
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    EXPECT_EQ((uint) 7, local_prg_ptr->kmer_prg.nodes.size());

    //setup the Node
    pangenome::Node pn1(local_prg_ptr);
    std::vector<KmerNodePtr> kmp;
    pn1.add_path(kmp, 0);


    //do the tests
    EXPECT_EQ((uint) 7, pn1.kmer_prg_with_coverage.kmer_prg->nodes.size());
    const auto& nodes = pn1.kmer_prg_with_coverage.kmer_prg->nodes;
    kmp = {nodes[0], nodes[3], nodes[4], nodes[6]};
    pn1.add_path(kmp, 0);
    EXPECT_EQ((uint) 1, pn1.kmer_prg_with_coverage.get_covg(0, 0, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg_with_coverage.get_covg(1, 0, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg_with_coverage.get_covg(2, 0, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg_with_coverage.get_covg(3, 0, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg_with_coverage.get_covg(4, 0, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg_with_coverage.get_covg(5, 0, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg_with_coverage.get_covg(6, 0, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg_with_coverage.get_covg(0, 1, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg_with_coverage.get_covg(1, 1, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg_with_coverage.get_covg(2, 1, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg_with_coverage.get_covg(3, 1, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg_with_coverage.get_covg(4, 1, 0));
    EXPECT_EQ((uint) 0, pn1.kmer_prg_with_coverage.get_covg(5, 1, 0));
    EXPECT_EQ((uint) 1, pn1.kmer_prg_with_coverage.get_covg(6, 1, 0));
}

TEST(PangenomeNodeTest, get_read_overlap_coordinates) {
    auto local_prg_ptr { std::make_shared<LocalPRG>(3, "3", "") };
    auto pan_node_ptr = std::make_shared<pangenome::Node>(local_prg_ptr);
    pangenome::ReadPtr pr;
    MinimizerHits mhits;

    std::deque<Interval> d;
    prg::Path p;

    // read1
    Minimizer m1(0, 1, 6, 0); // kmer, start, end, strand
    d = {Interval(7, 8), Interval(10, 14)};
    p.initialize(d);
    MiniRecord mr1(3, p, 0, 0);
    mhits.add_hit(1, m1, mr1); // read 1

    Minimizer m2(0, 0, 5, 0);
    d = {Interval(6, 10), Interval(11, 12)};
    p.initialize(d);
    MiniRecord mr2(3, p, 0, 0);
    mhits.add_hit(1, m2, mr2);

    Minimizer m3(0, 0, 5, 0);
    d = {Interval(6, 10), Interval(12, 13)};
    p.initialize(d);
    MiniRecord mr3(3, p, 0, 0);
    mhits.add_hit(1, m3, mr3);

    pr = std::make_shared<pangenome::Read>(1);
    pr->add_hits(pan_node_ptr, mhits.hits);
    pan_node_ptr->reads.insert(pr);
    mhits.clear();

    //read 2
    Minimizer m4(0, 2, 7, 1);
    d = {Interval(6, 10), Interval(11, 12)};
    p.initialize(d);
    MiniRecord mr4(3, p, 0, 0);
    mhits.add_hit(2, m4, mr4);

    Minimizer m5(0, 5, 10, 1);
    d = {Interval(6, 10), Interval(12, 13)};
    p.initialize(d);
    MiniRecord mr5(3, p, 0, 0);
    mhits.add_hit(2, m5, mr5);

    pr = std::make_shared<pangenome::Read>(2);
    pr->add_hits(pan_node_ptr, mhits.hits);
    pan_node_ptr->reads.insert(pr);
    mhits.clear();

    std::vector<std::vector<uint32_t>> read_overlap_coordinates;
    pan_node_ptr->get_read_overlap_coordinates(read_overlap_coordinates);
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

class PangenomeNodeTest___construct_multisample_vcf___Fixture : public ::testing::Test {
protected:
    PangenomeNodeTest___construct_multisample_vcf___Fixture() :
            w(1), k(3), min_kmer_covg(0), index(std::make_shared<Index>()), sample_names({"sample1", "sample2", "sample3", "sample4"}),
            pangraph(sample_names), master_vcf(create_VCF_with_default_parameters()),
            nested_varsite_PRG(std::make_shared<LocalPRG>(0, "nested varsite" , "A 5 G 7 C 8 T 8 CT 7  6 G 5 T")),
            modified_PRG(std::make_shared<LocalPRG>(1, "modified" , "A 5 G 7 G 8 A 8 GA 7  6 G 5 T"))
    {

    }

    void SetUp() override {
        nested_varsite_kmer_coverage_graph_pointer = &nested_varsite_PRG->kmer_prg;
        nested_varsite_vcf_reference_path = {nested_varsite_PRG->prg.nodes[0], nested_varsite_PRG->prg.nodes[1], nested_varsite_PRG->prg.nodes[3], nested_varsite_PRG->prg.nodes[5], nested_varsite_PRG->prg.nodes[7]};

        modified_kmer_coverage_graph_pointer = &modified_PRG->kmer_prg;
        modified_vcf_reference_path = {modified_PRG->prg.nodes[0], modified_PRG->prg.nodes[1], modified_PRG->prg.nodes[3], modified_PRG->prg.nodes[5], modified_PRG->prg.nodes[7]};
    }

    void TearDown() override {
    }

    uint32_t w, k, min_kmer_covg;
    std::shared_ptr<Index> index;
    std::vector<std::string> sample_names;
    pangenome::Graph pangraph;
    VCF master_vcf;

    std::shared_ptr<LocalPRG> nested_varsite_PRG;
    KmerGraph* nested_varsite_kmer_coverage_graph_pointer;
    std::vector<LocalNodePtr> nested_varsite_vcf_reference_path;

    std::shared_ptr<LocalPRG> modified_PRG;
    KmerGraph* modified_kmer_coverage_graph_pointer;
    std::vector<LocalNodePtr> modified_vcf_reference_path;


    void test_format_fields_for_each_sample_and_allele_for_a_given_record(const VCFRecord &record, uint32_t number_of_samples, uint32_t number_of_alleles) {
        for (uint32_t sample_index = 0; sample_index < number_of_samples; ++sample_index) {
            for (uint32_t allele_index = 0; allele_index < number_of_alleles; ++allele_index) {
                EXPECT_EQ(0, record.sampleIndex_to_sampleInfo[sample_index].get_mean_forward_coverage(allele_index));
                EXPECT_EQ(0, record.sampleIndex_to_sampleInfo[sample_index].get_mean_reverse_coverage(allele_index));

                EXPECT_EQ(0, record.sampleIndex_to_sampleInfo[sample_index].get_median_forward_coverage(allele_index));
                EXPECT_EQ(0, record.sampleIndex_to_sampleInfo[sample_index].get_median_reverse_coverage(allele_index));

                EXPECT_EQ(0, record.sampleIndex_to_sampleInfo[sample_index].get_sum_forward_coverage(allele_index));
                EXPECT_EQ(0, record.sampleIndex_to_sampleInfo[sample_index].get_sum_reverse_coverage(allele_index));

                EXPECT_EQ(0, record.sampleIndex_to_sampleInfo[sample_index].get_gaps(allele_index));
            }
        }
    };

};

TEST_F(PangenomeNodeTest___construct_multisample_vcf___Fixture, construct_multisample_vcf_single_prg)
{
    nested_varsite_PRG->minimizer_sketch(index, w, k);

    //sample1
    std::vector<KmerNodePtr> sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[2], nested_varsite_kmer_coverage_graph_pointer->nodes[6], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_names[0], sample_kmer_path);

    //sample2 identical to sample1
    sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[2], nested_varsite_kmer_coverage_graph_pointer->nodes[6], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_names[1], sample_kmer_path);

    //sample3 with top path
    sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[1], nested_varsite_kmer_coverage_graph_pointer->nodes[5], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_names[2], sample_kmer_path);

    //sample4 with bottom path
    sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[4], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_names[3], sample_kmer_path);

    auto &pannode = *pangraph.nodes[nested_varsite_PRG->id];
    pannode.construct_multisample_vcf(master_vcf, nested_varsite_vcf_reference_path, nested_varsite_PRG, w);



    EXPECT_EQ((uint)2, master_vcf.get_VCF_size());
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
    uint16_t alt_gt = 1;
    uint16_t ref_gt = 0;

    EXPECT_EQ((uint)1, master_vcf.get_records()[0]->pos);
    EXPECT_EQ("GT", master_vcf.get_records()[0]->ref);
    EXPECT_EQ((uint)1, master_vcf.get_records()[0]->alts.size());
    EXPECT_EQ("G", master_vcf.get_records()[0]->alts[0]);
    EXPECT_EQ((uint)4, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_TRUE(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample4_index].is_gt_from_max_likelihood_path_valid()) ;
    EXPECT_FALSE(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample3_index].is_gt_from_max_likelihood_path_valid()) ;
    EXPECT_TRUE(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample2_index].is_gt_from_max_likelihood_path_valid()) ;
    EXPECT_TRUE(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample1_index].is_gt_from_max_likelihood_path_valid()) ;
    EXPECT_EQ(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample4_index].get_gt_from_max_likelihood_path(), alt_gt);
    EXPECT_EQ(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample2_index].get_gt_from_max_likelihood_path(), ref_gt);
    EXPECT_EQ(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample1_index].get_gt_from_max_likelihood_path(), ref_gt);

    test_format_fields_for_each_sample_and_allele_for_a_given_record(*(master_vcf.get_records()[0]), sample_names.size(), 2);

    EXPECT_EQ((uint)2, master_vcf.get_records()[1]->pos);
    EXPECT_EQ("T", master_vcf.get_records()[1]->ref);
    EXPECT_EQ((uint)2, master_vcf.get_records()[1]->alts.size());
    EXPECT_EQ("C", master_vcf.get_records()[1]->alts[0]);
    EXPECT_EQ("CT", master_vcf.get_records()[1]->alts[1]);
    EXPECT_EQ((uint)4, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample4_index].is_gt_from_max_likelihood_path_valid());
    EXPECT_TRUE(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample3_index].is_gt_from_max_likelihood_path_valid());
    EXPECT_TRUE(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample2_index].is_gt_from_max_likelihood_path_valid());
    EXPECT_TRUE(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample1_index].is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample3_index].get_gt_from_max_likelihood_path(), alt_gt);
    EXPECT_EQ(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample2_index].get_gt_from_max_likelihood_path(), ref_gt);
    EXPECT_EQ(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample1_index].get_gt_from_max_likelihood_path(), ref_gt);

    test_format_fields_for_each_sample_and_allele_for_a_given_record(*(master_vcf.get_records()[1]), sample_names.size(), 3);
}

//TEST_F(PangenomeNodeTest___construct_multisample_vcf___Fixture,construct_multisample_vcf_two_prg)
//{
//    nested_varsite_PRG->minimizer_sketch(index, w, k);
//    modified_PRG->minimizer_sketch(index, w, k);
//
//    //sample1
//    std::vector<KmerNodePtr> sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[2], nested_varsite_kmer_coverage_graph_pointer->nodes[6], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(nested_varsite_PRG);
//    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_names[0], sample_kmer_path);
//    sample_kmer_path = {modified_kmer_coverage_graph_pointer->nodes[0], modified_kmer_coverage_graph_pointer->nodes[1], modified_kmer_coverage_graph_pointer->nodes[5], modified_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(modified_PRG);
//    pangraph.add_hits_between_PRG_and_sample(modified_PRG->id, sample_names[0], sample_kmer_path);
//
//    //sample2 identical to sample1 in prg1, no prg2
//    sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[2], nested_varsite_kmer_coverage_graph_pointer->nodes[6], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(nested_varsite_PRG);
//    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_names[1], sample_kmer_path);
//
//    //sample3 with top path
//    sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[1], nested_varsite_kmer_coverage_graph_pointer->nodes[5], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(nested_varsite_PRG);
//    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_names[2], sample_kmer_path);
//    sample_kmer_path = {modified_kmer_coverage_graph_pointer->nodes[0], modified_kmer_coverage_graph_pointer->nodes[4], modified_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(modified_PRG);
//    pangraph.add_hits_between_PRG_and_sample(modified_PRG->id, sample_names[2], sample_kmer_path);
//
//
//    //sample4 with bottom path
//    sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[4], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(nested_varsite_PRG);
//    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_names[3], sample_kmer_path);
//    sample_kmer_path = {modified_kmer_coverage_graph_pointer->nodes[0], modified_kmer_coverage_graph_pointer->nodes[3], modified_kmer_coverage_graph_pointer->nodes[7], modified_kmer_coverage_graph_pointer->nodes[8], modified_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(modified_PRG);
//    pangraph.add_hits_between_PRG_and_sample(modified_PRG->id, sample_names[3], sample_kmer_path);
//
//
//    auto &pannode1 = *pangraph.nodes[nested_varsite_PRG->id];
//    pannode1.construct_multisample_vcf(master_vcf, nested_varsite_vcf_reference_path, nested_varsite_PRG, w, default_genotyping_options);
//    auto &pannode2 = *pangraph.nodes[modified_PRG->id];
//    pannode2.construct_multisample_vcf(master_vcf, modified_vcf_reference_path, modified_PRG, w, default_genotyping_options);
//
//    EXPECT_EQ((uint)4, master_vcf.get_VCF_size());
//    EXPECT_EQ((uint)4, master_vcf.samples.size());
//
//    //NB samples order changes to get index of each sample so can compare
//    //samples 1 and 2 are ref, sample 3 is top path and sample 4 is bottom path
//    auto iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample1");
//    auto sample1_index = std::distance(master_vcf.samples.begin(), iter);
//    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample2");
//    auto sample2_index = std::distance(master_vcf.samples.begin(), iter);
//    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample3");
//    auto sample3_index = std::distance(master_vcf.samples.begin(), iter);
//    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample4");
//    auto sample4_index = std::distance(master_vcf.samples.begin(), iter);
//    uint16_t alt_gt = 1;
//    uint16_t ref_gt = 0;
//    uint16_t alt2_gt = 2;
//
//    EXPECT_EQ((uint)1, master_vcf.get_records()[0]->pos);
//    EXPECT_EQ("GT", master_vcf.get_records()[0]->ref);
//    EXPECT_EQ((uint)1, master_vcf.get_records()[0]->alts.size());
//    EXPECT_EQ("G", master_vcf.get_records()[0]->alts[0]);
//    EXPECT_EQ((uint)4, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
//    EXPECT_TRUE(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample4_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_FALSE(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample3_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_TRUE(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample2_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_TRUE(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample1_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_EQ(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample4_index].get_gt_from_max_likelihood_path(), alt_gt);
//    EXPECT_EQ(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample2_index].get_gt_from_max_likelihood_path(), ref_gt);
//    EXPECT_EQ(master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample1_index].get_gt_from_max_likelihood_path(), ref_gt);
//    test_format_fields_for_each_sample_and_allele_for_a_given_record(*(master_vcf.get_records()[0]), 4, 2);
//
//
//    EXPECT_EQ((uint)2, master_vcf.get_records()[1]->pos);
//    EXPECT_EQ("T", master_vcf.get_records()[1]->ref);
//    EXPECT_EQ((uint)2, master_vcf.get_records()[1]->alts.size());
//    EXPECT_EQ("C", master_vcf.get_records()[1]->alts[0]);
//    EXPECT_EQ("CT", master_vcf.get_records()[1]->alts[1]);
//    EXPECT_EQ((uint)4, master_vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());
//    EXPECT_FALSE(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample4_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_TRUE(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample3_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_TRUE(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample2_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_TRUE(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample1_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_EQ(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample3_index].get_gt_from_max_likelihood_path(), alt_gt);
//    EXPECT_EQ(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample2_index].get_gt_from_max_likelihood_path(), ref_gt);
//    EXPECT_EQ(master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample1_index].get_gt_from_max_likelihood_path(), ref_gt);
//    test_format_fields_for_each_sample_and_allele_for_a_given_record(*(master_vcf.get_records()[1]), 4, 3);
//
//
//    EXPECT_EQ((uint)1, master_vcf.get_records()[2]->pos);
//    EXPECT_EQ("GA", master_vcf.get_records()[2]->ref);
//    EXPECT_EQ((uint)1, master_vcf.get_records()[2]->alts.size());
//    EXPECT_EQ("G", master_vcf.get_records()[2]->alts[0]);
//    EXPECT_EQ((uint)4, master_vcf.get_records()[2]->sampleIndex_to_sampleInfo.size());
//    EXPECT_FALSE(master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample4_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_TRUE(master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample3_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_FALSE(master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample2_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_FALSE(master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample1_index].is_gt_from_max_likelihood_path_valid()) ;
//    EXPECT_EQ(master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample3_index].get_gt_from_max_likelihood_path(), alt_gt);
//    test_format_fields_for_each_sample_and_allele_for_a_given_record(*(master_vcf.get_records()[2]), 4, 2);

//    EXPECT_EQ((uint)2, master_vcf.get_records()[3]->pos);
//    EXPECT_EQ("A", master_vcf.get_records()[3]->ref);
//    EXPECT_EQ((uint)2, master_vcf.get_records()[3]->alts.size());
//    EXPECT_EQ("G", master_vcf.get_records()[3]->alts[0]);
//    EXPECT_EQ("GA", master_vcf.get_records()[3]->alts[1]);
//    EXPECT_EQ((uint)4, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo.size());
//    EXPECT_FALSE(master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index].find("GT") == master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index].end()) ;
//    EXPECT_TRUE(master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample3_index].find("GT") == master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample3_index].end()) ;
//    EXPECT_TRUE(master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample2_index].find("GT") == master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample2_index].end()) ;
//    EXPECT_FALSE(master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index].find("GT") == master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index].end()) ;
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index]["GT"], alt_gt);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index]["GT"], alt2_gt);
//
//    for (const auto format : formats){
//        EXPECT_FALSE(master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index].find(format) == master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index].end());
//        EXPECT_TRUE(master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample2_index].find(format) == master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample2_index].end());
//        EXPECT_FALSE(master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample3_index].find(format) == master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample3_index].end());
//        EXPECT_FALSE(master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index].find(format) == master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index].end());
//        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index][format], no_covg3);
//        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample3_index][format], no_covg3);
//        EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index][format], no_covg3);
//    }
//}
//
//TEST(PangenomeNodeTest,construct_multisample_vcf_two_prg_with_covgs)
//{
//    uint32_t prg_id1 = 3, w=1, k=3, min_kmer_covg=0;
//    std::string prg_name1 = "nested varsite";
//    LocalPRG local_prg1(nested_varsite_PRG->id, prg_name1 , "A 5 G 7 C 8 T 8 CT 7  6 G 5 T");
//    uint32_t modified_PRG->id = 5;
//    std::string prg_name2 = "modified";
//    LocalPRG local_prg2(modified_PRG->id, prg_name2 , "A 5 G 7 G 8 A 8 GA 7  6 G 5 T");
//
//    auto index = std::make_shared<Index>();
//    local_prg1.minimizer_sketch(index, w, k);
//    auto nested_varsite_PRG = std::make_shared<LocalPRG>(local_prg1);
//    auto &kg1 = local_prg1.kmer_prg;
//    local_prg2.minimizer_sketch(index, w, k);
//    auto modified_PRG = std::make_shared<LocalPRG>(local_prg2);
//    auto &kg2 = local_prg2.kmer_prg;
//
//    pangenome::Graph pangraph({"sample1", "sample2", "sample3", "sample4"});
//
//    //sample1
//    std::string sample_name = "sample1";
//    std::vector<KmerNodePtr> sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[2], nested_varsite_kmer_coverage_graph_pointer->nodes[6], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(nested_varsite_PRG);
//    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_name, sample_kmer_path);
//    sample_kmer_path = {modified_kmer_coverage_graph_pointer->nodes[0], modified_kmer_coverage_graph_pointer->nodes[1], modified_kmer_coverage_graph_pointer->nodes[5], modified_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(modified_PRG);
//    pangraph.add_hits_between_PRG_and_sample(modified_PRG->id, sample_name, sample_kmer_path);
//
//    auto &pannode1 = *pangraph.nodes[nested_varsite_PRG->id];
//    auto &pannode2 = *pangraph.nodes[modified_PRG->id];
//
//    pannode1.kmer_prg_with_coverage.set_covg(0, 4, 0, 0);
//    pannode1.kmer_prg_with_coverage.set_covg(2, 4, 0, 0);
//    pannode1.kmer_prg_with_coverage.set_covg(6, 4, 0, 0);
//    pannode1.kmer_prg_with_coverage.set_covg(9, 4, 0, 0);
//    pannode2.kmer_prg_with_coverage.set_covg(0, 4, 0, 0);
//    pannode2.kmer_prg_with_coverage.set_covg(1, 4, 0, 0);
//    pannode2.kmer_prg_with_coverage.set_covg(5, 4, 0, 0);
//    pannode2.kmer_prg_with_coverage.set_covg(9, 4, 0, 0);
//
//    //sample2 identical to sample1 in prg1, no prg2
//    sample_name = "sample2";
//    sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[2], nested_varsite_kmer_coverage_graph_pointer->nodes[6], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(nested_varsite_PRG);
//    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_name, sample_kmer_path);
//    pannode1.kmer_prg_with_coverage.set_covg(0, 10, 0, 1);
//    pannode1.kmer_prg_with_coverage.set_covg(2, 10, 0, 1);
//    pannode1.kmer_prg_with_coverage.set_covg(6, 10, 0, 1);
//    pannode1.kmer_prg_with_coverage.set_covg(9, 10, 0, 1);
//
//    //sample3 with top path
//    sample_name = "sample3";
//    sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[1], nested_varsite_kmer_coverage_graph_pointer->nodes[5], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(nested_varsite_PRG);
//    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_name, sample_kmer_path);
//    sample_kmer_path = {modified_kmer_coverage_graph_pointer->nodes[0], modified_kmer_coverage_graph_pointer->nodes[4], modified_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(modified_PRG);
//    pangraph.add_hits_between_PRG_and_sample(modified_PRG->id, sample_name, sample_kmer_path);
//    pannode1.kmer_prg_with_coverage.set_covg(0, 2, 0, 2);
//    pannode1.kmer_prg_with_coverage.set_covg(1, 2, 0, 2);
//    pannode1.kmer_prg_with_coverage.set_covg(5, 2, 0, 2);
//    pannode1.kmer_prg_with_coverage.set_covg(9, 2, 0, 2);
//    pannode2.kmer_prg_with_coverage.set_covg(0, 2, 0, 2);
//    pannode2.kmer_prg_with_coverage.set_covg(4, 2, 0, 2);
//    pannode2.kmer_prg_with_coverage.set_covg(9, 2, 0, 2);
//
//    //sample4 with bottom path
//    sample_name = "sample4";
//    sample_kmer_path = {nested_varsite_kmer_coverage_graph_pointer->nodes[0], nested_varsite_kmer_coverage_graph_pointer->nodes[4], nested_varsite_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(nested_varsite_PRG);
//    pangraph.add_hits_between_PRG_and_sample(nested_varsite_PRG->id, sample_name, sample_kmer_path);
//    sample_kmer_path = {modified_kmer_coverage_graph_pointer->nodes[0], modified_kmer_coverage_graph_pointer->nodes[3], modified_kmer_coverage_graph_pointer->nodes[7], modified_kmer_coverage_graph_pointer->nodes[8], modified_kmer_coverage_graph_pointer->nodes[9]};
//    pangraph.add_node(modified_PRG);
//    pangraph.add_hits_between_PRG_and_sample(modified_PRG->id, sample_name, sample_kmer_path);
//    pannode1.kmer_prg_with_coverage.set_covg(0, 5, 0, 3);
//    pannode1.kmer_prg_with_coverage.set_covg(4, 5, 0, 3);
//    pannode1.kmer_prg_with_coverage.set_covg(9, 5, 0, 3);
//    pannode2.kmer_prg_with_coverage.set_covg(0, 5, 0, 3);
//    pannode2.kmer_prg_with_coverage.set_covg(3, 5, 0, 3);
//    pannode2.kmer_prg_with_coverage.set_covg(7, 5, 0, 3);
//    pannode2.kmer_prg_with_coverage.set_covg(8, 5, 0, 3);
//    pannode2.kmer_prg_with_coverage.set_covg(9, 5, 0, 3);
//
//    VCF master_vcf;
//    std::vector<LocalNodePtr> vcf_reference_path1 = {local_prg1.prg.nodes[0], local_prg1.prg.nodes[1], local_prg1.prg.nodes[3], local_prg1.prg.nodes[5], local_prg1.prg.nodes[7]};
//    std::vector<LocalNodePtr> vcf_reference_path2 = {local_prg2.prg.nodes[0], local_prg2.prg.nodes[1], local_prg2.prg.nodes[3], local_prg2.prg.nodes[5], local_prg2.prg.nodes[7]};
//
//    pannode1.construct_multisample_vcf(master_vcf, vcf_reference_path1, nested_varsite_PRG, w, min_kmer_covg);
//    pannode2.construct_multisample_vcf(master_vcf, vcf_reference_path2, modified_PRG, w, min_kmer_covg);
//
//    EXPECT_EQ((uint)4, master_vcf.get_VCF_size());
//    EXPECT_EQ((uint)4, master_vcf.samples.size());
//
//    //NB samples order changes to get index of each sample so can compare
//    //samples 1 and 2 are ref, sample 3 is top path and sample 4 is bottom path
//    auto sample1_index = master_vcf.get_sample_index("sample1");
//    auto sample2_index = master_vcf.get_sample_index("sample2");
//    auto sample3_index = master_vcf.get_sample_index("sample3");
//    auto sample4_index = master_vcf.get_sample_index("sample4");
//    std::vector<uint16_t> covgs_40 = {4,0};
//    std::vector<uint16_t> covgs_100 = {10,0};
//    std::vector<uint16_t> covgs_02 = {0,2};
//    std::vector<uint16_t> covgs_05 = {0,5};
//    std::vector<uint16_t> covgs_00 = {0,0};
//    std::vector<uint16_t> covgs_400 = {4,0,0};
//    std::vector<uint16_t> covgs_040 = {0,4,0};
//    std::vector<uint16_t> covgs_1000 = {10,0,0};
//    std::vector<uint16_t> covgs_020 = {0,2,0};
//    std::vector<uint16_t> covgs_005 = {0,0,5};
//    std::vector<uint16_t> covgs_000 = {0,0,0};
//
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample4_index]["MEAN_FWD_COVG"], covgs_05);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample2_index]["MEAN_FWD_COVG"], covgs_100);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_FWD_COVG"], covgs_40);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample4_index]["MEAN_REV_COVG"], covgs_00);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample2_index]["MEAN_REV_COVG"], covgs_00);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_REV_COVG"], covgs_00);
//
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample3_index]["MEAN_FWD_COVG"], covgs_020);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample2_index]["MEAN_FWD_COVG"], covgs_1000);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_FWD_COVG"], covgs_400);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample4_index]["MEAN_REV_COVG"], covgs_000);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample2_index]["MEAN_REV_COVG"], covgs_000);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_REV_COVG"], covgs_000);
//
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample3_index]["MEAN_FWD_COVG"], covgs_02);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample3_index]["MEAN_REV_COVG"], covgs_00);
//
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_FWD_COVG"], covgs_040);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index]["MEAN_FWD_COVG"], covgs_005);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_REV_COVG"], covgs_000);
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index]["MEAN_REV_COVG"], covgs_000);
//}
//
//TEST(PangenomeNodeTest, equals) {
//    auto l1 { std::make_shared<LocalPRG>(3, "3", "") };
//    auto l2 { std::make_shared<LocalPRG>(2, "2", "") };
//    pangenome::Node pn1(l1);
//    pangenome::Node pn2(l2);
//    pangenome::Node pn3(l2);
//
//    EXPECT_EQ(pn1, pn1);
//    EXPECT_EQ(pn2, pn2);
//    EXPECT_EQ(pn3, pn3);
//    EXPECT_EQ(pn2, pn3);
//    EXPECT_EQ(pn3, pn2);
//    EXPECT_EQ((pn1 == pn2), false);
//    EXPECT_EQ((pn1 == pn3), false);
//}
//
//TEST(PangenomeNodeTest, nequals) {
//    auto l1 { std::make_shared<LocalPRG>(3, "3", "") };
//    auto l2 { std::make_shared<LocalPRG>(2, "2", "") };
//    pangenome::Node pn1(l1);
//    pangenome::Node pn2(l2);
//    pangenome::Node pn3(l2);
//
//    EXPECT_EQ((pn1 != pn2), true);
//    EXPECT_EQ((pn2 != pn1), true);
//    EXPECT_EQ((pn1 != pn1), false);
//    EXPECT_EQ((pn2 != pn2), false);
//    EXPECT_EQ((pn3 != pn3), false);
//    EXPECT_EQ((pn2 != pn3), false);
//}
//
//TEST(PangenomeNodeTest, less) {
//    auto l1 { std::make_shared<LocalPRG>(3, "3", "") };
//    auto l2 { std::make_shared<LocalPRG>(2, "2", "") };
//    pangenome::Node pn1(l1);
//    pangenome::Node pn2(l2);
//    pangenome::Node pn3(l2);
//
//    EXPECT_EQ((pn1 < pn1), false);
//    EXPECT_EQ((pn2 < pn2), false);
//    EXPECT_EQ((pn3 < pn3), false);
//    EXPECT_EQ((pn1 < pn3), false);
//    EXPECT_EQ((pn1 < pn2), false);
//    EXPECT_EQ((pn2 < pn1), true);
//    EXPECT_EQ((pn3 < pn1), true);
//
//}
//
//TEST(ExtractReadsTest, get_read_overlap_coordinates) {
//    //
//    //  Read 0 has prg 3 sequence in interval (2,12] only
//    //  Read 1 has prg 3 sequence in interval (6,16] as well as noise
//    //  Read 2 has prg 3 sequence in interval (4,20] stretched out
//    //  Read 3 has prg 3 sequence in interval (4,14] but is missing bits
//    //  Read 4 doesn't have prg 3 sequence - on all hits are noise
//    //
//    uint32_t read_id = 0, knode_id = 0;
//
//    bool orientation(true);
//    deque<Interval> d;
//    prg::Path prg_path;
//    MinimizerHitPtr mh;
//
//    uint32_t prg_id = 3;
//    auto local_prg_ptr { std::make_shared<LocalPRG>(prg_id, "three", "") };
//    PanNodePtr pan_node = make_shared<pangenome::Node>(local_prg_ptr);
//    PanReadPtr pr = make_shared<pangenome::Read>(read_id);
//    set<MinimizerHitPtr, pComp> hits;
//
//
//    // READ 0
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m1(0, 2, 5, orientation); // kmer, start, end, strand
//    MiniRecord mr1(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m1, mr1);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m2(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr2(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m2, mr2);
//    hits.insert(mh);
//    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
//    prg_path.initialize(d);
//    Minimizer m3(0, 7, 10, orientation); // kmer, start, end, strand
//    MiniRecord mr3(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m3, mr3);
//    hits.insert(mh);
//
//    // hits on path
//    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
//    prg_path.initialize(d);
//    Minimizer m4(0, 3, 6, orientation); // kmer, start, end, strand
//    MiniRecord mr4(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m4, mr4);
//    hits.insert(mh);
//    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m5(0, 4, 7, orientation); // kmer, start, end, strand
//    MiniRecord mr5(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m5, mr5);
//    hits.insert(mh);
//    d = { Interval(16, 17), Interval(27, 29) };
//    prg_path.initialize(d);
//    Minimizer m6(0, 5, 8, orientation); // kmer, start, end, strand
//    MiniRecord mr6(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m6, mr6);
//    hits.insert(mh);
//    d = { Interval(27, 30) };
//    prg_path.initialize(d);
//    Minimizer m7(0, 6, 9, orientation); // kmer, start, end, strand
//    MiniRecord mr7(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m7, mr7);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // READ 1
//    read_id = 1;
//    pr = make_shared<pangenome::Read>(read_id);
//
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m8(0, 6, 9, orientation); // kmer, start, end, strand
//    MiniRecord mr8(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m8, mr8);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m9(0, 12, 15, orientation); // kmer, start, end, strand
//    MiniRecord mr9(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m9, mr9);
//    hits.insert(mh);
//    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
//    prg_path.initialize(d);
//    Minimizer m10(0, 11, 14, orientation); // kmer, start, end, strand
//    MiniRecord mr10(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m10, mr10);
//    hits.insert(mh);
//
//    // hits on path
//    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
//    prg_path.initialize(d);
//    Minimizer m11(0, 7, 10, orientation); // kmer, start, end, strand
//    MiniRecord mr11(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m11, mr11);
//    hits.insert(mh);
//    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m12(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr12(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m12, mr12);
//    hits.insert(mh);
//    d = { Interval(16, 17), Interval(27, 29) };
//    prg_path.initialize(d);
//    Minimizer m13(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr13(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m13, mr13);
//    hits.insert(mh);
//    d = { Interval(27, 30) };
//    prg_path.initialize(d);
//    Minimizer m14(0, 10, 13, orientation); // kmer, start, end, strand
//    MiniRecord mr14(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m14, mr14);
//    hits.insert(mh);
//
//    // noise
//    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m15(0, 1, 4, orientation); // kmer, start, end, strand
//    MiniRecord mr15(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m15, mr15);
//    hits.insert(mh);
//    Minimizer m16(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr16(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m16, mr16);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(31, 33) };
//    prg_path.initialize(d);
//    Minimizer m17(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr17(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m17, mr17);
//    hits.insert(mh);
//    d = { Interval(78, 81) };
//    prg_path.initialize(d);
//    Minimizer m18(0, 13, 16, orientation); // kmer, start, end, strand
//    MiniRecord mr18(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m18, mr18);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // READ 2
//    read_id = 2;
//    pr = make_shared<pangenome::Read>(read_id);
//
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m19(0, 4, 7, orientation); // kmer, start, end, strand
//    MiniRecord mr19(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m19, mr19);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m20(0, 17, 20, orientation); // kmer, start, end, strand
//    MiniRecord mr20(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m20, mr20);
//    hits.insert(mh);
//    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
//    prg_path.initialize(d);
//    Minimizer m21(0, 15, 18, orientation); // kmer, start, end, strand
//    MiniRecord mr21(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m21, mr21);
//    hits.insert(mh);
//
//    // hits on path
//    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
//    prg_path.initialize(d);
//    Minimizer m22(0, 5, 8, orientation); // kmer, start, end, strand
//    MiniRecord mr22(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m22, mr22);
//    hits.insert(mh);
//    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m23(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr23(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m23, mr23);
//    hits.insert(mh);
//    d = { Interval(16, 17), Interval(27, 29) };
//    prg_path.initialize(d);
//    Minimizer m24(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr24(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m24, mr24);
//    hits.insert(mh);
//    d = { Interval(27, 30) };
//    prg_path.initialize(d);
//    Minimizer m25(0, 10, 13, orientation); // kmer, start, end, strand
//    MiniRecord mr25(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m25, mr25);
//    hits.insert(mh);
//
//    // noise
//    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m26(0, 1, 4, orientation); // kmer, start, end, strand
//    MiniRecord mr26(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m26, mr26);
//    hits.insert(mh);
//    Minimizer m27(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr27(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m27, mr27);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(31, 33) };
//    prg_path.initialize(d);
//    Minimizer m28(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr28(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m28, mr28);
//    hits.insert(mh);
//    d = { Interval(78, 81) };
//    prg_path.initialize(d);
//    Minimizer m29(0, 13, 16, orientation); // kmer, start, end, strand
//    MiniRecord mr29(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m29, mr29);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // READ 3
//    read_id = 3;
//    pr = make_shared<pangenome::Read>(read_id);
//
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m30(0, 4, 7, orientation); // kmer, start, end, strand
//    MiniRecord mr30(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m30, mr30);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m31(0, 10, 13, orientation); // kmer, start, end, strand
//    MiniRecord mr31(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m31, mr31);
//    hits.insert(mh);
//    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
//    prg_path.initialize(d);
//    Minimizer m32(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr32(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m32, mr32);
//    hits.insert(mh);
//
//    // hits on path
//    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m33(0, 6, 9, orientation); // kmer, start, end, strand
//    MiniRecord mr33(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m33, mr33);
//    hits.insert(mh);
//    d = { Interval(16, 17), Interval(27, 29) };
//    prg_path.initialize(d);
//    Minimizer m34(0, 7, 10, orientation); // kmer, start, end, strand
//    MiniRecord mr34(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m34, mr34);
//    hits.insert(mh);
//
//
//    // noise
//    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m35(0, 1, 4, orientation); // kmer, start, end, strand
//    MiniRecord mr35(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m35, mr35);
//    hits.insert(mh);
//    Minimizer m36(0, 7, 10, orientation); // kmer, start, end, strand
//    MiniRecord mr36(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m36, mr36);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // READ 4
//    read_id = 4;
//    pr = make_shared<pangenome::Read>(read_id);
//
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m37(0, 4, 7, orientation); // kmer, start, end, strand
//    MiniRecord mr37(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m37, mr37);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m38(0, 17, 20, orientation); // kmer, start, end, strand
//    MiniRecord mr38(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m38, mr38);
//    hits.insert(mh);
//
//    // noise
//    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m39(0, 1, 4, orientation); // kmer, start, end, strand
//    MiniRecord mr39(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m39, mr39);
//    hits.insert(mh);
//    Minimizer m40(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr40(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m40, mr40);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(31, 33) };
//    prg_path.initialize(d);
//    Minimizer m41(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr41(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m41, mr41);
//    hits.insert(mh);
//    d = { Interval(78, 81) };
//    prg_path.initialize(d);
//    Minimizer m42(0, 13, 16, orientation); // kmer, start, end, strand
//    MiniRecord mr42(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m42, mr42);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // RUN GET_READ_OVERLAPS
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp {//l3.prg.nodes[0],
//            l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6], l3.prg.nodes[7]//, l3.prg.nodes[9]
//    };
//    // A G C T CGG  TAT
//    const std::set<ReadCoordinate> expected_overlaps {{ 0, 3, 9,  1 },
//                                                      { 1, 7, 13, 1 },
//                                                      { 2, 5, 13, 1 },
//                                                      { 3, 6, 10, 1 }};
//
//    prg::Path local_path;
//    for (const auto &node : lmp) {
//        local_path.add_end_interval(node->pos);
//    }
//    const auto overlaps { pan_node->get_read_overlap_coordinates(local_path) };
//
//    EXPECT_ITERABLE_EQ(std::set<ReadCoordinate>, expected_overlaps, overlaps);
//}
//
//
//TEST(ExtractReadsTest, get_read_overlap_coordinates_no_duplicates) {
//    //
//    //  Read 0 has prg 3 sequence in interval (2,12] only
//    //  Read 1 has prg 3 sequence in interval (6,16] as well as noise
//    //  Read 2 has prg 3 sequence in interval (4,20] stretched out
//    //  Read 3 has prg 3 sequence in interval (4,14] but is missing bits
//    //  Read 4 doesn't have prg 3 sequence - on all hits are noise
//    //  Read 5 is a duplicate of Read 0
//    //
//    uint32_t read_id = 0, knode_id = 0;
//    bool orientation(true);
//    deque<Interval> d;
//    prg::Path prg_path;
//    MinimizerHitPtr mh;
//
//    uint32_t prg_id = 3;
//    auto local_prg_ptr { std::make_shared<LocalPRG>(prg_id, "three", "") };
//    PanNodePtr pan_node = make_shared<pangenome::Node>(local_prg_ptr);
//    PanReadPtr pr = make_shared<pangenome::Read>(read_id);
//
//    set<MinimizerHitPtr, pComp> hits;
//
//
//    // READ 0
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m1(0, 2, 5, orientation); // kmer, start, end, strand
//    MiniRecord mr1(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m1, mr1);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m2(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr2(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m2, mr2);
//    hits.insert(mh);
//    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
//    prg_path.initialize(d);
//    Minimizer m3(0, 7, 10, orientation); // kmer, start, end, strand
//    MiniRecord mr3(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m3, mr3);
//    hits.insert(mh);
//
//    // hits on path
//    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
//    prg_path.initialize(d);
//    Minimizer m4(0, 3, 6, orientation); // kmer, start, end, strand
//    MiniRecord mr4(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m4, mr4);
//    hits.insert(mh);
//    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m5(0, 4, 7, orientation); // kmer, start, end, strand
//    MiniRecord mr5(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m5, mr5);
//    hits.insert(mh);
//    d = { Interval(16, 17), Interval(27, 29) };
//    prg_path.initialize(d);
//    Minimizer m6(0, 5, 8, orientation); // kmer, start, end, strand
//    MiniRecord mr6(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m6, mr6);
//    hits.insert(mh);
//    d = { Interval(27, 30) };
//    prg_path.initialize(d);
//    Minimizer m7(0, 6, 9, orientation); // kmer, start, end, strand
//    MiniRecord mr7(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m7, mr7);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // READ 1
//    read_id = 1;
//    pr = make_shared<pangenome::Read>(read_id);
//
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m8(0, 6, 9, orientation); // kmer, start, end, strand
//    MiniRecord mr8(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m8, mr8);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m9(0, 12, 15, orientation); // kmer, start, end, strand
//    MiniRecord mr9(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m9 ,mr9);
//    hits.insert(mh);
//    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
//    prg_path.initialize(d);
//    Minimizer m10(0, 11, 14, orientation); // kmer, start, end, strand
//    MiniRecord mr10(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m10, mr10);
//    hits.insert(mh);
//
//    // hits on path
//    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
//    prg_path.initialize(d);
//    Minimizer m11(0, 7, 10, orientation); // kmer, start, end, strand
//    MiniRecord mr11(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m11, mr11);
//    hits.insert(mh);
//    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m12(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr12(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m12, mr12);
//    hits.insert(mh);
//    d = { Interval(16, 17), Interval(27, 29) };
//    prg_path.initialize(d);
//    Minimizer m13(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr13(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m13, mr13);
//    hits.insert(mh);
//    d = { Interval(27, 30) };
//    prg_path.initialize(d);
//    Minimizer m14(0, 10, 13, orientation); // kmer, start, end, strand
//    MiniRecord mr14(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m14, mr14);
//    hits.insert(mh);
//
//    // noise
//    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m15(0, 1, 4, orientation); // kmer, start, end, strand
//    MiniRecord mr15(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m15, mr15);
//    hits.insert(mh);
//    Minimizer m16(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr16(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m16, mr16);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(31, 33) };
//    prg_path.initialize(d);
//    Minimizer m17(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr17(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m17, mr17);
//    hits.insert(mh);
//    d = { Interval(78, 81) };
//    prg_path.initialize(d);
//    Minimizer m18(0, 13, 16, orientation); // kmer, start, end, strand
//    MiniRecord mr18(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m18, mr18);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // READ 2
//    read_id = 2;
//    pr = make_shared<pangenome::Read>(read_id);
//
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m19(0, 4, 7, orientation); // kmer, start, end, strand
//    MiniRecord mr19(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m19, mr19);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m20(0, 17, 20, orientation); // kmer, start, end, strand
//    MiniRecord mr20(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m20, mr20);
//    hits.insert(mh);
//    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
//    prg_path.initialize(d);
//    Minimizer m21(0, 15, 18, orientation); // kmer, start, end, strand
//    MiniRecord mr21(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m21, mr21);
//    hits.insert(mh);
//
//    // hits on path
//    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
//    prg_path.initialize(d);
//    Minimizer m22(0, 5, 8, orientation); // kmer, start, end, strand
//    MiniRecord mr22(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m22, mr22);
//    hits.insert(mh);
//    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m23(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr23(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m23, mr23);
//    hits.insert(mh);
//    d = { Interval(16, 17), Interval(27, 29) };
//    prg_path.initialize(d);
//    Minimizer m24(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr24(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m24, mr24);
//    hits.insert(mh);
//    d = { Interval(27, 30) };
//    prg_path.initialize(d);
//    Minimizer m25(0, 10, 13, orientation); // kmer, start, end, strand
//    MiniRecord mr25(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m25, mr25);
//    hits.insert(mh);
//
//    // noise
//    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m26(0, 1, 4, orientation); // kmer, start, end, strand
//    MiniRecord mr26(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m26, mr26);
//    hits.insert(mh);
//    Minimizer m27(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr27(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m27, mr27);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(31, 33) };
//    prg_path.initialize(d);
//    Minimizer m28(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr28(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m28, mr28);
//    hits.insert(mh);
//    d = { Interval(78, 81) };
//    prg_path.initialize(d);
//    Minimizer m29(0, 13, 16, orientation); // kmer, start, end, strand
//    MiniRecord mr29(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m29, mr29);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // READ 3
//    read_id = 3;
//    pr = make_shared<pangenome::Read>(read_id);
//
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m30(0, 4, 7, orientation); // kmer, start, end, strand
//    MiniRecord mr30(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m30, mr30);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m31(0, 10, 13, orientation); // kmer, start, end, strand
//    MiniRecord mr31(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m31, mr31);
//    hits.insert(mh);
//    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
//    prg_path.initialize(d);
//    Minimizer m32(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr32(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m32, mr32);
//    hits.insert(mh);
//
//    // hits on path
//    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m33(0, 6, 9, orientation); // kmer, start, end, strand
//    MiniRecord mr33(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m33, mr33);
//    hits.insert(mh);
//    d = { Interval(16, 17), Interval(27, 29) };
//    prg_path.initialize(d);
//    Minimizer m34(0, 7, 10, orientation); // kmer, start, end, strand
//    MiniRecord mr34(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m34, mr34);
//    hits.insert(mh);
//
//
//    // noise
//    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m35(0, 1, 4, orientation); // kmer, start, end, strand
//    MiniRecord mr35(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m35, mr35);
//    hits.insert(mh);
//    Minimizer m36(0, 7, 10, orientation); // kmer, start, end, strand
//    MiniRecord mr36(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m36, mr36);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // READ 4
//    read_id = 4;
//    pr = make_shared<pangenome::Read>(read_id);
//
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m37(0, 4, 7, orientation); // kmer, start, end, strand
//    MiniRecord mr37(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m37, mr37);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m38(0, 17, 20, orientation); // kmer, start, end, strand
//    MiniRecord mr38(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m38, mr38);
//    hits.insert(mh);
//
//    // noise
//    d = { Interval(7, 8), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m39(0, 1, 4, orientation); // kmer, start, end, strand
//    MiniRecord mr39(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m39, mr39);
//    hits.insert(mh);
//    Minimizer m40(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr40(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m40, mr40);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(31, 33) };
//    prg_path.initialize(d);
//    Minimizer m41(0, 9, 12, orientation); // kmer, start, end, strand
//    MiniRecord mr41(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m41, mr41);
//    hits.insert(mh);
//    d = { Interval(78, 81) };
//    prg_path.initialize(d);
//    Minimizer m42(0, 13, 16, orientation); // kmer, start, end, strand
//    MiniRecord mr42(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m42, mr42);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // READ 5
//    read_id = 0;
//    pr = make_shared<pangenome::Read>(read_id);
//
//    // hits overlapping edges of path
//    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
//    prg_path.initialize(d);
//    Minimizer m43(0, 2, 5, orientation); // kmer, start, end, strand
//    MiniRecord mr43(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m43, mr43);
//    hits.insert(mh);
//    d = { Interval(29, 30), Interval(33, 33), Interval(40, 42) };
//    prg_path.initialize(d);
//    Minimizer m44(0, 8, 11, orientation); // kmer, start, end, strand
//    MiniRecord mr44(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m44, mr44);
//    hits.insert(mh);
//    d = { Interval(28, 30), Interval(33, 33), Interval(40, 41) };
//    prg_path.initialize(d);
//    Minimizer m45(0, 7, 10, orientation); // kmer, start, end, strand
//    MiniRecord mr45(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m45, mr45);
//    hits.insert(mh);
//
//    // hits on path
//    d = { Interval(4, 5), Interval(8, 9), Interval(16, 17) };
//    prg_path.initialize(d);
//    Minimizer m46(0, 3, 6, orientation); // kmer, start, end, strand
//    MiniRecord mr46(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m46, mr46);
//    hits.insert(mh);
//    d = { Interval(8, 9), Interval(16, 17), Interval(27, 28) };
//    prg_path.initialize(d);
//    Minimizer m47(0, 4, 7, orientation); // kmer, start, end, strand
//    MiniRecord mr47(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m47, mr47);
//    hits.insert(mh);
//    d = { Interval(16, 17), Interval(27, 29) };
//    prg_path.initialize(d);
//    Minimizer m48(0, 5, 8, orientation); // kmer, start, end, strand
//    MiniRecord mr48(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m48, mr48);
//    hits.insert(mh);
//    d = { Interval(27, 30) };
//    prg_path.initialize(d);
//    Minimizer m49(0, 6, 9, orientation); // kmer, start, end, strand
//    MiniRecord mr49(prg_id, prg_path, knode_id, orientation);
//    mh = make_shared<MinimizerHit>(read_id, m49, mr49);
//    hits.insert(mh);
//
//    pr->add_hits(pan_node, hits);
//    pan_node->reads.insert(pr);
//    hits.clear();
//
//    // RUN GET_READ_OVERLAPS
//    LocalPRG l3(3, "nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
//    const std::vector<LocalNodePtr> lmp {//l3.prg.nodes[0],
//            l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4], l3.prg.nodes[6], l3.prg.nodes[7]//, l3.prg.nodes[9]
//    };
//    // A G C T CGG  TAT
//    const std::set<ReadCoordinate> expected_overlaps {{ 0, 3, 9,  1 },
//                                                      { 1, 7, 13, 1 },
//                                                      { 2, 5, 13, 1 },
//                                                      { 3, 6, 10, 1 }};
//
//    prg::Path local_path;
//    for (const auto &node : lmp) {
//        local_path.add_end_interval(node->pos);
//    }
//    const auto overlaps { pan_node->get_read_overlap_coordinates(local_path) };
//
//    EXPECT_ITERABLE_EQ(std::set<ReadCoordinate>, expected_overlaps, overlaps);
//}