#include <cstdint>
#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "pangenome/pannode.h"
#include "pangenome/pansample.h"
#include "pangenome/pangraph.h"
#include "pangenome/panread.h"
#include "minihit.h"
#include "localPRG.h"
#include "test_helpers.h"
#include "forward_declarations.h"

using namespace pangenome;

TEST(PangenomeNodeTest, create)
{
    auto local_graph_ptr { std::make_shared<LocalPRG>(4, "3", "") };
    pangenome::Node pan_node(local_graph_ptr, 3);
    uint32_t j = 3;
    EXPECT_EQ(j, pan_node.node_id);
    EXPECT_EQ((uint)4, pan_node.prg_id);
    EXPECT_EQ("3", pan_node.name);
    EXPECT_EQ((uint)0, pan_node.covg);
    EXPECT_EQ((uint)0, pan_node.reads.size());
    EXPECT_EQ((uint)0, pan_node.samples.size());
}

TEST(PangenomeNodeTest, get_name)
{
    auto l1 { std::make_shared<LocalPRG>(3, "3", "") };
    auto l2 { std::make_shared<LocalPRG>(2, "2", "") };
    pangenome::Node pn1(l1);
    pangenome::Node pn2(l2);
    pangenome::Node pn3(l2, 4);

    EXPECT_EQ(pn1.get_name(), "3");
    EXPECT_EQ(pn2.get_name(), "2");
    EXPECT_EQ(pn3.get_name(), "2.4");
}

TEST(PangenomeNodeTest, add_path)
{
    // setup the KmerGraph
    auto local_prg_ptr { std::make_shared<LocalPRG>(3, "3", "") };
    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13) };
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = { Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = { Interval(0, 1), Interval(19, 20), Interval(23, 24) };
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    d = { Interval(24, 24) };
    p.initialize(d);
    local_prg_ptr->kmer_prg.add_node(p);
    EXPECT_EQ((uint)7, local_prg_ptr->kmer_prg.nodes.size());

    // setup the Node
    pangenome::Node pn1(local_prg_ptr);
    std::vector<KmerNodePtr> kmp;
    pn1.add_path(kmp, 0);

    // do the tests
    EXPECT_EQ((uint)7, pn1.kmer_prg_with_coverage.kmer_prg->nodes.size());
    const auto& nodes = pn1.kmer_prg_with_coverage.kmer_prg->nodes;
    kmp = { nodes[0], nodes[3], nodes[4], nodes[6] };
    pn1.add_path(kmp, 0);
    EXPECT_EQ((uint)1, pn1.kmer_prg_with_coverage.get_reverse_covg(0, 0));
    EXPECT_EQ((uint)0, pn1.kmer_prg_with_coverage.get_reverse_covg(1, 0));
    EXPECT_EQ((uint)0, pn1.kmer_prg_with_coverage.get_reverse_covg(2, 0));
    EXPECT_EQ((uint)1, pn1.kmer_prg_with_coverage.get_reverse_covg(3, 0));
    EXPECT_EQ((uint)1, pn1.kmer_prg_with_coverage.get_reverse_covg(4, 0));
    EXPECT_EQ((uint)0, pn1.kmer_prg_with_coverage.get_reverse_covg(5, 0));
    EXPECT_EQ((uint)1, pn1.kmer_prg_with_coverage.get_reverse_covg(6, 0));
    EXPECT_EQ((uint)1, pn1.kmer_prg_with_coverage.get_forward_covg(0, 0));
    EXPECT_EQ((uint)0, pn1.kmer_prg_with_coverage.get_forward_covg(1, 0));
    EXPECT_EQ((uint)0, pn1.kmer_prg_with_coverage.get_forward_covg(2, 0));
    EXPECT_EQ((uint)1, pn1.kmer_prg_with_coverage.get_forward_covg(3, 0));
    EXPECT_EQ((uint)1, pn1.kmer_prg_with_coverage.get_forward_covg(4, 0));
    EXPECT_EQ((uint)0, pn1.kmer_prg_with_coverage.get_forward_covg(5, 0));
    EXPECT_EQ((uint)1, pn1.kmer_prg_with_coverage.get_forward_covg(6, 0));
}

class PangenomeNodeTest___construct_multisample_vcf___Fixture : public ::testing::Test {
protected:
    PangenomeNodeTest___construct_multisample_vcf___Fixture()
        : w(1)
        , k(3)
        , min_kmer_covg(0)
        , index(std::make_shared<Index>())
        , sample_names({ "sample1", "sample2", "sample3", "sample4" })
        , pangraph(sample_names)
        , master_vcf(create_VCF_with_default_parameters(0))
        , nested_varsite_PRG(std::make_shared<LocalPRG>(
              0, "nested varsite", "A 5 G 7 C 8 T 8 CT 7  6 G 5 T"))
        , modified_PRG(std::make_shared<LocalPRG>(
              1, "modified", "A 5 G 7 G 8 A 8 GA 7  6 G 5 T"))
    {
    }

    void SetUp() override
    {
        nested_varsite_kmer_coverage_graph_pointer = &nested_varsite_PRG->kmer_prg;
        nested_varsite_vcf_reference_path = { nested_varsite_PRG->prg.nodes[0],
            nested_varsite_PRG->prg.nodes[1], nested_varsite_PRG->prg.nodes[3],
            nested_varsite_PRG->prg.nodes[5], nested_varsite_PRG->prg.nodes[7] };

        modified_kmer_coverage_graph_pointer = &modified_PRG->kmer_prg;
        modified_vcf_reference_path = { modified_PRG->prg.nodes[0],
            modified_PRG->prg.nodes[1], modified_PRG->prg.nodes[3],
            modified_PRG->prg.nodes[5], modified_PRG->prg.nodes[7] };
    }

    void TearDown() override { }

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

    void assert_that_all_info_for_each_sample_and_allele_for_the_given_record_is_zero(
        const VCFRecord& record, uint32_t number_of_samples, uint32_t number_of_alleles)
    {
        EXPECT_EQ(number_of_samples, record.sampleIndex_to_sampleInfo.size());

        for (uint32_t sample_index = 0; sample_index < number_of_samples;
             ++sample_index) {
            EXPECT_EQ(number_of_alleles,
                record.sampleIndex_to_sampleInfo[sample_index].get_number_of_alleles());
            for (uint32_t allele_index = 0; allele_index < number_of_alleles;
                 ++allele_index) {
                EXPECT_EQ(0,
                    record.sampleIndex_to_sampleInfo[sample_index]
                        .get_mean_forward_coverage(allele_index));
                EXPECT_EQ(0,
                    record.sampleIndex_to_sampleInfo[sample_index]
                        .get_mean_reverse_coverage(allele_index));

                EXPECT_EQ(0,
                    record.sampleIndex_to_sampleInfo[sample_index]
                        .get_median_forward_coverage(allele_index));
                EXPECT_EQ(0,
                    record.sampleIndex_to_sampleInfo[sample_index]
                        .get_median_reverse_coverage(allele_index));

                EXPECT_EQ(0,
                    record.sampleIndex_to_sampleInfo[sample_index]
                        .get_sum_forward_coverage(allele_index));
                EXPECT_EQ(0,
                    record.sampleIndex_to_sampleInfo[sample_index]
                        .get_sum_reverse_coverage(allele_index));

                EXPECT_EQ(0,
                    record.sampleIndex_to_sampleInfo[sample_index].get_gaps(
                        allele_index));
            }
        }
    };
};

TEST_F(PangenomeNodeTest___construct_multisample_vcf___Fixture,
    construct_multisample_vcf_single_prg)
{
    nested_varsite_PRG->minimizer_sketch(index.get(), w, k);

    // sample1
    std::vector<KmerNodePtr> sample_kmer_path
        = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
              nested_varsite_kmer_coverage_graph_pointer->nodes[2],
              nested_varsite_kmer_coverage_graph_pointer->nodes[6],
              nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[0], sample_kmer_path);

    // sample2 identical to sample1
    sample_kmer_path = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
        nested_varsite_kmer_coverage_graph_pointer->nodes[2],
        nested_varsite_kmer_coverage_graph_pointer->nodes[6],
        nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[1], sample_kmer_path);

    // sample3 with top path
    sample_kmer_path = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
        nested_varsite_kmer_coverage_graph_pointer->nodes[1],
        nested_varsite_kmer_coverage_graph_pointer->nodes[5],
        nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[2], sample_kmer_path);

    // sample4 with bottom path
    sample_kmer_path = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
        nested_varsite_kmer_coverage_graph_pointer->nodes[4],
        nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[3], sample_kmer_path);

    auto& pannode = *pangraph.nodes[nested_varsite_PRG->id];
    pannode.construct_multisample_vcf(
        master_vcf, nested_varsite_vcf_reference_path, nested_varsite_PRG, w);

    EXPECT_EQ((uint)2, master_vcf.get_VCF_size());
    EXPECT_EQ((uint)4, master_vcf.samples.size());

    // NB samples order changes to get index of each sample so can compare
    // samples 1 and 2 are ref, sample 3 is top path and sample 4 is bottom path
    auto iter
        = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample1");
    auto sample1_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample2");
    auto sample2_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample3");
    auto sample3_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample4");
    auto sample4_index = std::distance(master_vcf.samples.begin(), iter);
    uint16_t alt_gt = 1;
    uint16_t ref_gt = 0;

    EXPECT_EQ((uint)1, master_vcf.get_records()[0]->get_pos());
    EXPECT_EQ("GT", master_vcf.get_records()[0]->get_ref());
    EXPECT_EQ((uint)1, master_vcf.get_records()[0]->get_alts().size());
    EXPECT_EQ("G", master_vcf.get_records()[0]->get_alts()[0]);
    EXPECT_EQ((uint)4, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());

    EXPECT_TRUE(master_vcf.get_records()[0]
                    ->sampleIndex_to_sampleInfo[sample1_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[0]
                  ->sampleIndex_to_sampleInfo[sample1_index]
                  .get_gt_from_max_likelihood_path(),
        ref_gt);
    EXPECT_TRUE(master_vcf.get_records()[0]
                    ->sampleIndex_to_sampleInfo[sample2_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[0]
                  ->sampleIndex_to_sampleInfo[sample2_index]
                  .get_gt_from_max_likelihood_path(),
        ref_gt);
    EXPECT_FALSE(master_vcf.get_records()[0]
                     ->sampleIndex_to_sampleInfo[sample3_index]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_TRUE(master_vcf.get_records()[0]
                    ->sampleIndex_to_sampleInfo[sample4_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[0]
                  ->sampleIndex_to_sampleInfo[sample4_index]
                  .get_gt_from_max_likelihood_path(),
        alt_gt);

    assert_that_all_info_for_each_sample_and_allele_for_the_given_record_is_zero(
        *(master_vcf.get_records()[0]), sample_names.size(), 2);

    EXPECT_EQ((uint)2, master_vcf.get_records()[1]->get_pos());
    EXPECT_EQ("T", master_vcf.get_records()[1]->get_ref());
    EXPECT_EQ((uint)2, master_vcf.get_records()[1]->get_alts().size());
    EXPECT_EQ("C", master_vcf.get_records()[1]->get_alts()[0]);
    EXPECT_EQ("CT", master_vcf.get_records()[1]->get_alts()[1]);
    EXPECT_EQ((uint)4, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());

    EXPECT_TRUE(master_vcf.get_records()[1]
                    ->sampleIndex_to_sampleInfo[sample1_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[1]
                  ->sampleIndex_to_sampleInfo[sample1_index]
                  .get_gt_from_max_likelihood_path(),
        ref_gt);
    EXPECT_TRUE(master_vcf.get_records()[1]
                    ->sampleIndex_to_sampleInfo[sample2_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[1]
                  ->sampleIndex_to_sampleInfo[sample2_index]
                  .get_gt_from_max_likelihood_path(),
        ref_gt);
    EXPECT_TRUE(master_vcf.get_records()[1]
                    ->sampleIndex_to_sampleInfo[sample3_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[1]
                  ->sampleIndex_to_sampleInfo[sample3_index]
                  .get_gt_from_max_likelihood_path(),
        alt_gt);
    EXPECT_FALSE(master_vcf.get_records()[1]
                     ->sampleIndex_to_sampleInfo[sample4_index]
                     .is_gt_from_max_likelihood_path_valid());

    assert_that_all_info_for_each_sample_and_allele_for_the_given_record_is_zero(
        *(master_vcf.get_records()[1]), sample_names.size(), 3);
}

TEST_F(PangenomeNodeTest___construct_multisample_vcf___Fixture,
    construct_multisample_vcf_two_prg)
{
    nested_varsite_PRG->minimizer_sketch(index.get(), w, k);
    modified_PRG->minimizer_sketch(index.get(), w, k);

    // sample1
    std::vector<KmerNodePtr> sample_kmer_path
        = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
              nested_varsite_kmer_coverage_graph_pointer->nodes[2],
              nested_varsite_kmer_coverage_graph_pointer->nodes[6],
              nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[0], sample_kmer_path);
    sample_kmer_path = { modified_kmer_coverage_graph_pointer->nodes[0],
        modified_kmer_coverage_graph_pointer->nodes[1],
        modified_kmer_coverage_graph_pointer->nodes[5],
        modified_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(modified_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        modified_PRG->id, sample_names[0], sample_kmer_path);

    // sample2 identical to sample1 in prg1, no prg2
    sample_kmer_path = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
        nested_varsite_kmer_coverage_graph_pointer->nodes[2],
        nested_varsite_kmer_coverage_graph_pointer->nodes[6],
        nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[1], sample_kmer_path);

    // sample3 with top path
    sample_kmer_path = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
        nested_varsite_kmer_coverage_graph_pointer->nodes[1],
        nested_varsite_kmer_coverage_graph_pointer->nodes[5],
        nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[2], sample_kmer_path);
    sample_kmer_path = { modified_kmer_coverage_graph_pointer->nodes[0],
        modified_kmer_coverage_graph_pointer->nodes[4],
        modified_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(modified_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        modified_PRG->id, sample_names[2], sample_kmer_path);

    // sample4 with bottom path
    sample_kmer_path = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
        nested_varsite_kmer_coverage_graph_pointer->nodes[4],
        nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[3], sample_kmer_path);
    sample_kmer_path = { modified_kmer_coverage_graph_pointer->nodes[0],
        modified_kmer_coverage_graph_pointer->nodes[3],
        modified_kmer_coverage_graph_pointer->nodes[7],
        modified_kmer_coverage_graph_pointer->nodes[8],
        modified_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(modified_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        modified_PRG->id, sample_names[3], sample_kmer_path);

    auto& pannode1 = *pangraph.nodes[nested_varsite_PRG->id];
    pannode1.construct_multisample_vcf(
        master_vcf, nested_varsite_vcf_reference_path, nested_varsite_PRG, w);
    auto& pannode2 = *pangraph.nodes[modified_PRG->id];
    pannode2.construct_multisample_vcf(
        master_vcf, modified_vcf_reference_path, modified_PRG, w);

    EXPECT_EQ((uint)4, master_vcf.get_VCF_size());
    EXPECT_EQ((uint)4, master_vcf.samples.size());

    // NB samples order changes to get index of each sample so can compare
    // samples 1 and 2 are ref, sample 3 is top path and sample 4 is bottom path
    auto iter
        = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample1");
    auto sample1_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample2");
    auto sample2_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample3");
    auto sample3_index = std::distance(master_vcf.samples.begin(), iter);
    iter = std::find(master_vcf.samples.begin(), master_vcf.samples.end(), "sample4");
    auto sample4_index = std::distance(master_vcf.samples.begin(), iter);
    uint16_t alt_gt = 1;
    uint16_t ref_gt = 0;
    uint16_t alt2_gt = 2;

    EXPECT_EQ((uint)1, master_vcf.get_records()[0]->get_pos());
    EXPECT_EQ("GT", master_vcf.get_records()[0]->get_ref());
    EXPECT_EQ((uint)1, master_vcf.get_records()[0]->get_alts().size());
    EXPECT_EQ("G", master_vcf.get_records()[0]->get_alts()[0]);
    EXPECT_EQ((uint)4, master_vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());

    EXPECT_TRUE(master_vcf.get_records()[0]
                    ->sampleIndex_to_sampleInfo[sample1_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[0]
                  ->sampleIndex_to_sampleInfo[sample1_index]
                  .get_gt_from_max_likelihood_path(),
        ref_gt);
    EXPECT_TRUE(master_vcf.get_records()[0]
                    ->sampleIndex_to_sampleInfo[sample2_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[0]
                  ->sampleIndex_to_sampleInfo[sample2_index]
                  .get_gt_from_max_likelihood_path(),
        ref_gt);
    EXPECT_FALSE(master_vcf.get_records()[0]
                     ->sampleIndex_to_sampleInfo[sample3_index]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_TRUE(master_vcf.get_records()[0]
                    ->sampleIndex_to_sampleInfo[sample4_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[0]
                  ->sampleIndex_to_sampleInfo[sample4_index]
                  .get_gt_from_max_likelihood_path(),
        alt_gt);

    assert_that_all_info_for_each_sample_and_allele_for_the_given_record_is_zero(
        *(master_vcf.get_records()[0]), 4, 2);

    EXPECT_EQ((uint)2, master_vcf.get_records()[1]->get_pos());
    EXPECT_EQ("T", master_vcf.get_records()[1]->get_ref());
    EXPECT_EQ((uint)2, master_vcf.get_records()[1]->get_alts().size());
    EXPECT_EQ("C", master_vcf.get_records()[1]->get_alts()[0]);
    EXPECT_EQ("CT", master_vcf.get_records()[1]->get_alts()[1]);
    EXPECT_EQ((uint)4, master_vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());

    EXPECT_TRUE(master_vcf.get_records()[1]
                    ->sampleIndex_to_sampleInfo[sample1_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[1]
                  ->sampleIndex_to_sampleInfo[sample1_index]
                  .get_gt_from_max_likelihood_path(),
        ref_gt);
    EXPECT_TRUE(master_vcf.get_records()[1]
                    ->sampleIndex_to_sampleInfo[sample2_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[1]
                  ->sampleIndex_to_sampleInfo[sample2_index]
                  .get_gt_from_max_likelihood_path(),
        ref_gt);
    EXPECT_TRUE(master_vcf.get_records()[1]
                    ->sampleIndex_to_sampleInfo[sample3_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[1]
                  ->sampleIndex_to_sampleInfo[sample3_index]
                  .get_gt_from_max_likelihood_path(),
        alt_gt);
    EXPECT_FALSE(master_vcf.get_records()[1]
                     ->sampleIndex_to_sampleInfo[sample4_index]
                     .is_gt_from_max_likelihood_path_valid());

    assert_that_all_info_for_each_sample_and_allele_for_the_given_record_is_zero(
        *(master_vcf.get_records()[1]), 4, 3);

    EXPECT_EQ((uint)1, master_vcf.get_records()[2]->get_pos());
    EXPECT_EQ("GA", master_vcf.get_records()[2]->get_ref());
    EXPECT_EQ((uint)1, master_vcf.get_records()[2]->get_alts().size());
    EXPECT_EQ("G", master_vcf.get_records()[2]->get_alts()[0]);
    EXPECT_EQ((uint)4, master_vcf.get_records()[2]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(master_vcf.get_records()[2]
                     ->sampleIndex_to_sampleInfo[sample1_index]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_FALSE(master_vcf.get_records()[2]
                     ->sampleIndex_to_sampleInfo[sample2_index]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_TRUE(master_vcf.get_records()[2]
                    ->sampleIndex_to_sampleInfo[sample3_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[2]
                  ->sampleIndex_to_sampleInfo[sample3_index]
                  .get_gt_from_max_likelihood_path(),
        alt_gt);
    EXPECT_FALSE(master_vcf.get_records()[2]
                     ->sampleIndex_to_sampleInfo[sample4_index]
                     .is_gt_from_max_likelihood_path_valid());

    assert_that_all_info_for_each_sample_and_allele_for_the_given_record_is_zero(
        *(master_vcf.get_records()[2]), 4, 2);

    EXPECT_EQ((uint)2, master_vcf.get_records()[3]->get_pos());
    EXPECT_EQ("A", master_vcf.get_records()[3]->get_ref());
    EXPECT_EQ((uint)2, master_vcf.get_records()[3]->get_alts().size());
    EXPECT_EQ("G", master_vcf.get_records()[3]->get_alts()[0]);
    EXPECT_EQ("GA", master_vcf.get_records()[3]->get_alts()[1]);
    EXPECT_EQ((uint)4, master_vcf.get_records()[3]->sampleIndex_to_sampleInfo.size());

    EXPECT_TRUE(master_vcf.get_records()[3]
                    ->sampleIndex_to_sampleInfo[sample1_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[3]
                  ->sampleIndex_to_sampleInfo[sample1_index]
                  .get_gt_from_max_likelihood_path(),
        alt_gt);
    EXPECT_FALSE(master_vcf.get_records()[3]
                     ->sampleIndex_to_sampleInfo[sample2_index]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_FALSE(master_vcf.get_records()[3]
                     ->sampleIndex_to_sampleInfo[sample3_index]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_TRUE(master_vcf.get_records()[3]
                    ->sampleIndex_to_sampleInfo[sample4_index]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(master_vcf.get_records()[3]
                  ->sampleIndex_to_sampleInfo[sample4_index]
                  .get_gt_from_max_likelihood_path(),
        alt2_gt);

    assert_that_all_info_for_each_sample_and_allele_for_the_given_record_is_zero(
        *(master_vcf.get_records()[3]), 4, 3);
}

TEST_F(PangenomeNodeTest___construct_multisample_vcf___Fixture,
    construct_multisample_vcf_two_prg_with_covgs)
{
    nested_varsite_PRG->minimizer_sketch(index.get(), w, k);
    modified_PRG->minimizer_sketch(index.get(), w, k);

    // sample1
    std::vector<KmerNodePtr> sample_kmer_path
        = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
              nested_varsite_kmer_coverage_graph_pointer->nodes[2],
              nested_varsite_kmer_coverage_graph_pointer->nodes[6],
              nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[0], sample_kmer_path);
    sample_kmer_path = { modified_kmer_coverage_graph_pointer->nodes[0],
        modified_kmer_coverage_graph_pointer->nodes[1],
        modified_kmer_coverage_graph_pointer->nodes[5],
        modified_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(modified_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        modified_PRG->id, sample_names[0], sample_kmer_path);

    auto& pannode1 = *pangraph.nodes[nested_varsite_PRG->id];
    auto& pannode2 = *pangraph.nodes[modified_PRG->id];

    pannode1.kmer_prg_with_coverage.set_forward_covg(0, 4, 0);
    pannode1.kmer_prg_with_coverage.set_forward_covg(2, 4, 0);
    pannode1.kmer_prg_with_coverage.set_forward_covg(6, 4, 0);
    pannode1.kmer_prg_with_coverage.set_forward_covg(9, 4, 0);
    pannode2.kmer_prg_with_coverage.set_forward_covg(0, 4, 0);
    pannode2.kmer_prg_with_coverage.set_forward_covg(1, 4, 0);
    pannode2.kmer_prg_with_coverage.set_forward_covg(5, 4, 0);
    pannode2.kmer_prg_with_coverage.set_forward_covg(9, 4, 0);

    // sample2 identical to sample1 in prg1, no prg2
    sample_kmer_path = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
        nested_varsite_kmer_coverage_graph_pointer->nodes[2],
        nested_varsite_kmer_coverage_graph_pointer->nodes[6],
        nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[1], sample_kmer_path);
    pannode1.kmer_prg_with_coverage.set_forward_covg(0, 10, 1);
    pannode1.kmer_prg_with_coverage.set_forward_covg(2, 10, 1);
    pannode1.kmer_prg_with_coverage.set_forward_covg(6, 10, 1);
    pannode1.kmer_prg_with_coverage.set_forward_covg(9, 10, 1);

    // sample3 with top path
    sample_kmer_path = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
        nested_varsite_kmer_coverage_graph_pointer->nodes[1],
        nested_varsite_kmer_coverage_graph_pointer->nodes[5],
        nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[2], sample_kmer_path);
    sample_kmer_path = { modified_kmer_coverage_graph_pointer->nodes[0],
        modified_kmer_coverage_graph_pointer->nodes[4],
        modified_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(modified_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        modified_PRG->id, sample_names[2], sample_kmer_path);
    pannode1.kmer_prg_with_coverage.set_forward_covg(0, 2, 2);
    pannode1.kmer_prg_with_coverage.set_forward_covg(1, 2, 2);
    pannode1.kmer_prg_with_coverage.set_forward_covg(5, 2, 2);
    pannode1.kmer_prg_with_coverage.set_forward_covg(9, 2, 2);
    pannode2.kmer_prg_with_coverage.set_forward_covg(0, 2, 2);
    pannode2.kmer_prg_with_coverage.set_forward_covg(4, 2, 2);
    pannode2.kmer_prg_with_coverage.set_forward_covg(9, 2, 2);

    // sample4 with bottom path
    sample_kmer_path = { nested_varsite_kmer_coverage_graph_pointer->nodes[0],
        nested_varsite_kmer_coverage_graph_pointer->nodes[4],
        nested_varsite_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(nested_varsite_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        nested_varsite_PRG->id, sample_names[3], sample_kmer_path);
    sample_kmer_path = { modified_kmer_coverage_graph_pointer->nodes[0],
        modified_kmer_coverage_graph_pointer->nodes[3],
        modified_kmer_coverage_graph_pointer->nodes[7],
        modified_kmer_coverage_graph_pointer->nodes[8],
        modified_kmer_coverage_graph_pointer->nodes[9] };
    pangraph.add_node(modified_PRG);
    pangraph.add_hits_between_PRG_and_sample(
        modified_PRG->id, sample_names[3], sample_kmer_path);
    pannode1.kmer_prg_with_coverage.set_forward_covg(0, 5, 3);
    pannode1.kmer_prg_with_coverage.set_forward_covg(4, 5, 3);
    pannode1.kmer_prg_with_coverage.set_forward_covg(9, 5, 3);
    pannode2.kmer_prg_with_coverage.set_forward_covg(0, 5, 3);
    pannode2.kmer_prg_with_coverage.set_forward_covg(3, 5, 3);
    pannode2.kmer_prg_with_coverage.set_forward_covg(7, 5, 3);
    pannode2.kmer_prg_with_coverage.set_forward_covg(8, 5, 3);
    pannode2.kmer_prg_with_coverage.set_forward_covg(9, 5, 3);

    VCF master_vcf = create_VCF_with_default_parameters(0);
    std::vector<LocalNodePtr> vcf_reference_path1 = { nested_varsite_PRG->prg.nodes[0],
        nested_varsite_PRG->prg.nodes[1], nested_varsite_PRG->prg.nodes[3],
        nested_varsite_PRG->prg.nodes[5], nested_varsite_PRG->prg.nodes[7] };
    std::vector<LocalNodePtr> vcf_reference_path2 = { modified_PRG->prg.nodes[0],
        modified_PRG->prg.nodes[1], modified_PRG->prg.nodes[3],
        modified_PRG->prg.nodes[5], modified_PRG->prg.nodes[7] };

    pannode1.construct_multisample_vcf(
        master_vcf, vcf_reference_path1, nested_varsite_PRG, w);
    pannode2.construct_multisample_vcf(
        master_vcf, vcf_reference_path2, modified_PRG, w);

    EXPECT_EQ((uint)4, master_vcf.get_VCF_size());
    EXPECT_EQ((uint)4, master_vcf.samples.size());

    // NB samples order changes to get index of each sample so can compare
    // samples 1 and 2 are ref, sample 3 is top path and sample 4 is bottom path
    auto sample1_index = master_vcf.get_sample_index("sample1");
    auto sample2_index = master_vcf.get_sample_index("sample2");
    auto sample3_index = master_vcf.get_sample_index("sample3");
    auto sample4_index = master_vcf.get_sample_index("sample4");
    std::vector<uint16_t> covgs_40 = { 4, 0 };
    std::vector<uint16_t> covgs_100 = { 10, 0 };
    std::vector<uint16_t> covgs_02 = { 0, 2 };
    std::vector<uint16_t> covgs_05 = { 0, 5 };
    std::vector<uint16_t> covgs_00 = { 0, 0 };
    std::vector<uint16_t> covgs_400 = { 4, 0, 0 };
    std::vector<uint16_t> covgs_040 = { 0, 4, 0 };
    std::vector<uint16_t> covgs_1000 = { 10, 0, 0 };
    std::vector<uint16_t> covgs_020 = { 0, 2, 0 };
    std::vector<uint16_t> covgs_005 = { 0, 0, 5 };
    std::vector<uint16_t> covgs_000 = { 0, 0, 0 };

    auto check_coverages = [](const SampleInfo& sample_info, size_t number_of_alleles,
                               const std::vector<uint16_t>& mean_fwd_coverages,
                               const std::vector<uint16_t>& mean_rev_coverages) {
        EXPECT_EQ(number_of_alleles, sample_info.get_number_of_alleles());

        for (size_t i = 0; i < number_of_alleles; ++i) {
            EXPECT_EQ(mean_fwd_coverages[i], sample_info.get_mean_forward_coverage(i));
            EXPECT_EQ(mean_rev_coverages[i], sample_info.get_mean_reverse_coverage(i));
        }
    };

    // sample1
    // EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_FWD_COVG"],
    // covgs_40); EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_REV_COVG"],
    // covgs_00);
    check_coverages(
        master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample1_index], 2,
        covgs_40, covgs_00);

    // EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_FWD_COVG"],
    // covgs_400); EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_REV_COVG"],
    // covgs_000);
    check_coverages(
        master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample1_index], 3,
        covgs_400, covgs_000);

    check_coverages(
        master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample1_index], 2,
        covgs_00, covgs_00);

    // EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_FWD_COVG"],
    // covgs_040); EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index]["MEAN_REV_COVG"],
    // covgs_000);
    check_coverages(
        master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample1_index], 3,
        covgs_040, covgs_000);

    // sample2
    // EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample2_index]["MEAN_FWD_COVG"],
    // covgs_100); EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample2_index]["MEAN_REV_COVG"],
    // covgs_00);
    check_coverages(
        master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample2_index], 2,
        covgs_100, covgs_00);

    // EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample2_index]["MEAN_FWD_COVG"],
    // covgs_1000); EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample2_index]["MEAN_REV_COVG"],
    // covgs_000);
    check_coverages(
        master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample2_index], 3,
        covgs_1000, covgs_000);

    check_coverages(
        master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample2_index], 2,
        covgs_00, covgs_00);
    check_coverages(
        master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample2_index], 3,
        covgs_000, covgs_000);

    // sample3
    check_coverages(
        master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample3_index], 2,
        covgs_00, covgs_00);

    // EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample3_index]["MEAN_FWD_COVG"],
    // covgs_020);
    check_coverages(
        master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample3_index], 3,
        covgs_020, covgs_000);

    // EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample3_index]["MEAN_FWD_COVG"],
    // covgs_02); EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample3_index]["MEAN_REV_COVG"],
    // covgs_00);
    check_coverages(
        master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample3_index], 2,
        covgs_02, covgs_00);

    check_coverages(
        master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample3_index], 3,
        covgs_000, covgs_000);

    // sample4
    // EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample4_index]["MEAN_FWD_COVG"],
    // covgs_05); EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample4_index]["MEAN_REV_COVG"],
    // covgs_00);
    check_coverages(
        master_vcf.get_records()[0]->sampleIndex_to_sampleInfo[sample4_index], 2,
        covgs_05, covgs_00);

    // EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample4_index]["MEAN_REV_COVG"],
    // covgs_000);
    check_coverages(
        master_vcf.get_records()[1]->sampleIndex_to_sampleInfo[sample4_index], 3,
        covgs_000, covgs_000);

    check_coverages(
        master_vcf.get_records()[2]->sampleIndex_to_sampleInfo[sample4_index], 2,
        covgs_00, covgs_00);

    // EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index]["MEAN_FWD_COVG"],
    // covgs_005); EXPECT_ITERABLE_EQ(std::vector<uint16_t>,
    // master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index]["MEAN_REV_COVG"],
    // covgs_000);
    check_coverages(
        master_vcf.get_records()[3]->sampleIndex_to_sampleInfo[sample4_index], 3,
        covgs_005, covgs_000);
}

TEST(PangenomeNodeTest, equals)
{
    auto l1 { std::make_shared<LocalPRG>(3, "3", "") };
    auto l2 { std::make_shared<LocalPRG>(2, "2", "") };
    pangenome::Node pn1(l1);
    pangenome::Node pn2(l2);
    pangenome::Node pn3(l2);

    EXPECT_EQ(pn1, pn1);
    EXPECT_EQ(pn2, pn2);
    EXPECT_EQ(pn3, pn3);
    EXPECT_EQ(pn2, pn3);
    EXPECT_EQ(pn3, pn2);
    EXPECT_EQ((pn1 == pn2), false);
    EXPECT_EQ((pn1 == pn3), false);
}

TEST(PangenomeNodeTest, nequals)
{
    auto l1 { std::make_shared<LocalPRG>(3, "3", "") };
    auto l2 { std::make_shared<LocalPRG>(2, "2", "") };
    pangenome::Node pn1(l1);
    pangenome::Node pn2(l2);
    pangenome::Node pn3(l2);

    EXPECT_EQ((pn1 != pn2), true);
    EXPECT_EQ((pn2 != pn1), true);
    EXPECT_EQ((pn1 != pn1), false);
    EXPECT_EQ((pn2 != pn2), false);
    EXPECT_EQ((pn3 != pn3), false);
    EXPECT_EQ((pn2 != pn3), false);
}

TEST(PangenomeNodeTest, less)
{
    auto l1 { std::make_shared<LocalPRG>(3, "3", "") };
    auto l2 { std::make_shared<LocalPRG>(2, "2", "") };
    pangenome::Node pn1(l1);
    pangenome::Node pn2(l2);
    pangenome::Node pn3(l2);

    EXPECT_EQ((pn1 < pn1), false);
    EXPECT_EQ((pn2 < pn2), false);
    EXPECT_EQ((pn3 < pn3), false);
    EXPECT_EQ((pn1 < pn3), false);
    EXPECT_EQ((pn1 < pn2), false);
    EXPECT_EQ((pn2 < pn1), true);
    EXPECT_EQ((pn3 < pn1), true);
}
