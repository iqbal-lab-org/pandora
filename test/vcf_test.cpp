#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "test_macro.cpp"
#include "vcf.h"
#include "vcfrecord.h"
#include "interval.h"
#include "localnode.h"
#include <stdint.h>
#include <iostream>
#include "test_helpers.h"
#include "utils.h"

using namespace std;
using ::testing::_;
using ::testing::Field;
using ::testing::InSequence;
using ::testing::Property;
using ::testing::Return;
using ::testing::ReturnRef;

TEST(VCFTest, add_record_with_values)
{

    VCF vcf = create_VCF_with_default_parameters();
    EXPECT_EQ((uint)0, vcf.get_VCF_size());
    vcf.add_record("chrom1", 5, "A", "G");
    EXPECT_EQ((uint)1, vcf.get_VCF_size());
}

TEST(VCFTest, add_record_twice_with_values)
{

    VCF vcf = create_VCF_with_default_parameters();
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 5, "A", "G");
    EXPECT_EQ((uint)1, vcf.get_VCF_size());
}

TEST(VCFTest, add_two_records_with_values)
{

    VCF vcf = create_VCF_with_default_parameters();
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    EXPECT_EQ((uint)2, vcf.get_VCF_size());
}

TEST(VCFTest, add_two_records_and_a_repeat_with_values)
{

    VCF vcf = create_VCF_with_default_parameters();
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 5, "A", "G");
    EXPECT_EQ((uint)2, vcf.get_VCF_size());
}

TEST(VCFTest, add_record_by_record)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr = VCFRecord(&vcf, "chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_or_update_record_restricted_to_the_given_samples(vr, empty);
    EXPECT_EQ((uint)1, vcf.get_VCF_size());
}

TEST(VCFTest, add_record_by_record_and_values)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr = VCFRecord(&vcf, "chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_or_update_record_restricted_to_the_given_samples(vr, empty);
    vcf.add_record("chrom1", 79, "C", "G");
    EXPECT_EQ((uint)1, vcf.get_VCF_size());
}

TEST(VCFTest, add_record_by_values_and_record)
{
    VCF vcf = create_VCF_with_default_parameters();
    vcf.add_record("chrom1", 79, "C", "G");
    VCFRecord vr = VCFRecord(&vcf, "chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_or_update_record_restricted_to_the_given_samples(vr, empty);
    EXPECT_EQ((uint)1, vcf.get_VCF_size());
}

TEST(VCFTest, add_record_by_record_returned_by_reference)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr = VCFRecord(&vcf, "chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    VCFRecord& ref_vr
        = vcf.add_or_update_record_restricted_to_the_given_samples(vr, empty);
    EXPECT_EQ(ref_vr.get_chrom(), "chrom1");
    EXPECT_EQ(ref_vr.get_pos(), (uint)79);
}

TEST(VCFTest, add_samples_empty)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    std::vector<std::string> samples;
    vcf.add_samples(samples);
    EXPECT_EQ(vcf.samples.size(), (uint)0);
    EXPECT_EQ(vcf.get_VCF_size(), (uint)0);
}

TEST(VCFTest, add_samples_simple)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    std::vector<std::string> samples = { "hello", "there", "people" };
    vcf.add_samples(samples);
    EXPECT_ITERABLE_EQ(std::vector<std::string>, samples, vcf.samples);
    EXPECT_EQ(vcf.get_VCF_size(), (uint)0);
}

TEST(VCFTest, add_samples_with_record)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 5, "A", "G");

    std::vector<std::string> samples = { "hello", "there", "people" };
    vcf.add_samples(samples);

    std::vector<std::string> exp_samples = { "sample", "hello", "there", "people" };

    EXPECT_ITERABLE_EQ(std::vector<std::string>, exp_samples, vcf.samples);
    EXPECT_EQ(vcf.get_VCF_size(), (uint)1);
    EXPECT_EQ(
        vcf.get_records()[0]->sampleIndex_to_sampleInfo.size(), exp_samples.size());
}

TEST(VCFTest, add_a_new_record_discovered_in_a_sample_and_genotype_it)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 46, "T", "TA");
    uint j = 1;
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint16_t)1,
        vcf.get_records()[1]
            ->sampleIndex_to_sampleInfo[0]
            .get_gt_from_max_likelihood_path());
    EXPECT_EQ(j, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[0]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(j, vcf.get_records()[2]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[2]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(j, vcf.get_records()[3]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[3]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());

    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 79, "C", "C");
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint16_t)1,
        vcf.get_records()[1]
            ->sampleIndex_to_sampleInfo[0]
            .get_gt_from_max_likelihood_path());
    EXPECT_EQ(j, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[0]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(j, vcf.get_records()[2]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint16_t)0,
        vcf.get_records()[2]
            ->sampleIndex_to_sampleInfo[0]
            .get_gt_from_max_likelihood_path());
    EXPECT_EQ(j, vcf.get_records()[3]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint16_t)0,
        vcf.get_records()[3]
            ->sampleIndex_to_sampleInfo[0]
            .get_gt_from_max_likelihood_path());
}

TEST(VCFTest, add_record_by_record_with_existing_sample)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord(&vcf, "chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    VCFRecord& ref_vr
        = vcf.add_or_update_record_restricted_to_the_given_samples(vr, empty);
    EXPECT_EQ(ref_vr.get_chrom(), "chrom1");
    EXPECT_EQ(ref_vr.get_pos(), (uint)79);
    EXPECT_EQ(ref_vr.sampleIndex_to_sampleInfo.size(), (uint)1);
}

TEST(VCFTest, add_record_by_record_with_same_existing_sample)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord(&vcf, "chrom1", 79, "C", "G");
    vr.sampleIndex_to_sampleInfo[0].set_gt_from_max_likelihood_path(1);
    std::vector<std::string> samples = { "sample" };
    VCFRecord& ref_vr
        = vcf.add_or_update_record_restricted_to_the_given_samples(vr, samples);
    EXPECT_EQ(ref_vr.get_chrom(), "chrom1");
    EXPECT_EQ(ref_vr.get_pos(), (uint)79);
    EXPECT_EQ(ref_vr.sampleIndex_to_sampleInfo.size(), (uint)1);
    EXPECT_ITERABLE_EQ(vector<std::string>, samples, vcf.samples);
    EXPECT_TRUE(
        ref_vr.sampleIndex_to_sampleInfo[0].is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(ref_vr.sampleIndex_to_sampleInfo[0].get_gt_from_max_likelihood_path(),
        (uint16_t)1);
}

TEST(VCFTest, add_record_by_record_with_different_existing_sample)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord(&vcf, "chrom1", 79, "C", "G");
    vr.sampleIndex_to_sampleInfo[0].set_gt_from_max_likelihood_path(1);
    std::vector<std::string> samples = { "sample1" };
    VCFRecord& ref_vr
        = vcf.add_or_update_record_restricted_to_the_given_samples(vr, samples);
    EXPECT_EQ(ref_vr.get_chrom(), "chrom1");
    EXPECT_EQ(ref_vr.get_pos(), (uint)79);
    EXPECT_EQ(ref_vr.sampleIndex_to_sampleInfo.size(), (uint)2);
    std::vector<std::string> exp_samples = { "sample", "sample1" };
    EXPECT_ITERABLE_EQ(vector<std::string>, exp_samples, vcf.samples);
    EXPECT_FALSE(
        ref_vr.sampleIndex_to_sampleInfo[0].is_gt_from_max_likelihood_path_valid());
    EXPECT_TRUE(
        ref_vr.sampleIndex_to_sampleInfo[1].is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(ref_vr.sampleIndex_to_sampleInfo[1].get_gt_from_max_likelihood_path(),
        (uint16_t)1);
}

TEST(VCFTest, set_sample_gt_to_ref_allele_for_records_in_the_interval)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_record("chrom2", 30, "C", "A");

    vcf.set_sample_gt_to_ref_allele_for_records_in_the_interval(
        "sample", "chrom1", 15, 78);
    EXPECT_EQ((uint)1, vcf.samples.size());
    EXPECT_EQ((uint)5, vcf.get_VCF_size());
    EXPECT_EQ((uint)1, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[0]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ((uint)1, vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint16_t)0,
        vcf.get_records()[1]
            ->sampleIndex_to_sampleInfo[0]
            .get_gt_from_max_likelihood_path());
    EXPECT_EQ((uint)1, vcf.get_records()[2]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[2]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ((uint)1, vcf.get_records()[3]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[3]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ((uint)1, vcf.get_records()[4]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[4]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());

    vcf.set_sample_gt_to_ref_allele_for_records_in_the_interval(
        "sample2", "chrom1", 5, 46);
    EXPECT_EQ((uint)2, vcf.samples.size());
    EXPECT_EQ((uint)5, vcf.get_VCF_size());
    EXPECT_EQ((uint)2, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint16_t)0,
        vcf.get_records()[0]
            ->sampleIndex_to_sampleInfo[1]
            .get_gt_from_max_likelihood_path());
    EXPECT_EQ((uint)2, vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[1]
                     ->sampleIndex_to_sampleInfo[1]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ((uint)2, vcf.get_records()[2]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[2]
                     ->sampleIndex_to_sampleInfo[1]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ((uint)2, vcf.get_records()[3]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[3]
                     ->sampleIndex_to_sampleInfo[1]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ((uint)2, vcf.get_records()[4]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[4]
                     ->sampleIndex_to_sampleInfo[1]
                     .is_gt_from_max_likelihood_path_valid());
}

TEST(VCFTest, reorder_add_record_and_sample)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample1", "chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample2", "chrom1", 79, "C", "C");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample1", "chrom1", 79, "C", "A");

    vcf.sort_records();

    EXPECT_EQ((uint)2, vcf.samples.size());
    EXPECT_EQ((uint)4, vcf.get_VCF_size());
    EXPECT_EQ((uint)2, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vcf.get_records()[2]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vcf.get_records()[3]->sampleIndex_to_sampleInfo.size());
    EXPECT_FALSE(vcf.get_records()[0]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ((uint16_t)1,
        vcf.get_records()[1]
            ->sampleIndex_to_sampleInfo[0]
            .get_gt_from_max_likelihood_path());
    EXPECT_EQ((uint16_t)1,
        vcf.get_records()[2]
            ->sampleIndex_to_sampleInfo[0]
            .get_gt_from_max_likelihood_path());
    EXPECT_FALSE(vcf.get_records()[3]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_FALSE(vcf.get_records()[0]
                     ->sampleIndex_to_sampleInfo[1]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_FALSE(vcf.get_records()[1]
                     ->sampleIndex_to_sampleInfo[1]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ((uint16_t)0,
        vcf.get_records()[2]
            ->sampleIndex_to_sampleInfo[1]
            .get_gt_from_max_likelihood_path());
    EXPECT_EQ((uint16_t)0,
        vcf.get_records()[3]
            ->sampleIndex_to_sampleInfo[1]
            .get_gt_from_max_likelihood_path());
}

TEST(VCFTest, append_vcf_simple_case)
{
    VCF vcf = create_VCF_with_default_parameters();
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    VCF new_vcf = create_VCF_with_default_parameters();
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom2", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_record("chrom2", 79, "C", "A");

    vcf.append_vcf(new_vcf);
    EXPECT_EQ((uint)8, vcf.get_VCF_size());
    for (uint i = 0; i < 4; ++i) {
        EXPECT_EQ(vcf.get_records()[i]->get_chrom(), "chrom1");
    }
    for (uint i = 4; i < 8; ++i) {
        EXPECT_EQ(vcf.get_records()[i]->get_chrom(), "chrom2");
    }
    EXPECT_EQ((uint)5, vcf.get_records()[4]->get_pos());
    EXPECT_EQ("TA", vcf.get_records()[5]->get_alts()[0]);
    EXPECT_EQ((uint)79, vcf.get_records()[6]->get_pos());
    EXPECT_EQ("A", vcf.get_records()[7]->get_alts()[0]);
}

TEST(VCFTest, append_vcf_some_duplicate_records)
{
    VCF vcf = create_VCF_with_default_parameters();
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    VCF new_vcf = create_VCF_with_default_parameters();
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom1", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_record("chrom1", 79, "C", "A");

    vcf.append_vcf(new_vcf);
    EXPECT_EQ((uint)6, vcf.get_VCF_size());
    for (uint i = 0; i < 4; ++i) {
        EXPECT_EQ(vcf.get_records()[i]->get_chrom(), "chrom1");
    }
    for (uint i = 4; i < 6; ++i) {
        EXPECT_EQ(vcf.get_records()[i]->get_chrom(), "chrom2");
    }
    EXPECT_EQ((uint)5, vcf.get_records()[4]->get_pos());
    EXPECT_EQ((uint)79, vcf.get_records()[5]->get_pos());
}

TEST(VCFTest, append_vcf_one_sample)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 79, "C", "G");

    VCF new_vcf = create_VCF_with_default_parameters(0);
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom1", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_record("chrom1", 79, "C", "A");

    vcf.append_vcf(new_vcf);
    EXPECT_EQ((uint)1, vcf.samples.size());
    EXPECT_EQ("sample", vcf.samples[0]);
    EXPECT_EQ((uint)1, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)1, vcf.get_records()[5]->sampleIndex_to_sampleInfo.size());
    bool valid_gt = vcf.get_records()[2]
                        ->sampleIndex_to_sampleInfo[0]
                        .is_gt_from_max_likelihood_path_valid();
    EXPECT_TRUE(valid_gt);
    EXPECT_EQ((uint)1,
        vcf.get_records()[2]
            ->sampleIndex_to_sampleInfo[0]
            .get_gt_from_max_likelihood_path());
    valid_gt = vcf.get_records()[0]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[1]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[4]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[3]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[5]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
}

TEST(VCFTest, append_vcf_one_sample_in_new_vcf)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    VCF new_vcf = create_VCF_with_default_parameters(0);
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom1", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_record("chrom1", 79, "C", "A");
    new_vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom2", 5, "A", "G");

    vcf.append_vcf(new_vcf);
    EXPECT_EQ((uint)1, vcf.samples.size());
    EXPECT_EQ("sample", vcf.samples[0]);
    EXPECT_EQ((uint)1, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)1, vcf.get_records()[5]->sampleIndex_to_sampleInfo.size());
    bool valid_gt = vcf.get_records()[4]
                        ->sampleIndex_to_sampleInfo[0]
                        .is_gt_from_max_likelihood_path_valid();
    EXPECT_TRUE(valid_gt);
    EXPECT_EQ((uint)1,
        vcf.get_records()[4]
            ->sampleIndex_to_sampleInfo[0]
            .get_gt_from_max_likelihood_path());
    valid_gt = vcf.get_records()[0]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[1]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[2]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[3]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[5]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
}

TEST(VCFTest, append_vcf_shared_sample)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 46, "T", "TA");

    VCF new_vcf = create_VCF_with_default_parameters(0);
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom1", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_record("chrom1", 79, "C", "A");
    new_vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 46, "T", "TA");

    vcf.append_vcf(new_vcf);
    EXPECT_EQ((uint)1, vcf.samples.size());
    EXPECT_EQ("sample", vcf.samples[0]);
    EXPECT_EQ((uint)1, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)1, vcf.get_records()[5]->sampleIndex_to_sampleInfo.size());
    bool valid_gt = vcf.get_records()[1]
                        ->sampleIndex_to_sampleInfo[0]
                        .is_gt_from_max_likelihood_path_valid();
    EXPECT_TRUE(valid_gt);
    EXPECT_EQ((uint)1,
        vcf.get_records()[1]
            ->sampleIndex_to_sampleInfo[0]
            .get_gt_from_max_likelihood_path());
    valid_gt = vcf.get_records()[0]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[4]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[2]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[3]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
    valid_gt = vcf.get_records()[5]
                   ->sampleIndex_to_sampleInfo[0]
                   .is_gt_from_max_likelihood_path_valid();
    EXPECT_FALSE(valid_gt);
}

TEST(VCFTest, append_vcf_shared_samples_different_order)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 46, "T", "TA");

    VCF new_vcf = create_VCF_with_default_parameters(0);
    new_vcf.add_record("chrom1", 79, "C", "A");
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom1", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample1", "chrom1", 46, "T", "T");
    new_vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample1", "chrom1", 79, "C", "A");

    vcf.append_vcf(new_vcf);

    EXPECT_EQ((uint)2, vcf.samples.size());
    vector<string> v = { "sample", "sample1" };
    EXPECT_ITERABLE_EQ(vector<string>, v, vcf.samples);
    EXPECT_EQ((uint)2, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vcf.get_records()[2]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vcf.get_records()[3]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vcf.get_records()[4]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vcf.get_records()[5]->sampleIndex_to_sampleInfo.size());

    uint16_t alt_gt = 1;
    uint16_t ref_gt = 0;

    EXPECT_FALSE(vcf.get_records()[0]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_FALSE(vcf.get_records()[0]
                     ->sampleIndex_to_sampleInfo[1]
                     .is_gt_from_max_likelihood_path_valid());

    EXPECT_TRUE(vcf.get_records()[1]
                    ->sampleIndex_to_sampleInfo[0]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_TRUE(vcf.get_records()[1]
                    ->sampleIndex_to_sampleInfo[1]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(vcf.get_records()[1]
                  ->sampleIndex_to_sampleInfo[0]
                  .get_gt_from_max_likelihood_path(),
        alt_gt);
    EXPECT_EQ(vcf.get_records()[1]
                  ->sampleIndex_to_sampleInfo[1]
                  .get_gt_from_max_likelihood_path(),
        ref_gt);

    EXPECT_FALSE(vcf.get_records()[2]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_FALSE(vcf.get_records()[2]
                     ->sampleIndex_to_sampleInfo[1]
                     .is_gt_from_max_likelihood_path_valid());

    EXPECT_FALSE(vcf.get_records()[3]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_TRUE(vcf.get_records()[3]
                    ->sampleIndex_to_sampleInfo[1]
                    .is_gt_from_max_likelihood_path_valid());
    EXPECT_EQ(vcf.get_records()[3]
                  ->sampleIndex_to_sampleInfo[1]
                  .get_gt_from_max_likelihood_path(),
        alt_gt);

    EXPECT_FALSE(vcf.get_records()[4]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_FALSE(vcf.get_records()[4]
                     ->sampleIndex_to_sampleInfo[1]
                     .is_gt_from_max_likelihood_path_valid());

    EXPECT_FALSE(vcf.get_records()[5]
                     ->sampleIndex_to_sampleInfo[0]
                     .is_gt_from_max_likelihood_path_valid());
    EXPECT_FALSE(vcf.get_records()[5]
                     ->sampleIndex_to_sampleInfo[1]
                     .is_gt_from_max_likelihood_path_valid());
}

TEST(VCFTest, sort_records)
{
    VCF vcf = create_VCF_with_default_parameters();
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_record("chrom2", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom2", 79, "C", "G");
    vcf.sort_records();

    EXPECT_EQ((uint)6, vcf.get_VCF_size());
    for (uint i = 0; i < 4; ++i) {
        EXPECT_EQ("chrom1", vcf.get_records()[i]->get_chrom());
    }
    for (uint i = 4; i < 6; ++i) {
        EXPECT_EQ("chrom2", vcf.get_records()[i]->get_chrom());
    }
    EXPECT_EQ((uint)5, vcf.get_records()[0]->get_pos());
    EXPECT_EQ((uint)5, vcf.get_records()[4]->get_pos());
    EXPECT_EQ((uint)46, vcf.get_records()[1]->get_pos());
    EXPECT_EQ((uint)79, vcf.get_records()[2]->get_pos());
    EXPECT_EQ((uint)79, vcf.get_records()[3]->get_pos());
    EXPECT_EQ((uint)79, vcf.get_records()[5]->get_pos());
    EXPECT_EQ("G", vcf.get_records()[3]->get_alts()[0]);
    EXPECT_EQ("G", vcf.get_records()[5]->get_alts()[0]);
}

TEST(VCFTest, pos_in_range)
{
    VCF vcf = create_VCF_with_default_parameters();
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom2", 20, "A", "G");
    vcf.add_record("chrom2", 79, "C", "G");
    // vcf.sort();

    EXPECT_TRUE(vcf.pos_in_range(4, 6, "chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(5, 6, "chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(4, 5, "chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(4, 6, "chrom2"));

    EXPECT_TRUE(vcf.pos_in_range(45, 47, "chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(46, 47, "chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(45, 46, "chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(45, 47, "chrom2"));

    EXPECT_TRUE(vcf.pos_in_range(78, 80, "chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(79, 80, "chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(78, 79, "chrom1"));
    EXPECT_TRUE(vcf.pos_in_range(78, 80, "chrom2"));
}

class VCFTest___genotype___Fixture : public ::testing::Test {
public:
    class VCFMock : public VCF {
    public:
        using VCF::records;
        using VCF::VCF;
        MOCK_METHOD(void, make_gt_compatible, (), (override));
    };

    class VCFRecordMock : public VCFRecord {
    public:
        using VCFRecord::VCFRecord;
        MOCK_METHOD(void, genotype_from_coverage, (), (override));
        MOCK_METHOD(void,
            genotype_from_coverage_using_maximum_likelihood_path_as_reference, (),
            (override));
        MOCK_METHOD(bool, is_SNP, (), (const override));
    };

    VCFTest___genotype___Fixture()
        : default_vcf(&default_genotyping_options)
        , snps_only_vcf(&genotyping_options_snps_only)
        , non_snp_vcf_record_ptr(new VCFRecordMock(&default_vcf))
        , snp_vcf_record_ptr(new VCFRecordMock(&default_vcf))
    {
    }

    void SetUp() override
    {
        ON_CALL(*non_snp_vcf_record_ptr, is_SNP).WillByDefault(Return(false));
        ON_CALL(*snp_vcf_record_ptr, is_SNP).WillByDefault(Return(true));

        default_vcf.records.push_back(non_snp_vcf_record_ptr);
        default_vcf.records.push_back(snp_vcf_record_ptr);

        snps_only_vcf.records.push_back(non_snp_vcf_record_ptr);
        snps_only_vcf.records.push_back(snp_vcf_record_ptr);
    }

    void TearDown() override {}

    VCFMock default_vcf;
    VCFMock snps_only_vcf;
    std::shared_ptr<VCFRecordMock> non_snp_vcf_record_ptr;
    std::shared_ptr<VCFRecordMock> snp_vcf_record_ptr;
    static GenotypingOptions genotyping_options_snps_only;
};
GenotypingOptions VCFTest___genotype___Fixture::genotyping_options_snps_only(
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 0.01, 0, 0, 0, 0, 0, 0, true);

TEST_F(VCFTest___genotype___Fixture, no_options_set___expects_death)
{
    EXPECT_DEATH(default_vcf.genotype(false, false), "");
}

TEST_F(VCFTest___genotype___Fixture, two_options_set___expects_death)
{
    EXPECT_DEATH(default_vcf.genotype(true, true), "");
}

TEST_F(VCFTest___genotype___Fixture, genotype_all_records___global_genotyping)
{
    {
        InSequence seq;
        EXPECT_CALL(*non_snp_vcf_record_ptr,
            genotype_from_coverage_using_maximum_likelihood_path_as_reference)
            .Times(1);
        EXPECT_CALL(*snp_vcf_record_ptr,
            genotype_from_coverage_using_maximum_likelihood_path_as_reference)
            .Times(1);
    }

    EXPECT_CALL(*non_snp_vcf_record_ptr, genotype_from_coverage).Times(0);
    EXPECT_CALL(*snp_vcf_record_ptr, genotype_from_coverage).Times(0);
    EXPECT_CALL(default_vcf, make_gt_compatible).Times(0);

    default_vcf.genotype(true, false);
}

TEST_F(VCFTest___genotype___Fixture, genotype_all_records___local_genotyping)
{
    {
        InSequence seq;
        EXPECT_CALL(*non_snp_vcf_record_ptr, genotype_from_coverage).Times(1);
        EXPECT_CALL(*snp_vcf_record_ptr, genotype_from_coverage).Times(1);
        EXPECT_CALL(default_vcf, make_gt_compatible).Times(1);
    }

    EXPECT_CALL(*non_snp_vcf_record_ptr,
        genotype_from_coverage_using_maximum_likelihood_path_as_reference)
        .Times(0);
    EXPECT_CALL(*snp_vcf_record_ptr,
        genotype_from_coverage_using_maximum_likelihood_path_as_reference)
        .Times(0);

    default_vcf.genotype(false, true);
}

TEST_F(VCFTest___genotype___Fixture, genotype_snp_records_only___global_genotyping)
{
    {
        InSequence seq;
        EXPECT_CALL(*snp_vcf_record_ptr,
            genotype_from_coverage_using_maximum_likelihood_path_as_reference)
            .Times(1);
    }

    EXPECT_CALL(*non_snp_vcf_record_ptr,
        genotype_from_coverage_using_maximum_likelihood_path_as_reference)
        .Times(0);
    EXPECT_CALL(*snp_vcf_record_ptr, genotype_from_coverage).Times(0);
    EXPECT_CALL(*non_snp_vcf_record_ptr, genotype_from_coverage).Times(0);
    EXPECT_CALL(snps_only_vcf, make_gt_compatible).Times(0);

    snps_only_vcf.genotype(true, false);
}

TEST_F(VCFTest___genotype___Fixture, genotype_snp_records_only___local_genotyping)
{
    {
        InSequence seq;
        EXPECT_CALL(*snp_vcf_record_ptr, genotype_from_coverage).Times(1);
        EXPECT_CALL(snps_only_vcf, make_gt_compatible).Times(1);
    }

    EXPECT_CALL(*non_snp_vcf_record_ptr, genotype_from_coverage).Times(0);
    EXPECT_CALL(*snp_vcf_record_ptr,
        genotype_from_coverage_using_maximum_likelihood_path_as_reference)
        .Times(0);
    EXPECT_CALL(*non_snp_vcf_record_ptr,
        genotype_from_coverage_using_maximum_likelihood_path_as_reference)
        .Times(0);

    snps_only_vcf.genotype(false, true);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OLD GENOTYPING TEST THAT WILL BE READDED AS INTEGRATION TEST
// INTEGRATION TEST BECAUSE THEY TEST THE WHOLE GENOTYPING PIPELINE
// TODO: READD
// FROM VCF::genotype() to VCF_Record::genotype_from_coverage() to
// SampleIndexToSampleInfo::genotype_from_coverage() to
// SampleInfo::genotype_from_coverage()
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TEST(VCFTest, genotype) {
//    // not a snp site
//    // missing count data
//    // not confident
//    // confident and have right gt
//    // confident and have wrong gt
//    // one sample needs regenotyping and the other doesn't
//    VCF vcf = create_VCF_with_default_parameters();
//
//    vcf.add_record("chrom2", 79, "C", "G");
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 2,
//    "T", "TA"); vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample",
//    "chrom1", 5, "A", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    79, "C", "A");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2",
//    20, "A", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2",
//    79, "C", "C");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2",
//    80, "A", "C");
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1",
//    2, "T", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1",
//    5, "A", "A");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1",
//    79, "C", "A");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2",
//    20, "A", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2",
//    79, "C", "C");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2",
//    80, "A", "A");
//
//    vcf.sort_records();
//    std::vector<float> f = {0.0, 0.0};
//
//    // record 0, not a snp site
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 10};
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 20};
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"] = {1, 15};
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"] = {2, 24};
//    vcf.get_records()[0]->set_format(0,"GAPS", f);
//    vcf.get_records()[0]->set_format(1,"GAPS", f);
//
//
//    // record 1, different genotypes but both correct
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 10};
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 20};
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"] = {10, 1};
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"] = {21, 2};
//    vcf.get_records()[1]->set_format(0,"GAPS", f);
//    vcf.get_records()[1]->set_format(1,"GAPS", f);
//
//    // record 2, same genotypes first correct
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 10};
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 20};
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"] = {10, 1};
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"] = {21, 2};
//    vcf.get_records()[2]->set_format(0,"GAPS", f);
//    vcf.get_records()[2]->set_format(1,"GAPS", f);
//
//    // record 3, same genotypes both wrong
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {20, 1};
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {21, 2};
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"] = {10, 1};
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"] = {21, 2};
//    vcf.get_records()[3]->set_format(0,"GAPS", f);
//    vcf.get_records()[3]->set_format(1,"GAPS", f);
//
//    // record 4, missing count data for first sample
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 10};
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {20};
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"] = {10, 1};
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"] = {21, 2};
//    vcf.get_records()[4]->set_format(0,"GAPS", f);
//    vcf.get_records()[4]->set_format(1,"GAPS", f);
//
//    // record 5, not confident for second sample
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 10};
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 20};
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"] = {2, 1};
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"] = {4, 2};
//    vcf.get_records()[5]->set_format(0,"GAPS", f);
//    vcf.get_records()[5]->set_format(1,"GAPS", f);
//
//    vcf.genotype({30, 30}, 0.01, 30, 0, 1, 0, 0, true);
//
//    // not genotyped first record
//    EXPECT_EQ((uint16_t) 1,
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint16_t)
//    1, vcf.get_records()[0]->sampleIndex_to_sampleInfo[1]["GT"][0]); bool
//    found_confidence =
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].find("GT_CONF") !=
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_confidence);
//    found_confidence =
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1].find("GT_CONF") !=
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1].end();
//    EXPECT_FALSE(found_confidence);
//
//    // both correct
//    EXPECT_EQ((uint) 2, vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());
//    bool valid_gt = vcf.get_records()[1]->sampleIndex_to_sampleInfo[0].find("GT") !=
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0].end(); EXPECT_TRUE(valid_gt);
//    valid_gt = vcf.get_records()[1]->sampleIndex_to_sampleInfo[1].find("GT") !=
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[1].end(); EXPECT_TRUE(valid_gt);
//    EXPECT_EQ((uint) 1,
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["GT"].size()); EXPECT_EQ((uint)
//    1, vcf.get_records()[1]->sampleIndex_to_sampleInfo[1]["GT"].size());
//    EXPECT_EQ((uint16_t) 1,
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint16_t)
//    0, vcf.get_records()[1]->sampleIndex_to_sampleInfo[1]["GT"][0]);
//    // first correct
//    EXPECT_EQ((uint16_t) 1,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint16_t)
//    0, vcf.get_records()[2]->sampleIndex_to_sampleInfo[1]["GT"][0]);
//    // both wrong
//    EXPECT_EQ((uint16_t) 0,
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint16_t)
//    0, vcf.get_records()[3]->sampleIndex_to_sampleInfo[1]["GT"][0]);
//    // first missing data
//    EXPECT_EQ((uint) 0,
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint16_t) 0,
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[1]["GT"][0]);
//    // second not confident
//    EXPECT_EQ((uint16_t) 1,
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint) 0,
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[1]["GT"].size());
//}
//
// TEST(VCFTest, genotype_with_all_sites) {
//    // not a snp site
//    // missing count data
//    // not confident
//    // confident and have right gt
//    // confident and have wrong gt
//    // one sample needs regenotyping and the other doesn't
//    VCF vcf = create_VCF_with_default_parameters();
//
//    vcf.add_record("chrom2", 79, "CC", "GC");
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 2,
//    "T", "TA"); vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample",
//    "chrom1", 5, "AC", "GC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    79, "CC", "AC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2",
//    20, "AC", "GC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2",
//    79, "CC", "CC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2",
//    80, "AC", "CC");
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1",
//    2, "T", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1",
//    5, "AC", "AC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1",
//    79, "CC", "AC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2",
//    20, "AC", "GC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2",
//    79, "CC", "CC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2",
//    80, "AC", "AC");
//
//    vcf.sort_records();
//    std::vector<float> f = {0.0, 0.0};
//
//    // record 0, not a snp site
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(0);
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(1);
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(10);
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(20);
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(15);
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(24);
//    vcf.get_records()[0]->set_format(0,"GAPS", f);
//    vcf.get_records()[0]->set_format(1,"GAPS", f);
//
//    // record 1, different genotypes but both correct
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(0);
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(1);
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(10);
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(20);
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(10);
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(21);
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.get_records()[1]->set_format(0,"GAPS", f);
//    vcf.get_records()[1]->set_format(1,"GAPS", f);
//
//    // record 2, same genotypes first correct
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(0);
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(1);
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(10);
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(20);
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(10);
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(21);
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.get_records()[2]->set_format(0,"GAPS", f);
//    vcf.get_records()[2]->set_format(1,"GAPS", f);
//
//    // record 3, same genotypes both wrong
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(20);
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(21);
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(1);
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(2);
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(10);
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(21);
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.get_records()[3]->set_format(0,"GAPS", f);
//    vcf.get_records()[3]->set_format(1,"GAPS", f);
//
//    // record 4, missing count data for first sample
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(0);
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(10);
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(20);
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(10);
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(21);
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.get_records()[4]->set_format(0,"GAPS", f);
//    vcf.get_records()[4]->set_format(1,"GAPS", f);
//
//    // record 5, not confident for second sample
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(0);
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(1);
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(10);
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(20);
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(2);
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(4);
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.get_records()[5]->set_format(0,"GAPS", f);
//    vcf.get_records()[5]->set_format(1,"GAPS", f);
//
//    bool snps_only = false;
//    vcf.genotype({30, 30}, 0.01, 30, 0, 1, 0, 0, snps_only);
//
//    // first record now genotyped
//    EXPECT_EQ((uint16_t) 1,
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint16_t)
//    1, vcf.get_records()[0]->sampleIndex_to_sampleInfo[1]["GT"][0]); bool
//    found_confidence =
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].find("GT_CONF") !=
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].end();
//    EXPECT_TRUE(found_confidence);
//    found_confidence =
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1].find("GT_CONF") !=
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1].end();
//    EXPECT_TRUE(found_confidence);
//    // both correct
//    EXPECT_EQ((uint) 2, vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());
//    bool valid_gt = vcf.get_records()[1]->sampleIndex_to_sampleInfo[0].find("GT") !=
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0].end(); EXPECT_TRUE(valid_gt);
//    valid_gt = vcf.get_records()[1]->sampleIndex_to_sampleInfo[1].find("GT") !=
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[1].end(); EXPECT_TRUE(valid_gt);
//    EXPECT_EQ((uint) 1,
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["GT"].size()); EXPECT_EQ((uint)
//    1, vcf.get_records()[1]->sampleIndex_to_sampleInfo[1]["GT"].size());
//    EXPECT_EQ((uint16_t) 1,
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint16_t)
//    0, vcf.get_records()[1]->sampleIndex_to_sampleInfo[1]["GT"][0]);
//    // first correct
//    EXPECT_EQ((uint16_t) 1,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint16_t)
//    0, vcf.get_records()[2]->sampleIndex_to_sampleInfo[1]["GT"][0]);
//    // both wrong
//    EXPECT_EQ((uint16_t) 0,
//    vcf.get_records()[3]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint16_t)
//    0, vcf.get_records()[3]->sampleIndex_to_sampleInfo[1]["GT"][0]);
//    // first missing data
//    EXPECT_EQ((uint) 0,
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint16_t) 0,
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[1]["GT"][0]);
//    // second not confident
//    EXPECT_EQ((uint16_t) 1,
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint) 0,
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[1]["GT"].size());
//
//}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// REASON COMMENTED OUT: THIS METHOD DOES NOT EXIST ANYMORE
// TEST(VCFTest, clean) {
//    VCF vcf = create_VCF_with_default_parameters();
//
//    VCFRecord dummy;
//    std::vector<std::string> empty = {};
//    vcf.add_record(dummy, empty);
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 2,
//    "T", "TA"); vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample",
//    "chrom1", 5, "A", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    79, "C", "A"); vcf.get_records()[2]->clear(); EXPECT_EQ((uint) 5,
//    vcf.get_VCF_size());
//
//    vcf.clean();
//    EXPECT_EQ((uint) 3, vcf.get_VCF_size());
//    EXPECT_EQ((uint) 79, vcf.get_records()[0]->get_pos());
//    EXPECT_EQ((uint) 1, vcf.get_records()[0]->get_alts().size());
//    EXPECT_EQ("G", vcf.get_records()[0]->get_alts()[0]);
//    EXPECT_EQ((uint) 5, vcf.get_records()[1]->get_pos());
//    EXPECT_EQ((uint) 79, vcf.get_records()[2]->get_pos());
//    EXPECT_EQ((uint) 1, vcf.get_records()[2]->get_alts().size());
//    EXPECT_EQ("A", vcf.get_records()[2]->get_alts()[0]);
//}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VCFTest___merge_multi_allelic_core___Fixture : public ::testing::Test {
public:
    class VCFVisibilityMock : public VCF {
    public:
        using VCF::merge_multi_allelic_core;
        using VCF::VCF;
    };

    class VCFMock : public VCFVisibilityMock {
    public:
        using VCF::records;
        using VCFVisibilityMock::VCFVisibilityMock;
        MOCK_METHOD(void, add_samples, (const std::vector<std::string>& sample_names),
            (override));
        MOCK_METHOD(void, add_record_core, (const VCFRecord& vr), (override));
        MOCK_METHOD(void, sort_records, (), (override));
    };

    class VCFRecordMock : public VCFRecord {
    public:
        using VCFRecord::VCFRecord;
        MOCK_METHOD(bool, can_biallelic_record_be_merged_into_this,
            (const VCFRecord& vcf_record_to_be_merged_in, uint32_t max_allele_length),
            (const override));
        MOCK_METHOD(void, merge_record_into_this, (const VCFRecord& other), (override));
        MOCK_METHOD(
            std::shared_ptr<VCFRecord>, make_copy_as_shared_ptr, (), (const override));
    };

    VCFTest___merge_multi_allelic_core___Fixture()
        : vcf(&default_genotyping_options)
        , merged_vcf(&default_genotyping_options)
    {
    }

    void SetUp() override {}

    void TearDown() override {}

    VCFMock vcf;
    VCFMock merged_vcf;
};

TEST_F(VCFTest___merge_multi_allelic_core___Fixture, merged_VCF_is_not_initially_empty)
{
    VCFVisibilityMock vcf(&default_genotyping_options);
    VCF merged_vcf(vcf.genotyping_options);
    merged_vcf.add_record("1", 1, "A", "T");

    EXPECT_DEATH(vcf.merge_multi_allelic_core(merged_vcf, 10000), "");
}

TEST_F(VCFTest___merge_multi_allelic_core___Fixture, one_sized_VCF)
{
    VCFVisibilityMock one_sized_vcf(&default_genotyping_options);
    one_sized_vcf.add_record("1", 1, "A", "T");

    VCFVisibilityMock merged_vcf(one_sized_vcf.genotyping_options);
    one_sized_vcf.merge_multi_allelic_core(merged_vcf, 10000);

    EXPECT_EQ(one_sized_vcf, merged_vcf);
}

TEST_F(VCFTest___merge_multi_allelic_core___Fixture, three_records_to_be_merged)
{
    // expectations for the merged record
    auto merged_record_1
        = std::make_shared<VCFRecordMock>(&merged_vcf, "merged_1", 1, "A", "T");
    ;
    EXPECT_CALL(*merged_record_1,
        can_biallelic_record_be_merged_into_this(
            Property(&VCFRecord::get_chrom, "original_2"), _))
        .Times(1)
        .WillOnce(Return(true));
    EXPECT_CALL(*merged_record_1,
        merge_record_into_this(Property(&VCFRecord::get_chrom, "original_2")))
        .Times(1);
    EXPECT_CALL(*merged_record_1,
        can_biallelic_record_be_merged_into_this(
            Property(&VCFRecord::get_chrom, "original_3"), _))
        .Times(1)
        .WillOnce(Return(true));
    EXPECT_CALL(*merged_record_1,
        merge_record_into_this(Property(&VCFRecord::get_chrom, "original_3")))
        .Times(1);

    // expectations for the original records
    auto vcf_record_1
        = std::make_shared<VCFRecordMock>(&vcf, "original_1", 1, "A", "T");
    EXPECT_CALL(*vcf_record_1, make_copy_as_shared_ptr)
        .Times(1)
        .WillOnce(Return(merged_record_1));
    auto vcf_record_2
        = std::make_shared<VCFRecordMock>(&vcf, "original_2", 2, "A", "T");
    EXPECT_CALL(*vcf_record_2, make_copy_as_shared_ptr).Times(0);
    auto vcf_record_3
        = std::make_shared<VCFRecordMock>(&vcf, "original_3", 3, "A", "T");
    EXPECT_CALL(*vcf_record_3, make_copy_as_shared_ptr).Times(0);

    // expectations for the merged vcf
    EXPECT_CALL(merged_vcf, add_samples(vcf.samples)).Times(1);
    EXPECT_CALL(
        merged_vcf, add_record_core(Property(&VCFRecord::get_chrom, "merged_1")))
        .Times(1);
    EXPECT_CALL(merged_vcf, sort_records).Times(1);

    vcf.records.push_back(vcf_record_1);
    vcf.records.push_back(vcf_record_2);
    vcf.records.push_back(vcf_record_3);
    vcf.merge_multi_allelic_core(merged_vcf, 10000);
}

TEST_F(VCFTest___merge_multi_allelic_core___Fixture, three_non_mergeable_records)
{
    // expectations for the merged records
    auto merged_record_1
        = std::make_shared<VCFRecordMock>(&merged_vcf, "merged_1", 1, "A", "T");
    ;
    EXPECT_CALL(*merged_record_1,
        can_biallelic_record_be_merged_into_this(
            Property(&VCFRecord::get_chrom, "original_2"), _))
        .Times(1)
        .WillOnce(Return(false));
    EXPECT_CALL(*merged_record_1, merge_record_into_this).Times(0);
    auto merged_record_2
        = std::make_shared<VCFRecordMock>(&merged_vcf, "merged_2", 1, "A", "T");
    ;
    EXPECT_CALL(*merged_record_2,
        can_biallelic_record_be_merged_into_this(
            Property(&VCFRecord::get_chrom, "original_3"), _))
        .Times(1)
        .WillOnce(Return(false));
    EXPECT_CALL(*merged_record_2, merge_record_into_this).Times(0);
    auto merged_record_3
        = std::make_shared<VCFRecordMock>(&merged_vcf, "merged_3", 1, "A", "T");
    ;
    EXPECT_CALL(*merged_record_3, can_biallelic_record_be_merged_into_this).Times(0);
    EXPECT_CALL(*merged_record_3, merge_record_into_this).Times(0);

    // expectations for the original records
    auto vcf_record_1
        = std::make_shared<VCFRecordMock>(&vcf, "original_1", 1, "A", "T");
    EXPECT_CALL(*vcf_record_1, make_copy_as_shared_ptr)
        .Times(1)
        .WillOnce(Return(merged_record_1));
    auto vcf_record_2
        = std::make_shared<VCFRecordMock>(&vcf, "original_2", 2, "A", "T");
    EXPECT_CALL(*vcf_record_2, make_copy_as_shared_ptr)
        .Times(1)
        .WillOnce(Return(merged_record_2));
    auto vcf_record_3
        = std::make_shared<VCFRecordMock>(&vcf, "original_3", 3, "A", "T");
    EXPECT_CALL(*vcf_record_3, make_copy_as_shared_ptr)
        .Times(1)
        .WillOnce(Return(merged_record_3));

    // expectations for the merged vcf
    EXPECT_CALL(merged_vcf, add_samples(vcf.samples)).Times(1);
    EXPECT_CALL(
        merged_vcf, add_record_core(Property(&VCFRecord::get_chrom, "merged_1")))
        .Times(1);
    EXPECT_CALL(
        merged_vcf, add_record_core(Property(&VCFRecord::get_chrom, "merged_2")))
        .Times(1);
    EXPECT_CALL(
        merged_vcf, add_record_core(Property(&VCFRecord::get_chrom, "merged_3")))
        .Times(1);
    EXPECT_CALL(merged_vcf, sort_records).Times(1);

    vcf.records.push_back(vcf_record_1);
    vcf.records.push_back(vcf_record_2);
    vcf.records.push_back(vcf_record_3);
    vcf.merge_multi_allelic_core(merged_vcf, 10000);
}

TEST_F(VCFTest___merge_multi_allelic_core___Fixture,
    two_records_to_be_merged___followed_by_one_record_not_mergeable___followed_by_two_to_be_merged)
{
    // expectations for the merged records
    auto merged_record_1
        = std::make_shared<VCFRecordMock>(&merged_vcf, "merged_1", 1, "A", "T");
    ;
    EXPECT_CALL(*merged_record_1,
        can_biallelic_record_be_merged_into_this(
            Property(&VCFRecord::get_chrom, "original_2"), _))
        .Times(1)
        .WillOnce(Return(true));
    EXPECT_CALL(*merged_record_1,
        merge_record_into_this(Property(&VCFRecord::get_chrom, "original_2")))
        .Times(1);
    EXPECT_CALL(*merged_record_1,
        can_biallelic_record_be_merged_into_this(
            Property(&VCFRecord::get_chrom, "original_3"), _))
        .Times(1)
        .WillOnce(Return(false));
    auto merged_record_2
        = std::make_shared<VCFRecordMock>(&merged_vcf, "merged_2", 1, "A", "T");
    ;
    EXPECT_CALL(*merged_record_2,
        can_biallelic_record_be_merged_into_this(
            Property(&VCFRecord::get_chrom, "original_4"), _))
        .Times(1)
        .WillOnce(Return(false));
    EXPECT_CALL(*merged_record_2, merge_record_into_this).Times(0);
    auto merged_record_3
        = std::make_shared<VCFRecordMock>(&merged_vcf, "merged_3", 1, "A", "T");
    ;
    EXPECT_CALL(*merged_record_3,
        can_biallelic_record_be_merged_into_this(
            Property(&VCFRecord::get_chrom, "original_5"), _))
        .Times(1)
        .WillOnce(Return(true));
    EXPECT_CALL(*merged_record_3,
        merge_record_into_this(Property(&VCFRecord::get_chrom, "original_5")))
        .Times(1);

    // expectations for the original records
    auto vcf_record_1
        = std::make_shared<VCFRecordMock>(&vcf, "original_1", 1, "A", "T");
    EXPECT_CALL(*vcf_record_1, make_copy_as_shared_ptr)
        .Times(1)
        .WillOnce(Return(merged_record_1));
    auto vcf_record_2
        = std::make_shared<VCFRecordMock>(&vcf, "original_2", 2, "A", "T");
    EXPECT_CALL(*vcf_record_2, make_copy_as_shared_ptr).Times(0);
    auto vcf_record_3
        = std::make_shared<VCFRecordMock>(&vcf, "original_3", 3, "A", "T");
    EXPECT_CALL(*vcf_record_3, make_copy_as_shared_ptr)
        .Times(1)
        .WillOnce(Return(merged_record_2));
    auto vcf_record_4
        = std::make_shared<VCFRecordMock>(&vcf, "original_4", 4, "A", "T");
    EXPECT_CALL(*vcf_record_4, make_copy_as_shared_ptr)
        .Times(1)
        .WillOnce(Return(merged_record_3));
    auto vcf_record_5
        = std::make_shared<VCFRecordMock>(&vcf, "original_5", 5, "A", "T");
    EXPECT_CALL(*vcf_record_5, make_copy_as_shared_ptr).Times(0);

    // expectations for the merged vcf
    EXPECT_CALL(merged_vcf, add_samples(vcf.samples)).Times(1);
    EXPECT_CALL(
        merged_vcf, add_record_core(Property(&VCFRecord::get_chrom, "merged_1")))
        .Times(1);
    EXPECT_CALL(
        merged_vcf, add_record_core(Property(&VCFRecord::get_chrom, "merged_2")))
        .Times(1);
    EXPECT_CALL(
        merged_vcf, add_record_core(Property(&VCFRecord::get_chrom, "merged_3")))
        .Times(1);
    EXPECT_CALL(merged_vcf, sort_records).Times(1);

    vcf.records.push_back(vcf_record_1);
    vcf.records.push_back(vcf_record_2);
    vcf.records.push_back(vcf_record_3);
    vcf.records.push_back(vcf_record_4);
    vcf.records.push_back(vcf_record_5);
    vcf.merge_multi_allelic_core(merged_vcf, 10000);
}

TEST(VCFTest___merge_multi_allelic,
    vcf_with_two_samples_and_two_records_second_record_does_not_map_to_second_sample)
{
    // declares vcf with two samples
    VCF vcf = create_VCF_with_default_parameters(0);
    vcf.add_samples({ "sample1", "sample2" });

    // add two bi-allelic records to be merged
    vcf.add_record("chrom1", 5, "A", "C", "SVTYPE=SNP", "GRAPHTYPE=SIMPLE");
    vcf.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP", "GRAPHTYPE=SIMPLE");

    // causing conflict and using genotype from the coverages to solve
    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].set_gt_from_max_likelihood_path(
        1);
    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1].set_gt_from_max_likelihood_path(
        1);
    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0].set_gt_from_max_likelihood_path(
        1);
    vcf.get_records()[1]->sampleIndex_to_sampleInfo[1].set_gt_from_max_likelihood_path(
        1);

    // first record maps to both samples
    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].set_coverage_information(
        { { 1 }, { 2 } }, { { 0 }, { 0 } });
    vcf.get_records()[0]->sampleIndex_to_sampleInfo[1].set_coverage_information(
        { { 1 }, { 3 } }, { { 0 }, { 0 } });

    // second record maps to first sample only
    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0].set_coverage_information(
        { { 1 }, { 4 } }, { { 0 }, { 0 } });

    // do the merge
    vcf = vcf.merge_multi_allelic();

    EXPECT_EQ(1, vcf.get_VCF_size());
    EXPECT_EQ("chrom1", vcf.get_records()[0]->get_chrom());
    EXPECT_EQ(5, vcf.get_records()[0]->get_pos());
    EXPECT_EQ("A", vcf.get_records()[0]->get_ref());
    EXPECT_EQ(std::vector<std::string>({ "C", "G" }), vcf.get_records()[0]->get_alts());
    EXPECT_EQ(2, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
    EXPECT_EQ(1,
        vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].get_mean_forward_coverage(
            0));
    EXPECT_EQ(2,
        vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].get_mean_forward_coverage(
            1));
    EXPECT_EQ(4,
        vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].get_mean_forward_coverage(
            2));
    EXPECT_EQ(1,
        vcf.get_records()[0]->sampleIndex_to_sampleInfo[1].get_mean_forward_coverage(
            0));
    EXPECT_EQ(3,
        vcf.get_records()[0]->sampleIndex_to_sampleInfo[1].get_mean_forward_coverage(
            1));
    EXPECT_EQ(0,
        vcf.get_records()[0]->sampleIndex_to_sampleInfo[1].get_mean_forward_coverage(
            2));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OLD merge_multi_allelic TESTS - TO BE READDED AS INTEGRATION TESTS
// REASON COMMENTED OUT: These tests are covered by and redundant with
// sampleinfo_test.cpp::merge_other_sample_info_into_this tests
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TEST(VCFTest, merge_multi_allelic) {
//    VCF vcf = create_VCF_with_default_parameters();
//    // no gt
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 5, "A", "C");
//    // gt
//    vcf.add_record("chrom1", 46, "CTT", "A");
//    vcf.add_record("chrom1", 46, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    46, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    46, "CTT", "A");
//    // likelihoods too
//    vcf.add_record("chrom1", 76, "CTT", "A");
//    vcf.add_record("chrom1", 76, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    76, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    76, "CTT", "A");
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1,
//    &default_genotyping_options);
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1,
//    &default_genotyping_options);
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"] = {-50, -3};
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"] = {-50, -16};
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["GT_CONF"] = {47};
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["GT_CONF"] = {56};
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {2, 30};
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {2, 30};
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {2, 30};
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {2, 30};
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["GAPS"] = {4, 0};
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["GAPS"] = {4, 1};
//    // incompatible
//    vcf.add_record("chrom1", 85, "A", "G");
//    vcf.add_record("chrom1", 85, "T", "C");
//
//    vcf = vcf.merge_multi_allelic();
//
//    EXPECT_EQ((uint) 5, vcf.get_VCF_size());
//    EXPECT_EQ((uint) 5, vcf.get_records()[0]->get_pos());
//    EXPECT_EQ((uint) 2, vcf.get_records()[0]->get_alts().size());
//    EXPECT_EQ((uint) 1, vcf.get_records()[0]->sampleIndex_to_sampleInfo.size());
//    EXPECT_EQ((uint) 0, vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].size());
//
//    EXPECT_EQ((uint) 46, vcf.get_records()[1]->get_pos());
//    EXPECT_EQ((uint) 2, vcf.get_records()[1]->get_alts().size());
//    EXPECT_EQ((uint) 1, vcf.get_records()[1]->sampleIndex_to_sampleInfo.size());
//    bool valid_gt = vcf.get_records()[1]->sampleIndex_to_sampleInfo[0].find("GT") !=
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0].end(); EXPECT_TRUE(valid_gt);
//    EXPECT_EQ((uint) 0,
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0]["GT"].size());
//
//    EXPECT_EQ((uint) 76, vcf.get_records()[2]->get_pos());
//    EXPECT_EQ((uint) 2, vcf.get_records()[2]->get_alts().size());
//    EXPECT_EQ((uint) 1, vcf.get_records()[2]->sampleIndex_to_sampleInfo.size());
//    valid_gt = vcf.get_records()[2]->sampleIndex_to_sampleInfo[0].find("GT") !=
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0].end(); EXPECT_TRUE(valid_gt);
//    EXPECT_EQ((uint) 1, vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint) 3, vcf.get_records()[2]->sampleIndex_to_sampleInfo[0].size());
//    bool found = vcf.get_records()[2]->sampleIndex_to_sampleInfo[0].find("LIKELIHOOD")
//    != vcf.get_records()[2]->sampleIndex_to_sampleInfo[0].end(); EXPECT_TRUE(found);
//    EXPECT_EQ((uint) 3,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"].size());
//    EXPECT_EQ(-50.0,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"][0]);
//    EXPECT_EQ(-3.0,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"][1]);
//    EXPECT_EQ(-16.0,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"][2]);
//    EXPECT_EQ((uint) 3,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["GAPS"].size()); EXPECT_EQ(4,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["GAPS"][0]); EXPECT_EQ(0,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["GAPS"][1]); EXPECT_EQ(1,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["GAPS"][2]); found =
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0].find("GT_CONF") !=
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0].end(); EXPECT_TRUE(found);
//    EXPECT_EQ((uint) 1,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["GT_CONF"].size());
//    EXPECT_EQ(13.0, vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["GT_CONF"][0]);
//
//    EXPECT_EQ((uint) 85, vcf.get_records()[3]->get_pos());
//    EXPECT_EQ((uint) 1, vcf.get_records()[3]->get_alts().size());
//    EXPECT_EQ((uint) 85, vcf.get_records()[4]->get_pos());
//    EXPECT_EQ((uint) 1, vcf.get_records()[4]->get_alts().size());
//}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST(VCFTest, correct_dot_alleles)
{
    VCF vcf = create_VCF_with_default_parameters(0);
    // at start
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 0, ".", "TA");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom2", 0, "T", ".");
    // in middle
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 35, ".", "A");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom2", 35, "TA", ".");
    // multiple alts
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 44, "TA", "T");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom1", 44, "TA", ".");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom2", 44, ".", "T");
    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it(
        "sample", "chrom2", 44, ".", "TA");

    string vcf_ref = "TATATGTGTC"
                     "GCGACACTGC"
                     "ATGCATGCAT"
                     "AGTCCTAAAG"
                     "TCCTTAAACG"
                     "TTTATAGTCG";

    vcf = vcf.correct_dot_alleles(vcf_ref, "chrom1");
    vcf = vcf.correct_dot_alleles(vcf_ref, "chrom2");

    EXPECT_EQ(vcf.get_records()[0]->get_ref(), "T");
    EXPECT_EQ(vcf.get_records()[1]->get_ref(), "C");
    EXPECT_EQ(vcf.get_records()[2]->get_ref(), "TTA");
    EXPECT_EQ(vcf.get_records()[3]->get_ref(), "TA");
    EXPECT_EQ(vcf.get_records()[4]->get_ref(), "TA");
    EXPECT_EQ(vcf.get_records()[5]->get_ref(), "CTA");
    EXPECT_EQ(vcf.get_records()[6]->get_ref(), "T");
    EXPECT_EQ(vcf.get_records()[7]->get_ref(), "T");

    EXPECT_EQ(vcf.get_records()[0]->get_alts().size(), 1);
    EXPECT_EQ(vcf.get_records()[0]->get_alts()[0], "TAT");
    EXPECT_EQ(vcf.get_records()[1]->get_alts().size(), 1);
    EXPECT_EQ(vcf.get_records()[1]->get_alts()[0], "CA");
    EXPECT_EQ(vcf.get_records()[2]->get_alts().size(), 1);
    EXPECT_EQ(vcf.get_records()[2]->get_alts()[0], "T");
    EXPECT_EQ(vcf.get_records()[3]->get_alts().size(), 1);
    EXPECT_EQ(vcf.get_records()[3]->get_alts()[0], "T");
    EXPECT_EQ(vcf.get_records()[4]->get_alts().size(), 1);
    EXPECT_EQ(vcf.get_records()[4]->get_alts()[0], "A");
    EXPECT_EQ(vcf.get_records()[5]->get_alts().size(), 1);
    EXPECT_EQ(vcf.get_records()[5]->get_alts()[0], "C");
    EXPECT_EQ(vcf.get_records()[6]->get_alts().size(), 1);
    EXPECT_EQ(vcf.get_records()[6]->get_alts()[0], "TT");
    EXPECT_EQ(vcf.get_records()[7]->get_alts().size(), 1);
    EXPECT_EQ(vcf.get_records()[7]->get_alts()[0], "TTA");
}

class VCFTest___make_gt_compatible___Fixture : public ::testing::Test {
public:
    class VCFRecordMock : public VCFRecord {
    public:
        using VCFRecord::VCFRecord;
        MOCK_METHOD(
            void, solve_incompatible_gt_conflict_with, (VCFRecord & other), (override));
    };

    class VCFMock : public VCF {
    public:
        using VCF::records;
        using VCF::VCF;
        MOCK_METHOD(std::vector<VCFRecord*>,
            get_all_records_overlapping_the_given_record, (const VCFRecord& vcf_record),
            (override));
    };

    VCFTest___make_gt_compatible___Fixture()
        : vcf(&default_genotyping_options)
        , vcf_record_1(std::make_shared<VCFRecordMock>(&vcf, "1", 1, "A", "T"))
        , vcf_record_2(std::make_shared<VCFRecordMock>(&vcf, "1", 2, "A", "T"))
        , vcf_record_3(std::make_shared<VCFRecordMock>(&vcf, "1", 3, "A", "T"))
        , vcf_record_3B(std::make_shared<VCFRecordMock>(&vcf, "1", 3, "AC", "T"))
    {
    }

    void SetUp() override
    {
        vcf.records.push_back(vcf_record_1);
        vcf.records.push_back(vcf_record_2);
    }

    void TearDown() override {}

    VCFMock vcf;
    std::shared_ptr<VCFRecordMock> vcf_record_1;
    std::shared_ptr<VCFRecordMock> vcf_record_2;
    std::shared_ptr<VCFRecordMock> vcf_record_3;
    std::shared_ptr<VCFRecordMock> vcf_record_3B;
};

TEST_F(VCFTest___make_gt_compatible___Fixture,
    side_by_side_no_overlapping_records___no_conflict)
{
    {
        InSequence seq;
        EXPECT_CALL(vcf,
            get_all_records_overlapping_the_given_record(
                Property(&VCFRecord::get_pos, 1)))
            .Times(1)
            .WillOnce(Return(std::vector<VCFRecord*>({ vcf_record_1.get() })));
        EXPECT_CALL(vcf,
            get_all_records_overlapping_the_given_record(
                Property(&VCFRecord::get_pos, 2)))
            .Times(1)
            .WillOnce(Return(std::vector<VCFRecord*>({ vcf_record_2.get() })));
    }

    EXPECT_CALL(*vcf_record_1, solve_incompatible_gt_conflict_with).Times(0);
    EXPECT_CALL(*vcf_record_2, solve_incompatible_gt_conflict_with).Times(0);

    vcf.make_gt_compatible();
}

TEST_F(VCFTest___make_gt_compatible___Fixture, two_overlapping_records___conflict)
{
    {
        InSequence seq;

        EXPECT_CALL(vcf,
            get_all_records_overlapping_the_given_record(
                Property(&VCFRecord::get_pos, 1)))
            .Times(1)
            .WillOnce(Return(
                std::vector<VCFRecord*>({ vcf_record_1.get(), vcf_record_2.get() })));

        EXPECT_CALL(*vcf_record_1,
            solve_incompatible_gt_conflict_with(Property(&VCFRecord::get_pos, 2)))
            .Times(1);

        EXPECT_CALL(vcf,
            get_all_records_overlapping_the_given_record(
                Property(&VCFRecord::get_pos, 2)))
            .Times(1)
            .WillOnce(Return(
                std::vector<VCFRecord*>({ vcf_record_1.get(), vcf_record_2.get() })));
    }

    EXPECT_CALL(*vcf_record_1,
        solve_incompatible_gt_conflict_with(Property(&VCFRecord::get_pos, 1)))
        .Times(0);
    EXPECT_CALL(*vcf_record_2, solve_incompatible_gt_conflict_with).Times(0);

    vcf.make_gt_compatible();
}

TEST_F(VCFTest___make_gt_compatible___Fixture,
    two_overlapping_records_starting_on_the_same_pos___conflict)
{
    VCFMock vcf(&default_genotyping_options);
    vcf.records.push_back(vcf_record_3);
    vcf.records.push_back(vcf_record_3B);

    {
        InSequence seq;

        EXPECT_CALL(vcf,
            get_all_records_overlapping_the_given_record(
                Property(&VCFRecord::get_ref, "A")))
            .Times(1)
            .WillOnce(Return(
                std::vector<VCFRecord*>({ vcf_record_3.get(), vcf_record_3B.get() })));

        EXPECT_CALL(*vcf_record_3,
            solve_incompatible_gt_conflict_with(Property(&VCFRecord::get_ref, "AC")))
            .Times(1);

        EXPECT_CALL(vcf,
            get_all_records_overlapping_the_given_record(
                Property(&VCFRecord::get_ref, "AC")))
            .Times(1)
            .WillOnce(Return(std::vector<VCFRecord*>({ vcf_record_3B.get() })));
    }

    EXPECT_CALL(*vcf_record_3B,
        solve_incompatible_gt_conflict_with(Property(&VCFRecord::get_pos, 1)))
        .Times(0);

    vcf.make_gt_compatible();
}

TEST_F(VCFTest___make_gt_compatible___Fixture, several_records___conflict)
{
    vcf.records.push_back(vcf_record_3);
    {
        InSequence seq;

        EXPECT_CALL(vcf,
            get_all_records_overlapping_the_given_record(
                Property(&VCFRecord::get_pos, 1)))
            .Times(1)
            .WillOnce(Return(std::vector<VCFRecord*>(
                { vcf_record_1.get(), vcf_record_2.get(), vcf_record_3.get() })));

        EXPECT_CALL(*vcf_record_1,
            solve_incompatible_gt_conflict_with(Property(&VCFRecord::get_pos, 2)))
            .Times(1);

        EXPECT_CALL(*vcf_record_1,
            solve_incompatible_gt_conflict_with(Property(&VCFRecord::get_pos, 3)))
            .Times(1);

        EXPECT_CALL(vcf,
            get_all_records_overlapping_the_given_record(
                Property(&VCFRecord::get_pos, 2)))
            .Times(1)
            .WillOnce(Return(std::vector<VCFRecord*>(
                { vcf_record_1.get(), vcf_record_2.get(), vcf_record_3.get() })));

        EXPECT_CALL(*vcf_record_2,
            solve_incompatible_gt_conflict_with(Property(&VCFRecord::get_pos, 3)))
            .Times(1);

        EXPECT_CALL(vcf,
            get_all_records_overlapping_the_given_record(
                Property(&VCFRecord::get_pos, 3)))
            .Times(1)
            .WillOnce(Return(std::vector<VCFRecord*>(
                { vcf_record_1.get(), vcf_record_2.get(), vcf_record_3.get() })));
    }

    EXPECT_CALL(*vcf_record_1,
        solve_incompatible_gt_conflict_with(Property(&VCFRecord::get_pos, 1)))
        .Times(0);
    EXPECT_CALL(*vcf_record_2,
        solve_incompatible_gt_conflict_with(Property(&VCFRecord::get_pos, 1)))
        .Times(0);
    EXPECT_CALL(*vcf_record_2,
        solve_incompatible_gt_conflict_with(Property(&VCFRecord::get_pos, 2)))
        .Times(0);
    EXPECT_CALL(*vcf_record_3, solve_incompatible_gt_conflict_with).Times(0);

    vcf.make_gt_compatible();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OLD make_gt_compatible TEST THAT WILL BE READDED AS INTEGRATION TEST
// TODO: READD
// TODO: SPLIT THIS TEST INTO 5 OR 6
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TEST(VCFTest, make_gt_compatible) {
//    VCF vcf = create_VCF_with_default_parameters();
//    // no gt
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 5, "A", "C");
//    // gt incompatible no likelihoods
//    vcf.add_record("chrom1", 46, "CTT", "A");
//    vcf.add_record("chrom1", 46, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    46, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    46, "CTT", "A");
//    // gt incompatible, likelihoods too both alts
//    vcf.add_record("chrom1", 76, "CTT", "A");
//    vcf.add_record("chrom1", 76, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    76, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    76, "CTT", "A");
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1,
//    vcf.genotyping_options);
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1,
//    vcf.genotyping_options);
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"] = {-50, -3};
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"] = {-50, -16};
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["GT_CONF"] = {47};
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["GT_CONF"] = {56};
//    // gt incompatible one ref, ref correct
//    vcf.add_record("chrom1", 85, "A", "G");
//    vcf.add_record("chrom1", 85, "A", "C");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    85, "A", "A"); vcf.get_records()[6]->sampleIndex_to_sampleInfo[0]["GT"] = {1};
//    vcf.get_records()[6]->sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1);
//    vcf.get_records()[7]->sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1);
//    vcf.get_records()[6]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"] = {-5, -30};
//    vcf.get_records()[7]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"] = {-5, -16};
//    vcf.get_records()[6]->sampleIndex_to_sampleInfo[0]["GT_CONF"] = {47};
//    vcf.get_records()[7]->sampleIndex_to_sampleInfo[0]["GT_CONF"] = {56};
//    // gt incompatible one ref, ref wrong
//    vcf.add_record("chrom1", 95, "A", "G");
//    vcf.add_record("chrom1", 95, "A", "C");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    95, "A", "A"); vcf.get_records()[8]->sampleIndex_to_sampleInfo[0]["GT"] = {1};
//    vcf.get_records()[8]->sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1);
//    vcf.get_records()[9]->sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1);
//    vcf.get_records()[8]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"] = {-50, -3};
//    vcf.get_records()[9]->sampleIndex_to_sampleInfo[0]["LIKELIHOOD"] = {-50, -60};
//    vcf.get_records()[8]->sampleIndex_to_sampleInfo[0]["GT_CONF"] = {47};
//    vcf.get_records()[9]->sampleIndex_to_sampleInfo[0]["GT_CONF"] = {10};
//
//    vcf.make_gt_compatible();
//
//    bool valid_gt = vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].find("GT") !=
//    vcf.get_records()[0]->sampleIndex_to_sampleInfo[0].end(); EXPECT_FALSE(valid_gt);
//    valid_gt = vcf.get_records()[1]->sampleIndex_to_sampleInfo[0].find("GT") !=
//    vcf.get_records()[1]->sampleIndex_to_sampleInfo[0].end(); EXPECT_FALSE(valid_gt);
//    EXPECT_EQ((uint) 0,
//    vcf.get_records()[2]->sampleIndex_to_sampleInfo[0]["GT"].size()); EXPECT_EQ((uint)
//    0, vcf.get_records()[3]->sampleIndex_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint16_t) 1,
//    vcf.get_records()[4]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint) 0,
//    vcf.get_records()[5]->sampleIndex_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint16_t) 0,
//    vcf.get_records()[6]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint16_t)
//    0, vcf.get_records()[7]->sampleIndex_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 1,
//    vcf.get_records()[8]->sampleIndex_to_sampleInfo[0]["GT"][0]); EXPECT_EQ((uint16_t)
//    0, vcf.get_records()[9]->sampleIndex_to_sampleInfo[0]["GT"].size());
//}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VCFTest___get_all_records_overlapping_the_given_record___Fixture
    : public ::testing::Test {
public:
    VCFTest___get_all_records_overlapping_the_given_record___Fixture()
        : vcf(&default_genotyping_options)
        , vcf_record_5_to_10(&vcf, "1", 5, "AAAAA", "T")
        , vcf_record_7_to_8(&vcf, "1", 7, "A", "T")
        , vcf_record_9_to_12(&vcf, "1", 9, "AAA", "T")
        , vcf_record_10_to_12(&vcf, "1", 10, "AA", "T")
        , vcf_record_3_to_6(&vcf, "1", 3, "AAA", "T")
        , vcf_record_3_to_5(&vcf, "1", 3, "AA", "T")
        , vcf_record_3_to_13(&vcf, "1", 3, "AAAAAAAAAA", "T")
        , vcf_record_5_to_10_other_chrom(&vcf, "2", 5, "AAAAA", "T")
    {
    }

    void SetUp() override {}

    void TearDown() override {}

    VCF vcf;
    VCFRecord vcf_record_5_to_10;
    VCFRecord vcf_record_7_to_8;
    VCFRecord vcf_record_9_to_12;
    VCFRecord vcf_record_10_to_12;
    VCFRecord vcf_record_3_to_6;
    VCFRecord vcf_record_3_to_5;
    VCFRecord vcf_record_3_to_13;
    VCFRecord vcf_record_5_to_10_other_chrom;
};

TEST_F(VCFTest___get_all_records_overlapping_the_given_record___Fixture,
    querying_the_indexed_record___returns_itself)
{
    vcf.add_record(vcf_record_5_to_10);

    std::vector<VCFRecord*> expected
        = vcf.get_all_records_overlapping_the_given_record(vcf_record_5_to_10);

    EXPECT_EQ(expected.size(), 1);
    EXPECT_EQ(*(expected[0]), vcf_record_5_to_10);
}

TEST_F(VCFTest___get_all_records_overlapping_the_given_record___Fixture,
    querying_a_record_inside_indexed_record)
{
    vcf.add_record(vcf_record_5_to_10);

    std::vector<VCFRecord*> expected
        = vcf.get_all_records_overlapping_the_given_record(vcf_record_7_to_8);

    EXPECT_EQ(expected.size(), 1);
    EXPECT_EQ(*(expected[0]), vcf_record_5_to_10);
}

TEST_F(VCFTest___get_all_records_overlapping_the_given_record___Fixture,
    querying_a_record_overlapping_just_on_the_right_border)
{
    vcf.add_record(vcf_record_5_to_10);

    std::vector<VCFRecord*> expected
        = vcf.get_all_records_overlapping_the_given_record(vcf_record_9_to_12);

    EXPECT_EQ(expected.size(), 1);
    EXPECT_EQ(*(expected[0]), vcf_record_5_to_10);
}

TEST_F(VCFTest___get_all_records_overlapping_the_given_record___Fixture,
    querying_a_record_overlapping_just_after_the_right_border)
{
    vcf.add_record(vcf_record_5_to_10);

    std::vector<VCFRecord*> expected
        = vcf.get_all_records_overlapping_the_given_record(vcf_record_10_to_12);

    EXPECT_EQ(expected.size(), 0);
}

TEST_F(VCFTest___get_all_records_overlapping_the_given_record___Fixture,
    querying_a_record_overlapping_just_on_the_left_border)
{
    vcf.add_record(vcf_record_5_to_10);

    std::vector<VCFRecord*> expected
        = vcf.get_all_records_overlapping_the_given_record(vcf_record_3_to_6);

    EXPECT_EQ(expected.size(), 1);
    EXPECT_EQ(*(expected[0]), vcf_record_5_to_10);
}

TEST_F(VCFTest___get_all_records_overlapping_the_given_record___Fixture,
    querying_a_record_overlapping_just_before_the_left_border)
{
    vcf.add_record(vcf_record_5_to_10);

    std::vector<VCFRecord*> expected
        = vcf.get_all_records_overlapping_the_given_record(vcf_record_3_to_5);

    EXPECT_EQ(expected.size(), 0);
}

TEST_F(VCFTest___get_all_records_overlapping_the_given_record___Fixture,
    querying_a_record_envelopping_indexed_record)
{
    vcf.add_record(vcf_record_5_to_10);

    std::vector<VCFRecord*> expected
        = vcf.get_all_records_overlapping_the_given_record(vcf_record_3_to_13);

    EXPECT_EQ(expected.size(), 1);
    EXPECT_EQ(*(expected[0]), vcf_record_5_to_10);
}

TEST_F(VCFTest___get_all_records_overlapping_the_given_record___Fixture,
    querying_a_record_in_other_chrom)
{
    vcf.add_record(vcf_record_5_to_10);

    std::vector<VCFRecord*> expected = vcf.get_all_records_overlapping_the_given_record(
        vcf_record_5_to_10_other_chrom);

    EXPECT_EQ(expected.size(), 0);
}

TEST_F(VCFTest___get_all_records_overlapping_the_given_record___Fixture,
    indexing_all_records___querying_vcf_record_5_to_10___returns_sorted_records)
{
    vcf.add_record(vcf_record_5_to_10);
    vcf.add_record(vcf_record_7_to_8);
    vcf.add_record(vcf_record_9_to_12);
    vcf.add_record(vcf_record_10_to_12);
    vcf.add_record(vcf_record_3_to_6);
    vcf.add_record(vcf_record_3_to_5);
    vcf.add_record(vcf_record_3_to_13);
    vcf.add_record(vcf_record_5_to_10_other_chrom);

    std::vector<VCFRecord*> expected
        = vcf.get_all_records_overlapping_the_given_record(vcf_record_5_to_10);

    EXPECT_EQ(expected.size(), 5);
    EXPECT_EQ(*(expected[0]), vcf_record_3_to_6);
    EXPECT_EQ(*(expected[1]), vcf_record_3_to_13);
    EXPECT_EQ(*(expected[2]), vcf_record_5_to_10);
    EXPECT_EQ(*(expected[3]), vcf_record_7_to_8);
    EXPECT_EQ(*(expected[4]), vcf_record_9_to_12);
}

TEST(VCFTest, equals)
{
    VCF vcf = create_VCF_with_default_parameters();
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord(&vcf, "chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_or_update_record_restricted_to_the_given_samples(vr, empty);
    EXPECT_EQ(vcf, vcf);
    EXPECT_FALSE(vcf != vcf);

    // different order
    VCF vcf1 = create_VCF_with_default_parameters();
    vcf1.add_record("chrom1", 5, "A", "G");
    vcf1.add_or_update_record_restricted_to_the_given_samples(vr, empty);
    vcf1.add_record("chrom1", 46, "T", "TA");
    EXPECT_EQ(vcf1, vcf1);
    EXPECT_FALSE(vcf1 != vcf1);
    EXPECT_EQ(vcf, vcf1);
    EXPECT_FALSE(vcf != vcf1);
    EXPECT_EQ(vcf1, vcf);
    EXPECT_FALSE(vcf1 != vcf);

    // same length, one different
    VCF vcf2 = create_VCF_with_default_parameters();
    vcf2.add_record("chrom1", 10, "A", "G");
    vcf2.add_or_update_record_restricted_to_the_given_samples(vr, empty);
    vcf2.add_record("chrom1", 46, "T", "TA");
    EXPECT_EQ(vcf2, vcf2);
    EXPECT_FALSE(vcf2 != vcf2);
    EXPECT_NE(vcf, vcf2);
    EXPECT_FALSE(vcf == vcf2);
    EXPECT_NE(vcf2, vcf);
    EXPECT_FALSE(vcf2 == vcf);

    // different length
    VCF vcf3 = create_VCF_with_default_parameters();
    vcf3.add_record("chrom1", 5, "A", "G");
    vcf3.add_or_update_record_restricted_to_the_given_samples(vr, empty);
    vcf3.add_record("chrom1", 46, "T", "TA");
    vcf3.add_record("chrom1", 30, "G", "CC");
    EXPECT_EQ(vcf3, vcf3);
    EXPECT_FALSE(vcf3 != vcf3);
    EXPECT_NE(vcf, vcf3);
    EXPECT_FALSE(vcf == vcf3);
    EXPECT_NE(vcf3, vcf);
    EXPECT_FALSE(vcf3 == vcf);
}

class VCFTest___header___Fixture : public ::testing::Test {
protected:
    class VCF_Mock : public VCF {
    public:
        using VCF::VCF;
        MOCK_METHOD(std::string, get_current_date, (), (const override));
    };

    VCFTest___header___Fixture()
        : vcf(&default_genotyping_options)
    {
    }

    void SetUp() override
    {
        vcf.add_samples({ "sample_1", "sample_2" });
        vcf.add_record("chrom_Z", 0, "A", "T");
        vcf.add_record("chrom_G", 0, "A", "T");
        vcf.add_record("chrom_A", 0, "A", "T");

        EXPECT_CALL(vcf, get_current_date)
            .Times(1)
            .WillOnce(Return(std::string("dummy_date")));
    }

    void TearDown() override {}

    VCF_Mock vcf;
};

TEST_F(VCFTest___header___Fixture, header)
{
    std::string expected
        = "##fileformat=VCFv4.3\n"
          "##fileDate==dummy_date\n"
          "##ALT=<ID=SNP,Description=\"SNP\">\n"
          "##ALT=<ID=PH_SNPs,Description=\"Phased SNPs\">\n"
          "##ALT=<ID=INDEL,Description=\"Insertion-deletion\">\n"
          "##ALT=<ID=COMPLEX,Description=\"Complex variant, collection of SNPs and "
          "indels\">\n"
          "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of variant\">\n"
          "##ALT=<ID=SIMPLE,Description=\"Graph bubble is simple\">\n"
          "##ALT=<ID=NESTED,Description=\"Variation site was a nested feature in the "
          "graph\">\n"
          "##ALT=<ID=TOO_MANY_ALTS,Description=\"Variation site was a multinested "
          "feature with too many alts to include all in the VCF\">\n"
          "##INFO=<ID=GRAPHTYPE,Number=1,Type=String,Description=\"Type of graph "
          "feature\">\n"
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
          "##FORMAT=<ID=MEAN_FWD_COVG,Number=A,Type=Integer,Description=\"Mean forward "
          "coverage\">\n"
          "##FORMAT=<ID=MEAN_REV_COVG,Number=A,Type=Integer,Description=\"Mean reverse "
          "coverage\">\n"
          "##FORMAT=<ID=MED_FWD_COVG,Number=A,Type=Integer,Description=\"Med forward "
          "coverage\">\n"
          "##FORMAT=<ID=MED_REV_COVG,Number=A,Type=Integer,Description=\"Med reverse "
          "coverage\">\n"
          "##FORMAT=<ID=SUM_FWD_COVG,Number=A,Type=Integer,Description=\"Sum forward "
          "coverage\">\n"
          "##FORMAT=<ID=SUM_REV_COVG,Number=A,Type=Integer,Description=\"Sum reverse "
          "coverage\">\n"
          "##FORMAT=<ID=GAPS,Number=A,Type=Float,Description=\"Number of gap bases\">\n"
          "##FORMAT=<ID=LIKELIHOOD,Number=A,Type=Float,Description=\"Likelihood\">\n"
          "##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description=\"Genotype "
          "confidence\">\n"
          "##contig=<ID=chrom_A>\n"
          "##contig=<ID=chrom_G>\n"
          "##contig=<ID=chrom_Z>\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_1\tsample_2\n";

    std::string actual = vcf.header();

    EXPECT_EQ(actual, expected);
}

class VCFTest___to_string___Fixture : public ::testing::Test {
protected:
    class VCF_DummyHeader_Mock : public VCF {
    public:
        using VCF::VCF;
        MOCK_METHOD(std::string, header, (), (const override));
    };

    std::shared_ptr<VCF_DummyHeader_Mock> vcf_with_all_records;
    VCFRecord graph_type_is_simple_sv_is_snp;
    VCFRecord graph_type_is_nested_sv_is_snp;
    VCFRecord graph_type_has_too_many_alts_sv_is_snp;
    VCFRecord graph_type_is_simple_sv_is_indel;
    VCFRecord graph_type_is_simple_sv_is_ph_snps;
    VCFRecord graph_type_is_simple_sv_is_complex;
    VCFRecord record_with_dot_allele;
    std::vector<std::string> sample_names;

    VCFTest___to_string___Fixture()
        : vcf_with_all_records(
            std::make_shared<VCF_DummyHeader_Mock>(&default_genotyping_options))
        , graph_type_is_simple_sv_is_snp(VCFRecord(vcf_with_all_records.get(), "0", 0,
              "0", "0", "SVTYPE=SNP", "GRAPHTYPE=SIMPLE"))
        , graph_type_is_nested_sv_is_snp(VCFRecord(vcf_with_all_records.get(), "0", 1,
              "0", "0", "SVTYPE=SNP", "GRAPHTYPE=NESTED"))
        , graph_type_has_too_many_alts_sv_is_snp(VCFRecord(vcf_with_all_records.get(),
              "0", 2, "0", "0", "SVTYPE=SNP", "GRAPHTYPE=TOO_MANY_ALTS"))
        , graph_type_is_simple_sv_is_indel(VCFRecord(vcf_with_all_records.get(), "0", 3,
              "0", "0", "SVTYPE=INDEL", "GRAPHTYPE=SIMPLE"))
        , graph_type_is_simple_sv_is_ph_snps(VCFRecord(vcf_with_all_records.get(), "0",
              4, "0", "0", "SVTYPE=PH_SNPs", "GRAPHTYPE=SIMPLE"))
        , graph_type_is_simple_sv_is_complex(VCFRecord(vcf_with_all_records.get(), "0",
              5, "0", "0", "SVTYPE=COMPLEX", "GRAPHTYPE=SIMPLE"))
        , record_with_dot_allele(
              VCFRecord(vcf_with_all_records.get(), "0", 6, ".", ".", ".", "."))
    {
    }

    void SetUp() override
    {
        vcf_with_all_records->add_or_update_record_restricted_to_the_given_samples(
            graph_type_is_simple_sv_is_snp, sample_names);
        vcf_with_all_records->add_or_update_record_restricted_to_the_given_samples(
            graph_type_is_nested_sv_is_snp, sample_names);
        vcf_with_all_records->add_or_update_record_restricted_to_the_given_samples(
            graph_type_has_too_many_alts_sv_is_snp, sample_names);
        vcf_with_all_records->add_or_update_record_restricted_to_the_given_samples(
            graph_type_is_simple_sv_is_indel, sample_names);
        vcf_with_all_records->add_or_update_record_restricted_to_the_given_samples(
            graph_type_is_simple_sv_is_ph_snps, sample_names);
        vcf_with_all_records->add_or_update_record_restricted_to_the_given_samples(
            graph_type_is_simple_sv_is_complex, sample_names);
        vcf_with_all_records->add_or_update_record_restricted_to_the_given_samples(
            record_with_dot_allele, sample_names);

        EXPECT_CALL(*vcf_with_all_records, header)
            .Times(1)
            .WillOnce(Return(std::string("##Dummy_header;\n")));
    }

    void TearDown() override {}
};

TEST_F(VCFTest___to_string___Fixture, graph_type_is_simple_sv_is_snp)
{
    std::string actual = vcf_with_all_records->to_string(
        true, false, false, true, false, false, true, false, false, false);

    std::string expected = "##Dummy_header;\n0\t1\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE="
                           "SIMPLE\tGT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_"
                           "REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";

    EXPECT_EQ(actual, expected);
}

TEST_F(VCFTest___to_string___Fixture, graph_type_is_nested_sv_is_snp)
{
    std::string actual = vcf_with_all_records->to_string(
        true, false, false, false, true, false, true, false, false, false);

    std::string expected = "##Dummy_header;\n0\t2\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE="
                           "NESTED\tGT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_"
                           "REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";

    EXPECT_EQ(actual, expected);
}

TEST_F(VCFTest___to_string___Fixture, graph_type_has_too_many_alts_sv_is_snp)
{
    std::string actual = vcf_with_all_records->to_string(
        true, false, false, false, false, true, true, false, false, false);

    std::string expected = "##Dummy_header;\n0\t3\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE="
                           "TOO_MANY_ALTS\tGT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:"
                           "MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";

    EXPECT_EQ(actual, expected);
}

TEST_F(VCFTest___to_string___Fixture, graph_type_is_simple_sv_is_indel)
{
    std::string actual = vcf_with_all_records->to_string(
        true, false, false, true, false, false, false, true, false, false);

    std::string expected = "##Dummy_header;\n0\t4\t.\t0\t0\t.\t.\tSVTYPE=INDEL;"
                           "GRAPHTYPE=SIMPLE\tGT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_"
                           "COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";

    EXPECT_EQ(actual, expected);
}

TEST_F(VCFTest___to_string___Fixture, graph_type_is_simple_sv_is_ph_snps)
{
    std::string actual = vcf_with_all_records->to_string(
        true, false, false, true, false, false, false, false, true, false);

    std::string expected = "##Dummy_header;\n0\t5\t.\t0\t0\t.\t.\tSVTYPE=PH_SNPs;"
                           "GRAPHTYPE=SIMPLE\tGT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_"
                           "COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";

    EXPECT_EQ(actual, expected);
}

TEST_F(VCFTest___to_string___Fixture, graph_type_is_simple_sv_is_complex)
{
    std::string actual = vcf_with_all_records->to_string(
        true, false, false, true, false, false, false, false, false, true);

    std::string expected = "##Dummy_header;\n0\t6\t.\t0\t0\t.\t.\tSVTYPE=COMPLEX;"
                           "GRAPHTYPE=SIMPLE\tGT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_"
                           "COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";

    EXPECT_EQ(actual, expected);
}

TEST_F(VCFTest___to_string___Fixture, record_with_dot_allele)
{
    std::string actual = vcf_with_all_records->to_string(
        true, false, true, false, false, false, false, false, false, false);

    std::string expected
        = "##Dummy_header;\n0\t7\t.\t.\t.\t.\t.\t.;.\tGT:MEAN_FWD_COVG:MEAN_REV_COVG:"
          "MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";

    EXPECT_EQ(actual, expected);
}

TEST_F(VCFTest___to_string___Fixture, all_records_filtered_out)
{
    std::string actual = vcf_with_all_records->to_string(
        true, false, false, false, false, false, false, false, false, false);

    std::string expected = "##Dummy_header;\n";

    EXPECT_EQ(actual, expected);
}

TEST_F(VCFTest___to_string___Fixture, no_records_filtered_out)
{
    std::string actual = vcf_with_all_records->to_string(
        true, false, true, true, true, true, true, true, true, true);

    std::string expected = "##Dummy_header;\n";
    expected
        += "0\t1\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE=SIMPLE\tGT:MEAN_FWD_COVG:MEAN_"
           "REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";
    expected
        += "0\t2\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE=NESTED\tGT:MEAN_FWD_COVG:MEAN_"
           "REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";
    expected
        += "0\t3\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE=TOO_MANY_ALTS\tGT:MEAN_FWD_COVG:"
           "MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";
    expected
        += "0\t4\t.\t0\t0\t.\t.\tSVTYPE=INDEL;GRAPHTYPE=SIMPLE\tGT:MEAN_FWD_COVG:MEAN_"
           "REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";
    expected
        += "0\t5\t.\t0\t0\t.\t.\tSVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE\tGT:MEAN_FWD_COVG:"
           "MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";
    expected
        += "0\t6\t.\t0\t0\t.\t.\tSVTYPE=COMPLEX;GRAPHTYPE=SIMPLE\tGT:MEAN_FWD_COVG:"
           "MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";
    expected += "0\t7\t.\t.\t.\t.\t.\t.;.\tGT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:"
                "MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS\t\n";

    EXPECT_EQ(actual, expected);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OLD serialization TESTS
// WILL NOT BE READDED AS WE DO NOT NEED FULL SERIALIZATION OF VCFs (ONLY SAVE, WHICH
// JUST WRAPS VCF::to_string()), WHICH IS TESTED WITH A NEW TEST
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class VCFTest___serialization___Fixture : public ::testing::Test {
// protected:
//    VCF vcf_with_zero_records;
//    VCF vcf_with_one_record;
//    VCF vcf_with_three_records;
//    void SetUp() override {
//        {
//            vcf_with_one_record.add_record("chrom1", 5, "A", "G",
//            "GRAPHTYPE=SIMPLE;SVTYPE=SNP");
//        }
//
//        {
//            vcf_with_three_records.add_record("chrom1", 5, "A", "G",
//            "GRAPHTYPE=SIMPLE;SVTYPE=SNP");
//            vcf_with_three_records.add_record("chrom1", 46, "T", "TA",
//            "GRAPHTYPE=SIMPLE;SVTYPE=SNP"); VCFRecord vcf_record = VCFRecord("chrom1",
//            79, "C", "G", "GRAPHTYPE=SIMPLE;SVTYPE=SNP"); std::vector<std::string>
//            empty_sample_names = {}; vcf_with_three_records.add_record(vcf_record,
//            empty_sample_names);
//        }
//    }
//
//    void TearDown() override {
//    }
//};
//
// TEST_F(VCFTest___serialization___Fixture,
// save_vcf_with_zero_records___load_vcf___expect_equal_vcf) {
//    vcf_with_zero_records.save("vcf_serialization_test_zero.vcf");
//
//    VCF actual;
//    actual.load("vcf_serialization_test_zero.vcf");
//
//    VCF& expected = vcf_with_zero_records;
//    EXPECT_EQ(actual, expected);
//}
//
// TEST_F(VCFTest___serialization___Fixture,
// save_vcf_with_one_record___load_vcf___expect_equal_vcf) {
//    vcf_with_one_record.save("vcf_serialization_test_one.vcf");
//
//    VCF actual;
//    actual.load("vcf_serialization_test_one.vcf");
//
//    VCF& expected = vcf_with_one_record;
//    EXPECT_EQ(actual, expected);
//}
//
//
// TEST_F(VCFTest___serialization___Fixture,
// save_vcf_with_three_records___load_vcf___expect_equal_vcf) {
//    vcf_with_three_records.save("vcf_serialization_test_three.vcf");
//
//    VCF actual;
//    actual.load("vcf_serialization_test_three.vcf");
//
//    VCF& expected = vcf_with_three_records;
//    EXPECT_EQ(actual, expected);
//}
//
//
// TEST(VCFTest, filter) {
//    VCF vcf, vcf1, vcf2, vcf3, vcf4;
//    vcf.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
//    vcf.add_record("chrom1", 46, "T", "TA", "SVTYPE=INDEL;GRAPHTYPE=NESTED");
//    vcf.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE");
//    vcf.add_record("chrom1", 79, "CTT", "ATA", "SVTYPE=PH_SNPs;GRAPHTYPE=NESTED");
//    vcf.save("vcf_filter_test.vcf", false, true, false, false, true, false, true,
//    false);
//
//    vcf1.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
//    vcf1.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE");
//    vcf2.load("vcf_filter_test.vcf");
//    EXPECT_EQ(vcf2, vcf1);
//
//    vcf.save("vcf_filter_test.vcf", false, true, true, false, false, false, true,
//    false); vcf3.add_record("chrom1", 79, "CTT", "GTA",
//    "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE"); vcf3.add_record("chrom1", 79, "CTT", "ATA",
//    "SVTYPE=PH_SNPs;GRAPHTYPE=NESTED"); vcf4.load("vcf_filter_test.vcf");
//    EXPECT_EQ(vcf3, vcf4);
//}

class VCFTest___save___Fixture : public ::testing::Test {
protected:
    class VCF_Mock : public VCF {
    public:
        using VCF::VCF;
        MOCK_METHOD(std::string, to_string,
            (bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage,
                bool output_dot_allele, bool graph_is_simple, bool graph_is_nested,
                bool graph_has_too_many_alts, bool sv_type_is_snp,
                bool sv_type_is_indel, bool sv_type_is_ph_snps,
                bool sv_type_is_complex),
            (override));
    };

    VCFTest___save___Fixture()
        : vcf(&default_genotyping_options)
    {
    }

    VCF_Mock vcf;
};

// TODO : improve this test with dependency injection in the file handler
TEST_F(VCFTest___save___Fixture, save_true_then_false_flags)
{
    EXPECT_CALL(
        vcf, to_string(true, false, true, false, true, false, true, false, true, false))
        .Times(1);

    vcf.save(
        "/dev/null", true, false, true, false, true, false, true, false, true, false);
}

TEST_F(VCFTest___save___Fixture, save_false_then_true_flags)
{
    EXPECT_CALL(
        vcf, to_string(false, true, false, true, false, true, false, true, false, true))
        .Times(1);

    vcf.save(
        "/dev/null", false, true, false, true, false, true, false, true, false, true);
}

TEST(VCFTest, concatenate_VCFs)
{
    VCF::concatenate_VCFs({ "../../test/test_cases/concatenate_VCFs/fake_vcf1.vcf",
                              "../../test/test_cases/concatenate_VCFs/fake_vcf2.vcf",
                              "../../test/test_cases/concatenate_VCFs/fake_vcf3.vcf" },
        "concatenated_vcf.vcf");
    std::vector<std::string> actual
        = get_vector_of_strings_from_file("concatenated_vcf.vcf");

    std::vector<std::string> expected
        = { "#dummy_header", "line_1", "line_2", "line_3" };

    EXPECT_EQ(actual, expected);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PREVIOUS TESTS FROM VCF_RECORD FOLLOW
// COMMENTED OUT == NO NEED ANYMORE AND I PUT THE REASON WHY IT IS NOT NEEDED
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// REASON WHY THESE ARE COMMENTED OUT: NO NEED ANYMORE, FORMATS ARE NOT VARIABLE
// ANYMORE, WHICH RENDERS THE CODE SIMPLER
// TEST(VCFTest, add_formats) {
//    VCF vcf = create_VCF_with_default_parameters();
//    std::vector<std::string> formats = {"GT", "LIKELIHOOD", "GT_CONF",
//    "MEAN_FWD_COVG", "MEAN_REV_COVG", "GAPS"};
//
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1",
//    46, "CTT", "TA");
//
//    vcf.add_formats(formats);
//
//    for (const auto& record: vcf.records){
//        for (const auto &f: formats){
//            EXPECT_TRUE(std::find(record->format.begin(), record->format.end(),
//            f)!=record->format.end());
//        }
//    }
//}