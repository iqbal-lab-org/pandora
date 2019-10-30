//#include "gtest/gtest.h"
//#include "test_macro.cpp"
//#include "vcf.h"
//#include "vcfrecord.h"
//#include "interval.h"
//#include "localnode.h"
//#include <stdint.h>
//#include <iostream>
//
//
//using namespace std;
//
//TEST(VCFTest, add_record_with_values) {
//
//    VCF vcf;
//    EXPECT_EQ((uint) 0, vcf.records.size());
//    vcf.add_record("chrom1", 5, "A", "G");
//    EXPECT_EQ((uint) 1, vcf.records.size());
//}
//
//TEST(VCFTest, add_record_twice_with_values) {
//
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 5, "A", "G");
//    EXPECT_EQ((uint) 1, vcf.records.size());
//}
//
//TEST(VCFTest, add_two_records_with_values) {
//
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    EXPECT_EQ((uint) 2, vcf.records.size());
//}
//
//TEST(VCFTest, add_two_records_and_a_repeat_with_values) {
//
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 5, "A", "G");
//    EXPECT_EQ((uint) 2, vcf.records.size());
//}
//
//TEST(VCFTest, add_record_by_record) {
//    VCF vcf;
//    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
//    std::vector<std::string> empty = {};
//    vcf.add_record(vr, empty);
//    EXPECT_EQ((uint) 1, vcf.records.size());
//}
//
//TEST(VCFTest, add_record_by_record_and_values) {
//    VCF vcf;
//    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
//    std::vector<std::string> empty = {};
//    vcf.add_record(vr, empty);
//    vcf.add_record("chrom1", 79, "C", "G");
//    EXPECT_EQ((uint) 1, vcf.records.size());
//}
//
//TEST(VCFTest, add_record_by_values_and_record) {
//    VCF vcf;
//    vcf.add_record("chrom1", 79, "C", "G");
//    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
//    std::vector<std::string> empty = {};
//    vcf.add_record(vr, empty);
//    EXPECT_EQ((uint) 1, vcf.records.size());
//}
//
//TEST(VCFTest, add_record_by_record_returned_by_reference) {
//    VCF vcf;
//    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
//    std::vector<std::string> empty = {};
//    VCFRecord &ref_vr = vcf.add_record(vr, empty);
//    EXPECT_EQ(ref_vr.chrom, "chrom1");
//    EXPECT_EQ(ref_vr.pos, (uint) 79);
//}
//
//TEST(VCFTest, add_samples_empty) {
//    VCF vcf;
//    std::vector<std::string> samples;
//    vcf.add_samples(samples);
//    EXPECT_EQ(vcf.samples.size(), (uint)0);
//    EXPECT_EQ(vcf.records.size(), (uint)0);
//}
//
//TEST(VCFTest, add_samples_simple) {
//    VCF vcf;
//    std::vector<std::string> samples = {"hello", "there", "people"};
//    vcf.add_samples(samples);
//    EXPECT_ITERABLE_EQ(std::vector<std::string>, samples, vcf.samples);
//    EXPECT_EQ(vcf.records.size(), (uint)0);
//}
//
//TEST(VCFTest, add_samples_with_record) {
//    VCF vcf;
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 5, "A", "G");
//
//    std::vector<std::string> samples = {"hello", "there", "people"};
//    vcf.add_samples(samples);
//
//    std::vector<std::string> exp_samples = {"sample", "hello", "there", "people"};
//
//    EXPECT_ITERABLE_EQ(std::vector<std::string>, exp_samples, vcf.samples);
//    EXPECT_EQ(vcf.records.size(), (uint)1);
//    EXPECT_EQ(vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size(), exp_samples.size());
//}
//
//TEST(VCFTest, add_sample_gt) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_record("chrom1", 79, "C", "A");
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "T", "TA");
//    uint j = 1;
//    EXPECT_EQ(j, vcf.samples.size());
//    EXPECT_EQ(j, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ(j, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_EQ(j, vcf.records[2]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_EQ(j, vcf.records[3]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].end());
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 79, "C", "C");
//    EXPECT_EQ(j, vcf.samples.size());
//    EXPECT_EQ(j, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ(j, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_EQ(j, vcf.records[2]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint16_t) 0, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ(j, vcf.records[3]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint16_t) 0, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//}
//
//TEST(VCFTest, add_record_by_record_with_existing_sample) {
//    VCF vcf;
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "T", "TA");
//    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
//    std::vector<std::string> empty = {};
//    VCFRecord &ref_vr = vcf.add_record(vr, empty);
//    EXPECT_EQ(ref_vr.chrom, "chrom1");
//    EXPECT_EQ(ref_vr.pos, (uint) 79);
//    EXPECT_EQ(ref_vr.sampleIndex_to_format_to_sampleInfo.size(), (uint) 1);
//}
//
//TEST(VCFTest, add_record_by_record_with_same_existing_sample) {
//    VCF vcf;
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "T", "TA");
//    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
//    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
//    vr.sampleIndex_to_format_to_sampleInfo[0]["GT"] = {1};
//    std::vector<std::string> samples = {"sample"};
//    VCFRecord &ref_vr = vcf.add_record(vr, samples);
//    EXPECT_EQ(ref_vr.chrom, "chrom1");
//    EXPECT_EQ(ref_vr.pos, (uint) 79);
//    EXPECT_EQ(ref_vr.sampleIndex_to_format_to_sampleInfo.size(), (uint) 1);
//    EXPECT_ITERABLE_EQ(vector<std::string>, samples, vcf.samples);
//    EXPECT_FALSE(ref_vr.sampleIndex_to_format_to_sampleInfo[0].find("GT") == ref_vr.sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_EQ(ref_vr.sampleIndex_to_format_to_sampleInfo[0]["GT"][0], (uint16_t)1);
//}
//
//TEST(VCFTest, add_record_by_record_with_different_existing_sample) {
//    VCF vcf;
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "T", "TA");
//    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
//    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
//    vr.sampleIndex_to_format_to_sampleInfo[0]["GT"] = {1};
//    std::vector<std::string> samples = {"sample1"};
//    VCFRecord &ref_vr = vcf.add_record(vr, samples);
//    EXPECT_EQ(ref_vr.chrom, "chrom1");
//    EXPECT_EQ(ref_vr.pos, (uint) 79);
//    EXPECT_EQ(ref_vr.sampleIndex_to_format_to_sampleInfo.size(), (uint) 2);
//    std::vector<std::string> exp_samples = {"sample", "sample1"};
//    EXPECT_ITERABLE_EQ(vector<std::string>, exp_samples, vcf.samples);
//    EXPECT_TRUE(ref_vr.sampleIndex_to_format_to_sampleInfo[0].find("GT") == ref_vr.sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_FALSE(ref_vr.sampleIndex_to_format_to_sampleInfo[1].find("GT") == ref_vr.sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_EQ(ref_vr.sampleIndex_to_format_to_sampleInfo[1]["GT"][0], (uint16_t)1);
//}
//
//TEST(VCFTest, add_sample_ref_alleles) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_record("chrom1", 79, "C", "A");
//    vcf.add_record("chrom2", 30, "C", "A");
//
//    vcf.set_sample_gt_to_ref_allele_for_records_in_the_interval("sample", "chrom1", 15, 78);
//    EXPECT_EQ((uint) 1, vcf.samples.size());
//    EXPECT_EQ((uint) 5, vcf.records.size());
//    EXPECT_EQ((uint) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_EQ((uint) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint) 1, vcf.records[2]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_EQ((uint) 1, vcf.records[3]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_EQ((uint) 1, vcf.records[4]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].end());
//
//    vcf.set_sample_gt_to_ref_allele_for_records_in_the_interval("sample2", "chrom1", 5, 46);
//    EXPECT_EQ((uint) 2, vcf.samples.size());
//    EXPECT_EQ((uint) 5, vcf.records.size());
//    EXPECT_EQ((uint) 2, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint16_t) 0, vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    EXPECT_EQ((uint) 2, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1].find("GT") == vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1].end());
//    EXPECT_EQ((uint) 2, vcf.records[2]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1].find("GT") == vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1].end());
//    EXPECT_EQ((uint) 2, vcf.records[3]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1].find("GT") == vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1].end());
//    EXPECT_EQ((uint) 2, vcf.records[4]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1].find("GT") == vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1].end());
//}
//
//TEST(VCFTest, reorder_add_record_and_sample) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample1", "chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample2", "chrom1", 79, "C", "C");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample1", "chrom1", 79, "C", "A");
//
//    vcf.sort_records();
//
//    EXPECT_EQ((uint) 2, vcf.samples.size());
//    EXPECT_EQ((uint) 4, vcf.records.size());
//    EXPECT_EQ((uint) 2, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 2, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 2, vcf.records[2]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 2, vcf.records[3]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_TRUE(vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 1, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_TRUE(vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].find("GT") == vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_TRUE(vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1].find("GT") == vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1].end());
//    EXPECT_TRUE(vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1].find("GT") == vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1].end());
//    EXPECT_EQ((uint16_t) 0, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    EXPECT_EQ((uint16_t) 0, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//
//}
//
//
//TEST(VCFTest, append_vcf_simple_case) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_record("chrom1", 79, "C", "A");
//
//    VCF new_vcf;
//    new_vcf.add_record("chrom2", 5, "A", "G");
//    new_vcf.add_record("chrom2", 46, "T", "TA");
//    new_vcf.add_record("chrom2", 79, "C", "G");
//    new_vcf.add_record("chrom2", 79, "C", "A");
//
//    vcf.append_vcf(new_vcf);
//    EXPECT_EQ((uint) 8, vcf.records.size());
//    for (uint i = 0; i < 4; ++i) {
//        EXPECT_EQ(vcf.records[i]->chrom, "chrom1");
//    }
//    for (uint i = 4; i < 8; ++i) {
//        EXPECT_EQ(vcf.records[i]->chrom, "chrom2");
//    }
//    EXPECT_EQ((uint) 5, vcf.records[4]->pos);
//    EXPECT_EQ("TA", vcf.records[5]->alts[0]);
//    EXPECT_EQ((uint) 79, vcf.records[6]->pos);
//    EXPECT_EQ("A", vcf.records[7]->alts[0]);
//}
//
//TEST(VCFTest, append_vcf_some_duplicate_records) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_record("chrom1", 79, "C", "A");
//
//    VCF new_vcf;
//    new_vcf.add_record("chrom2", 5, "A", "G");
//    new_vcf.add_record("chrom1", 46, "T", "TA");
//    new_vcf.add_record("chrom2", 79, "C", "G");
//    new_vcf.add_record("chrom1", 79, "C", "A");
//
//    vcf.append_vcf(new_vcf);
//    EXPECT_EQ((uint) 6, vcf.records.size());
//    for (uint i = 0; i < 4; ++i) {
//        EXPECT_EQ(vcf.records[i]->chrom, "chrom1");
//    }
//    for (uint i = 4; i < 6; ++i) {
//        EXPECT_EQ(vcf.records[i]->chrom, "chrom2");
//    }
//    EXPECT_EQ((uint) 5, vcf.records[4]->pos);
//    EXPECT_EQ((uint) 79, vcf.records[5]->pos);
//}
//
//TEST(VCFTest, append_vcf_one_sample) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_record("chrom1", 79, "C", "A");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 79, "C", "G");
//
//    VCF new_vcf;
//    new_vcf.add_record("chrom2", 5, "A", "G");
//    new_vcf.add_record("chrom1", 46, "T", "TA");
//    new_vcf.add_record("chrom2", 79, "C", "G");
//    new_vcf.add_record("chrom1", 79, "C", "A");
//
//    vcf.append_vcf(new_vcf);
//    EXPECT_EQ((uint) 1, vcf.samples.size());
//    EXPECT_EQ("sample", vcf.samples[0]);
//    EXPECT_EQ((uint) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 1, vcf.records[5]->sampleIndex_to_format_to_sampleInfo.size());
//    bool found_gt = vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_TRUE(found_gt);
//    EXPECT_EQ((uint) 1, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    found_gt = vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//}
//
//TEST(VCFTest, append_vcf_one_sample_in_new_vcf) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_record("chrom1", 79, "C", "A");
//
//    VCF new_vcf;
//    new_vcf.add_record("chrom2", 5, "A", "G");
//    new_vcf.add_record("chrom1", 46, "T", "TA");
//    new_vcf.add_record("chrom2", 79, "C", "G");
//    new_vcf.add_record("chrom1", 79, "C", "A");
//    new_vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 5, "A", "G");
//
//    vcf.append_vcf(new_vcf);
//    EXPECT_EQ((uint) 1, vcf.samples.size());
//    EXPECT_EQ("sample", vcf.samples[0]);
//    EXPECT_EQ((uint) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 1, vcf.records[5]->sampleIndex_to_format_to_sampleInfo.size());
//    bool found_gt = vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_TRUE(found_gt);
//    EXPECT_EQ((uint) 1, vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    found_gt = vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//}
//
//TEST(VCFTest, append_vcf_shared_sample) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_record("chrom1", 79, "C", "A");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "T", "TA");
//
//    VCF new_vcf;
//    new_vcf.add_record("chrom2", 5, "A", "G");
//    new_vcf.add_record("chrom1", 46, "T", "TA");
//    new_vcf.add_record("chrom2", 79, "C", "G");
//    new_vcf.add_record("chrom1", 79, "C", "A");
//    new_vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "T", "TA");
//
//    vcf.append_vcf(new_vcf);
//    EXPECT_EQ((uint) 1, vcf.samples.size());
//    EXPECT_EQ("sample", vcf.samples[0]);
//    EXPECT_EQ((uint) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 1, vcf.records[5]->sampleIndex_to_format_to_sampleInfo.size());
//    bool found_gt = vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_TRUE(found_gt);
//    EXPECT_EQ((uint) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    found_gt = vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//}
//
//TEST(VCFTest, append_vcf_shared_samples_different_order) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_record("chrom1", 79, "C", "A");
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "T", "TA");
//
//    VCF new_vcf;
//    new_vcf.add_record("chrom1", 79, "C", "A");
//    new_vcf.add_record("chrom2", 5, "A", "G");
//    new_vcf.add_record("chrom1", 46, "T", "TA");
//    new_vcf.add_record("chrom2", 79, "C", "G");
//    new_vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample1", "chrom1", 46, "T", "T");
//    new_vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample1", "chrom1", 79, "C", "A");
//
//    vcf.append_vcf(new_vcf);
//
//    EXPECT_EQ((uint) 2, vcf.samples.size());
//    vector<string> v = {"sample", "sample1"};
//    EXPECT_ITERABLE_EQ(vector<string>, v, vcf.samples);
//    EXPECT_EQ((uint) 2, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 2, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 2, vcf.records[2]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 2, vcf.records[3]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 2, vcf.records[4]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 2, vcf.records[5]->sampleIndex_to_format_to_sampleInfo.size());
//
//    vector<uint16_t> alt_gt = {1};
//    vector<uint16_t> ref_gt = {0};
//
//    EXPECT_FALSE(vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_FALSE(vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1].find("GT") != vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1].end());
//
//    EXPECT_TRUE(vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_TRUE(vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1].end());
//    EXPECT_ITERABLE_EQ(vector<uint16_t>, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"], alt_gt);
//    EXPECT_ITERABLE_EQ(vector<uint16_t>, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["GT"], ref_gt);
//
//    EXPECT_FALSE(vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_FALSE(vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1].find("GT") != vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1].end());
//
//    EXPECT_FALSE(vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_TRUE(vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1].find("GT") != vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1].end());
//    EXPECT_ITERABLE_EQ(vector<uint16_t>, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["GT"], alt_gt);
//
//    EXPECT_FALSE(vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_FALSE(vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1].find("GT") != vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1].end());
//
//    EXPECT_FALSE(vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0].end());
//    EXPECT_FALSE(vcf.records[5]->sampleIndex_to_format_to_sampleInfo[1].find("GT") != vcf.records[5]->sampleIndex_to_format_to_sampleInfo[1].end());
//}
//
//TEST(VCFTest, sort_records) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_record("chrom1", 79, "C", "A");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "T", "TA");
//    vcf.add_record("chrom1", 79, "C", "A");
//    vcf.add_record("chrom2", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom2", 79, "C", "G");
//    vcf.sort_records();
//
//    EXPECT_EQ((uint) 6, vcf.records.size());
//    for (uint i = 0; i < 4; ++i) {
//        EXPECT_EQ("chrom1", vcf.records[i]->chrom);
//    }
//    for (uint i = 4; i < 6; ++i) {
//        EXPECT_EQ("chrom2", vcf.records[i]->chrom);
//    }
//    EXPECT_EQ((uint) 5, vcf.records[0]->pos);
//    EXPECT_EQ((uint) 5, vcf.records[4]->pos);
//    EXPECT_EQ((uint) 46, vcf.records[1]->pos);
//    EXPECT_EQ((uint) 79, vcf.records[2]->pos);
//    EXPECT_EQ((uint) 79, vcf.records[3]->pos);
//    EXPECT_EQ((uint) 79, vcf.records[5]->pos);
//    EXPECT_EQ("G", vcf.records[3]->alts[0]);
//    EXPECT_EQ("G", vcf.records[5]->alts[0]);
//}
//
//TEST(VCFTest, pos_in_range) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_record("chrom1", 79, "C", "A");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    vcf.add_record("chrom2", 20, "A", "G");
//    vcf.add_record("chrom2", 79, "C", "G");
//    //vcf.sort();
//
//    EXPECT_TRUE(vcf.pos_in_range(4, 6, "chrom1"));
//    EXPECT_FALSE(vcf.pos_in_range(5, 6, "chrom1"));
//    EXPECT_FALSE(vcf.pos_in_range(4, 5, "chrom1"));
//    EXPECT_FALSE(vcf.pos_in_range(4, 6, "chrom2"));
//
//    EXPECT_TRUE(vcf.pos_in_range(45, 47, "chrom1"));
//    EXPECT_FALSE(vcf.pos_in_range(46, 47, "chrom1"));
//    EXPECT_FALSE(vcf.pos_in_range(45, 46, "chrom1"));
//    EXPECT_FALSE(vcf.pos_in_range(45, 47, "chrom2"));
//
//    EXPECT_TRUE(vcf.pos_in_range(78, 80, "chrom1"));
//    EXPECT_FALSE(vcf.pos_in_range(79, 80, "chrom1"));
//    EXPECT_FALSE(vcf.pos_in_range(78, 79, "chrom1"));
//    EXPECT_TRUE(vcf.pos_in_range(78, 80, "chrom2"));
//
//}
//
//TEST(VCFTest, genotype) {
//    // not a snp site
//    // missing count data
//    // not confident
//    // confident and have right gt
//    // confident and have wrong gt
//    // one sample needs regenotyping and the other doesn't
//    VCF vcf;
//
//    vcf.add_record("chrom2", 79, "C", "G");
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 2, "T", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 5, "A", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 79, "C", "A");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 20, "A", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 79, "C", "C");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 80, "A", "C");
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1", 2, "T", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1", 5, "A", "A");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1", 79, "C", "A");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2", 20, "A", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2", 79, "C", "C");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2", 80, "A", "A");
//
//    vcf.sort_records();
//    std::vector<float> f = {0.0, 0.0};
//
//    // record 0, not a snp site
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 10};
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 20};
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"] = {1, 15};
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"] = {2, 24};
//    vcf.records[0]->set_format(0,"GAPS", f);
//    vcf.records[0]->set_format(1,"GAPS", f);
//
//
//    // record 1, different genotypes but both correct
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 10};
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 20};
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"] = {10, 1};
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"] = {21, 2};
//    vcf.records[1]->set_format(0,"GAPS", f);
//    vcf.records[1]->set_format(1,"GAPS", f);
//
//    // record 2, same genotypes first correct
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 10};
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 20};
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"] = {10, 1};
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"] = {21, 2};
//    vcf.records[2]->set_format(0,"GAPS", f);
//    vcf.records[2]->set_format(1,"GAPS", f);
//
//    // record 3, same genotypes both wrong
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {20, 1};
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {21, 2};
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"] = {10, 1};
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"] = {21, 2};
//    vcf.records[3]->set_format(0,"GAPS", f);
//    vcf.records[3]->set_format(1,"GAPS", f);
//
//    // record 4, missing count data for first sample
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 10};
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {20};
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"] = {10, 1};
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"] = {21, 2};
//    vcf.records[4]->set_format(0,"GAPS", f);
//    vcf.records[4]->set_format(1,"GAPS", f);
//
//    // record 5, not confident for second sample
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 10};
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 20};
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"] = {2, 1};
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"] = {4, 2};
//    vcf.records[5]->set_format(0,"GAPS", f);
//    vcf.records[5]->set_format(1,"GAPS", f);
//
//    vcf.genotype({30, 30}, 0.01, 30, 0, 1, 0, 0, true);
//
//    // not genotyped first record
//    EXPECT_EQ((uint16_t) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    bool found_confidence = vcf.records[0]->sampleIndex_to_format_to_sampleGenotypedInfo[0].find("GT_CONF") != vcf.records[0]->sampleIndex_to_format_to_sampleGenotypedInfo[0].end();
//    EXPECT_FALSE(found_confidence);
//    found_confidence = vcf.records[0]->sampleIndex_to_format_to_sampleGenotypedInfo[1].find("GT_CONF") != vcf.records[0]->sampleIndex_to_format_to_sampleGenotypedInfo[1].end();
//    EXPECT_FALSE(found_confidence);
//
//    // both correct
//    EXPECT_EQ((uint) 2, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
//    bool found_gt = vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_TRUE(found_gt);
//    found_gt = vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1].end();
//    EXPECT_TRUE(found_gt);
//    EXPECT_EQ((uint) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["GT"].size());
//    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    // first correct
//    EXPECT_EQ((uint16_t) 1, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 0, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    // both wrong
//    EXPECT_EQ((uint16_t) 0, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 0, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    // first missing data
//    EXPECT_EQ((uint) 0, vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint16_t) 0, vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    // second not confident
//    EXPECT_EQ((uint16_t) 1, vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint) 0, vcf.records[5]->sampleIndex_to_format_to_sampleInfo[1]["GT"].size());
//}
//
//TEST(VCFTest, genotype_with_all_sites) {
//    // not a snp site
//    // missing count data
//    // not confident
//    // confident and have right gt
//    // confident and have wrong gt
//    // one sample needs regenotyping and the other doesn't
//    VCF vcf;
//
//    vcf.add_record("chrom2", 79, "CC", "GC");
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 2, "T", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 5, "AC", "GC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 79, "CC", "AC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 20, "AC", "GC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 79, "CC", "CC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 80, "AC", "CC");
//
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1", 2, "T", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1", 5, "AC", "AC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom1", 79, "CC", "AC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2", 20, "AC", "GC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2", 79, "CC", "CC");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("asample", "chrom2", 80, "AC", "AC");
//
//    vcf.sort_records();
//    std::vector<float> f = {0.0, 0.0};
//
//    // record 0, not a snp site
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(0);
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(1);
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(10);
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(20);
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(15);
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(24);
//    vcf.records[0]->set_format(0,"GAPS", f);
//    vcf.records[0]->set_format(1,"GAPS", f);
//
//    // record 1, different genotypes but both correct
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(0);
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(1);
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(10);
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(20);
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(10);
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(21);
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.records[1]->set_format(0,"GAPS", f);
//    vcf.records[1]->set_format(1,"GAPS", f);
//
//    // record 2, same genotypes first correct
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(0);
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(1);
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(10);
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(20);
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(10);
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(21);
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.records[2]->set_format(0,"GAPS", f);
//    vcf.records[2]->set_format(1,"GAPS", f);
//
//    // record 3, same genotypes both wrong
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(20);
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(21);
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(1);
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(2);
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(10);
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(21);
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.records[3]->set_format(0,"GAPS", f);
//    vcf.records[3]->set_format(1,"GAPS", f);
//
//    // record 4, missing count data for first sample
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(0);
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(10);
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(20);
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(10);
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(21);
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.records[4]->set_format(0,"GAPS", f);
//    vcf.records[4]->set_format(1,"GAPS", f);
//
//    // record 5, not confident for second sample
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(0);
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(1);
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"].push_back(10);
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"].push_back(20);
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(2);
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(4);
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"].push_back(1);
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"].push_back(2);
//    vcf.records[5]->set_format(0,"GAPS", f);
//    vcf.records[5]->set_format(1,"GAPS", f);
//
//    bool snps_only = false;
//    vcf.genotype({30, 30}, 0.01, 30, 0, 1, 0, 0, snps_only);
//
//    // first record now genotyped
//    EXPECT_EQ((uint16_t) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    bool found_confidence = vcf.records[0]->sampleIndex_to_format_to_sampleGenotypedInfo[0].find("GT_CONF") != vcf.records[0]->sampleIndex_to_format_to_sampleGenotypedInfo[0].end();
//    EXPECT_TRUE(found_confidence);
//    found_confidence = vcf.records[0]->sampleIndex_to_format_to_sampleGenotypedInfo[1].find("GT_CONF") != vcf.records[0]->sampleIndex_to_format_to_sampleGenotypedInfo[1].end();
//    EXPECT_TRUE(found_confidence);
//    // both correct
//    EXPECT_EQ((uint) 2, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
//    bool found_gt = vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_TRUE(found_gt);
//    found_gt = vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1].end();
//    EXPECT_TRUE(found_gt);
//    EXPECT_EQ((uint) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["GT"].size());
//    EXPECT_EQ((uint16_t) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    // first correct
//    EXPECT_EQ((uint16_t) 1, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 0, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    // both wrong
//    EXPECT_EQ((uint16_t) 0, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 0, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    // first missing data
//    EXPECT_EQ((uint) 0, vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint16_t) 0, vcf.records[4]->sampleIndex_to_format_to_sampleInfo[1]["GT"][0]);
//    // second not confident
//    EXPECT_EQ((uint16_t) 1, vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint) 0, vcf.records[5]->sampleIndex_to_format_to_sampleInfo[1]["GT"].size());
//
//}
//
//TEST(VCFTest, clean) {
//    VCF vcf;
//
//    VCFRecord dummy;
//    std::vector<std::string> empty = {};
//    vcf.add_record(dummy, empty);
//    vcf.add_record("chrom1", 79, "C", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 2, "T", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 5, "A", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 79, "C", "A");
//    vcf.records[2]->clear();
//    EXPECT_EQ((uint) 5, vcf.records.size());
//
//    vcf.clean();
//    EXPECT_EQ((uint) 3, vcf.records.size());
//    EXPECT_EQ((uint) 79, vcf.records[0]->pos);
//    EXPECT_EQ((uint) 1, vcf.records[0]->alts.size());
//    EXPECT_EQ("G", vcf.records[0]->alts[0]);
//    EXPECT_EQ((uint) 5, vcf.records[1]->pos);
//    EXPECT_EQ((uint) 79, vcf.records[2]->pos);
//    EXPECT_EQ((uint) 1, vcf.records[2]->alts.size());
//    EXPECT_EQ("A", vcf.records[2]->alts[0]);
//}
//
//TEST(VCFTest, add_formats) {
//    VCF vcf;
//    std::vector<std::string> formats = {"GT", "LIKELIHOOD", "GT_CONF", "MEAN_FWD_COVG", "MEAN_REV_COVG", "GAPS"};
//
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "CTT", "TA");
//
//    vcf.add_formats(formats);
//
//    for (const auto& record: vcf.records){
//        for (const auto &f: formats){
//            EXPECT_TRUE(std::find(record->format.begin(), record->format.end(), f)!=record->format.end());
//        }
//    }
//}
//
//TEST(VCFTest, merge_multi_allelic) {
//    VCF vcf;
//    // no gt
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 5, "A", "C");
//    // gt
//    vcf.add_record("chrom1", 46, "CTT", "A");
//    vcf.add_record("chrom1", 46, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "CTT", "A");
//    // likelihoods too
//    vcf.add_record("chrom1", 76, "CTT", "A");
//    vcf.add_record("chrom1", 76, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 76, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 76, "CTT", "A");
//    vcf.records[4]->sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
//    vcf.records[5]->sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
//    vcf.records[4]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-50, -3};
//    vcf.records[5]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-50, -16};
//    vcf.records[4]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"] = {47};
//    vcf.records[5]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"] = {56};
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {2, 30};
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {2, 30};
//    vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {2, 30};
//    vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {2, 30};
//    vcf.records[4]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GAPS"] = {4, 0};
//    vcf.records[5]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GAPS"] = {4, 1};
//    // incompatible
//    vcf.add_record("chrom1", 85, "A", "G");
//    vcf.add_record("chrom1", 85, "T", "C");
//
//    vcf = vcf.merge_multi_allelic();
//    std::vector<std::string> formats = {"GT", "LIKELIHOOD", "GT_CONF", "MEAN_FWD_COVG", "MEAN_REV_COVG", "GAPS"};
//    vcf.add_formats(formats);
//
//    EXPECT_EQ((uint) 5, vcf.records.size());
//    EXPECT_EQ((uint) 5, vcf.records[0]->pos);
//    EXPECT_EQ((uint) 2, vcf.records[0]->alts.size());
//    EXPECT_EQ((uint) 1, vcf.records[0]->sampleIndex_to_format_to_sampleInfo.size());
//    EXPECT_EQ((uint) 0, vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].size());
//
//    EXPECT_EQ((uint) 46, vcf.records[1]->pos);
//    EXPECT_EQ((uint) 2, vcf.records[1]->alts.size());
//    EXPECT_EQ((uint) 1, vcf.records[1]->sampleIndex_to_format_to_sampleInfo.size());
//    bool found_gt = vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_TRUE(found_gt);
//    EXPECT_EQ((uint) 0, vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["GT"].size());
//
//    EXPECT_EQ((uint) 76, vcf.records[2]->pos);
//    EXPECT_EQ((uint) 2, vcf.records[2]->alts.size());
//    EXPECT_EQ((uint) 1, vcf.records[2]->sampleIndex_to_format_to_sampleInfo.size());
//    found_gt = vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_TRUE(found_gt);
//    EXPECT_EQ((uint) 1, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint) 3, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0].size());
//    bool found = vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0].find("LIKELIHOOD") != vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0].end();
//    EXPECT_TRUE(found);
//    EXPECT_EQ((uint) 3, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"].size());
//    EXPECT_EQ(-50.0, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][0]);
//    EXPECT_EQ(-3.0, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][1]);
//    EXPECT_EQ(-16.0, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][2]);
//    EXPECT_EQ((uint) 3, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GAPS"].size());
//    EXPECT_EQ(4, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GAPS"][0]);
//    EXPECT_EQ(0, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GAPS"][1]);
//    EXPECT_EQ(1, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GAPS"][2]);
//    found = vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0].find("GT_CONF") != vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0].end();
//    EXPECT_TRUE(found);
//    EXPECT_EQ((uint) 1, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"].size());
//    EXPECT_EQ(13.0, vcf.records[2]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"][0]);
//
//    EXPECT_EQ((uint) 85, vcf.records[3]->pos);
//    EXPECT_EQ((uint) 1, vcf.records[3]->alts.size());
//    EXPECT_EQ((uint) 85, vcf.records[4]->pos);
//    EXPECT_EQ((uint) 1, vcf.records[4]->alts.size());
//}
//
//
//
//TEST(VCFTest, merge_multi_allelic___vcf_with_two_samples_and_two_records_second_record_does_not_map_to_second_sample) {
//    // declares vcf with two samples
//    VCF vcf;
//    vcf.add_samples({"sample1", "sample2"});
//
//    // add two bi-allelic records to be merged
//    vcf.add_record("chrom1", 5, "A", "C", "SVTYPE=SNP", "GRAPHTYPE=SIMPLE");
//    vcf.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP", "GRAPHTYPE=SIMPLE");
//
//    // first record maps to both samples
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 2};
//    vcf.records[0]->sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"] = {1, 3};
//
//    // second record maps to first sample only
//    vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 4};
//
//    // do the merge
//    vcf = vcf.merge_multi_allelic();
//
//    // output the vcf just for us to look at it
//    std::vector<std::string> formats = {"MEAN_FWD_COVG"};
//    vcf.add_formats(formats);
//
//    std::cout << vcf.to_string() << std::endl;
//
//
//    /*
//Doubts come here:
//Current output:
//chrom1	6	.	A	G	.	.	SVTYPE=SNP	GT:MEAN_FWD_COVG	.:.	.:1,3
//chrom1	6	.	A	G,C	.	.	SVTYPE=SNP	GT:MEAN_FWD_COVG	.:1,2,4	.:1,3
//
//Bugs (I think):
//1. First record should be removed, as it got merged into the second one
//2. ".:1,3" -> ".:1,3,0" - we should add a coverage of 0 if the record does not map to a sample
//     * Or should it be ".:1,3" -> ".:1,3,." : a '.' instead of 0?
//
//Expected output (? - could you please fill - it can be the str representation of the record?):
//     */
//}
//
//TEST(VCFTest, correct_dot_alleles) {
//    VCF vcf;
//    // at start
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 0, ".", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 0, "T", ".");
//    // in middle
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 35, ".", "A");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 35, "TA", ".");
//    // multiple alts
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 44, "TA", "T");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 44, "TA", ".");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 44, ".", "T");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom2", 44, ".", "TA");
//
//    string vcf_ref = "TATATGTGTC"
//            "GCGACACTGC"
//            "ATGCATGCAT"
//            "AGTCCTAAAG"
//            "TCCTTAAACG"
//            "TTTATAGTCG";
//
//    vcf.correct_dot_alleles(vcf_ref, "chrom1");
//    vcf.correct_dot_alleles(vcf_ref, "chrom2");
//
//    EXPECT_EQ(vcf.records[0]->ref, "T");
//    EXPECT_EQ(vcf.records[1]->ref, "C");
//    EXPECT_EQ(vcf.records[2]->ref, "TTA");
//    EXPECT_EQ(vcf.records[3]->ref, "TA");
//    EXPECT_EQ(vcf.records[4]->ref, "TA");
//    EXPECT_EQ(vcf.records[5]->ref, "CTA");
//    EXPECT_EQ(vcf.records[6]->ref, "T");
//    EXPECT_EQ(vcf.records[7]->ref, "T");
//
//    EXPECT_EQ(vcf.records[0]->alts.size(), 1);
//    EXPECT_EQ(vcf.records[0]->alts[0], "TAT");
//    EXPECT_EQ(vcf.records[1]->alts.size(), 1);
//    EXPECT_EQ(vcf.records[1]->alts[0], "CA");
//    EXPECT_EQ(vcf.records[2]->alts.size(), 1);
//    EXPECT_EQ(vcf.records[2]->alts[0], "T");
//    EXPECT_EQ(vcf.records[3]->alts.size(), 1);
//    EXPECT_EQ(vcf.records[3]->alts[0], "T");
//    EXPECT_EQ(vcf.records[4]->alts.size(), 1);
//    EXPECT_EQ(vcf.records[4]->alts[0], "A");
//    EXPECT_EQ(vcf.records[5]->alts.size(), 1);
//    EXPECT_EQ(vcf.records[5]->alts[0], "C");
//    EXPECT_EQ(vcf.records[6]->alts.size(), 1);
//    EXPECT_EQ(vcf.records[6]->alts[0], "TT");
//    EXPECT_EQ(vcf.records[7]->alts.size(), 1);
//    EXPECT_EQ(vcf.records[7]->alts[0], "TTA");
//}
//
//TEST(VCFTest, make_gt_compatible) {
//    VCF vcf;
//    // no gt
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 5, "A", "C");
//    // gt incompatible no likelihoods
//    vcf.add_record("chrom1", 46, "CTT", "A");
//    vcf.add_record("chrom1", 46, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 46, "CTT", "A");
//    // gt incompatible, likelihoods too both alts
//    vcf.add_record("chrom1", 76, "CTT", "A");
//    vcf.add_record("chrom1", 76, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 76, "CTT", "TA");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 76, "CTT", "A");
//    vcf.records[4]->sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
//    vcf.records[5]->sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
//    vcf.records[4]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-50, -3};
//    vcf.records[5]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-50, -16};
//    vcf.records[4]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"] = {47};
//    vcf.records[5]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"] = {56};
//    // gt incompatible one ref, ref correct
//    vcf.add_record("chrom1", 85, "A", "G");
//    vcf.add_record("chrom1", 85, "A", "C");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 85, "A", "A");
//    vcf.records[6]->sampleIndex_to_format_to_sampleInfo[0]["GT"] = {1};
//    vcf.records[6]->sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
//    vcf.records[7]->sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
//    vcf.records[6]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-5, -30};
//    vcf.records[7]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-5, -16};
//    vcf.records[6]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"] = {47};
//    vcf.records[7]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"] = {56};
//    // gt incompatible one ref, ref wrong
//    vcf.add_record("chrom1", 95, "A", "G");
//    vcf.add_record("chrom1", 95, "A", "C");
//    vcf.add_a_new_record_discovered_in_a_sample_and_genotype_it("sample", "chrom1", 95, "A", "A");
//    vcf.records[8]->sampleIndex_to_format_to_sampleInfo[0]["GT"] = {1};
//    vcf.records[8]->sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
//    vcf.records[9]->sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
//    vcf.records[8]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-50, -3};
//    vcf.records[9]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-50, -60};
//    vcf.records[8]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"] = {47};
//    vcf.records[9]->sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"] = {10};
//
//    vcf.make_gt_compatible();
//
//    bool found_gt = vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[0]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    found_gt = vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].find("GT") != vcf.records[1]->sampleIndex_to_format_to_sampleInfo[0].end();
//    EXPECT_FALSE(found_gt);
//    EXPECT_EQ((uint) 0, vcf.records[2]->sampleIndex_to_format_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint) 0, vcf.records[3]->sampleIndex_to_format_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint16_t) 1, vcf.records[4]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint) 0, vcf.records[5]->sampleIndex_to_format_to_sampleInfo[0]["GT"].size());
//    EXPECT_EQ((uint16_t) 0, vcf.records[6]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 0, vcf.records[7]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 1, vcf.records[8]->sampleIndex_to_format_to_sampleInfo[0]["GT"][0]);
//    EXPECT_EQ((uint16_t) 0, vcf.records[9]->sampleIndex_to_format_to_sampleInfo[0]["GT"].size());
//}
//
//TEST(VCFTest, equals) {
//    VCF vcf;
//    vcf.add_record("chrom1", 5, "A", "G");
//    vcf.add_record("chrom1", 46, "T", "TA");
//    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
//    std::vector<std::string> empty = {};
//    vcf.add_record(vr, empty);
//    EXPECT_EQ(vcf, vcf);
//
//    // different order
//    VCF vcf1;
//    vcf1.add_record("chrom1", 5, "A", "G");
//    vcf1.add_record(vr, empty);
//    vcf1.add_record("chrom1", 46, "T", "TA");
//    EXPECT_EQ(vcf1, vcf1);
//    EXPECT_EQ(vcf, vcf1);
//    EXPECT_EQ(vcf1, vcf);
//
//    // same length, one different
//    VCF vcf2;
//    vcf2.add_record("chrom1", 10, "A", "G");
//    vcf2.add_record(vr, empty);
//    vcf2.add_record("chrom1", 46, "T", "TA");
//    EXPECT_EQ(vcf2, vcf2);
//    EXPECT_EQ((vcf == vcf2), false);
//    EXPECT_EQ((vcf2 == vcf), false);
//
//    // different length
//    VCF vcf3;
//    vcf3.add_record("chrom1", 5, "A", "G");
//    vcf3.add_record(vr, empty);
//    vcf3.add_record("chrom1", 46, "T", "TA");
//    vcf3.add_record("chrom1", 30, "G", "CC");
//    EXPECT_EQ(vcf3, vcf3);
//    EXPECT_EQ((vcf == vcf3), false);
//    EXPECT_EQ((vcf3 == vcf), false);
//}
//
//class VCFTest___serialization___Fixture : public ::testing::Test {
//protected:
//    VCF vcf_with_zero_records;
//    VCF vcf_with_one_record;
//    VCF vcf_with_three_records;
//    void SetUp() override {
//        {
//            vcf_with_one_record.add_record("chrom1", 5, "A", "G", "GRAPHTYPE=SIMPLE;SVTYPE=SNP");
//        }
//
//        {
//            vcf_with_three_records.add_record("chrom1", 5, "A", "G", "GRAPHTYPE=SIMPLE;SVTYPE=SNP");
//            vcf_with_three_records.add_record("chrom1", 46, "T", "TA", "GRAPHTYPE=SIMPLE;SVTYPE=SNP");
//            VCFRecord vcf_record = VCFRecord("chrom1", 79, "C", "G", "GRAPHTYPE=SIMPLE;SVTYPE=SNP");
//            std::vector<std::string> empty_sample_names = {};
//            vcf_with_three_records.add_record(vcf_record, empty_sample_names);
//        }
//    }
//
//    void TearDown() override {
//    }
//};
//
//TEST_F(VCFTest___serialization___Fixture, save_vcf_with_zero_records___load_vcf___expect_equal_vcf) {
//    vcf_with_zero_records.save("vcf_serialization_test_zero.vcf");
//
//    // TODO: use factory pattern instead
//    VCF actual;
//    actual.load("vcf_serialization_test_zero.vcf");
//
//    VCF& expected = vcf_with_zero_records;
//    EXPECT_EQ(actual, expected);
//}
//
//TEST_F(VCFTest___serialization___Fixture, save_vcf_with_one_record___load_vcf___expect_equal_vcf) {
//    vcf_with_one_record.save("vcf_serialization_test_one.vcf");
//
//    // TODO: use factory pattern instead
//    VCF actual;
//    actual.load("vcf_serialization_test_one.vcf");
//
//    VCF& expected = vcf_with_one_record;
//    EXPECT_EQ(actual, expected);
//}
//
//
//TEST_F(VCFTest___serialization___Fixture, save_vcf_with_three_records___load_vcf___expect_equal_vcf) {
//    vcf_with_three_records.save("vcf_serialization_test_three.vcf");
//
//    // TODO: use factory pattern instead
//    VCF actual;
//    actual.load("vcf_serialization_test_three.vcf");
//
//    VCF& expected = vcf_with_three_records;
//    EXPECT_EQ(actual, expected);
//}
//
//
//class VCFTest___to_string___Fixture : public ::testing::Test {
//protected:
//    class VCF_DummyHeader_Mock : public VCF {
//    public:
//        using VCF::VCF;
//        const std::string dummy_header {"##Dummy_header;\n"};
//        virtual std::string header() const {
//            return dummy_header;
//        }
//    };
//
//    std::shared_ptr<VCF> vcf_with_all_records;
//    VCFRecord graph_type_is_simple_sv_is_snp;
//    VCFRecord graph_type_is_nested_sv_is_snp;
//    VCFRecord graph_type_has_too_many_alts_sv_is_snp;
//    VCFRecord graph_type_is_simple_sv_is_indel;
//    VCFRecord graph_type_is_simple_sv_is_ph_snps;
//    VCFRecord graph_type_is_simple_sv_is_complex;
//    VCFRecord record_with_dot_allele;
//    std::vector<std::string> sample_names;
//    void SetUp() override {
//        graph_type_is_simple_sv_is_snp = VCFRecord("0", 0, "0", "0", "SVTYPE=SNP", "GRAPHTYPE=SIMPLE");
//        graph_type_is_nested_sv_is_snp = VCFRecord("0", 1, "0", "0", "SVTYPE=SNP", "GRAPHTYPE=NESTED");
//        graph_type_has_too_many_alts_sv_is_snp = VCFRecord("0", 2, "0", "0", "SVTYPE=SNP", "GRAPHTYPE=TOO_MANY_ALTS");
//        graph_type_is_simple_sv_is_indel = VCFRecord("0", 3, "0", "0", "SVTYPE=INDEL", "GRAPHTYPE=SIMPLE");
//        graph_type_is_simple_sv_is_ph_snps = VCFRecord("0", 4, "0", "0", "SVTYPE=PH_SNPs", "GRAPHTYPE=SIMPLE");
//        graph_type_is_simple_sv_is_complex = VCFRecord("0", 5, "0", "0", "SVTYPE=COMPLEX", "GRAPHTYPE=SIMPLE");
//        vcf_with_all_records = std::make_shared<VCF_DummyHeader_Mock>();
//        record_with_dot_allele = VCFRecord("0", 6, ".", ".", ".", ".");
//        vcf_with_all_records->add_record(graph_type_is_simple_sv_is_snp, sample_names);
//        vcf_with_all_records->add_record(graph_type_is_nested_sv_is_snp, sample_names);
//        vcf_with_all_records->add_record(graph_type_has_too_many_alts_sv_is_snp, sample_names);
//        vcf_with_all_records->add_record(graph_type_is_simple_sv_is_indel, sample_names);
//        vcf_with_all_records->add_record(graph_type_is_simple_sv_is_ph_snps, sample_names);
//        vcf_with_all_records->add_record(graph_type_is_simple_sv_is_complex, sample_names);
//        vcf_with_all_records->add_record(record_with_dot_allele, sample_names);
//    }
//
//    void TearDown() override {
//    }
//};
//
//TEST_F(VCFTest___to_string___Fixture, graph_type_is_simple_sv_is_snp) {
//    std::string actual = vcf_with_all_records->to_string(false, true, false, false, true, false, false, false);
//
//    std::string expected = "##Dummy_header;\n0\t1\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE=SIMPLE\tGT\n";
//
//    EXPECT_EQ(actual, expected);
//}
//
//TEST_F(VCFTest___to_string___Fixture, graph_type_is_nested_sv_is_snp) {
//    std::string actual = vcf_with_all_records->to_string(false, false, true, false, true, false, false, false);
//
//    std::string expected = "##Dummy_header;\n0\t2\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE=NESTED\tGT\n";
//
//    EXPECT_EQ(actual, expected);
//}
//
//TEST_F(VCFTest___to_string___Fixture, graph_type_has_too_many_alts_sv_is_snp) {
//    std::string actual = vcf_with_all_records->to_string(false, false, false, true, true, false, false, false);
//
//    std::string expected = "##Dummy_header;\n0\t3\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE=TOO_MANY_ALTS\tGT\n";
//
//    EXPECT_EQ(actual, expected);
//}
//
//TEST_F(VCFTest___to_string___Fixture, graph_type_is_simple_sv_is_indel) {
//    std::string actual = vcf_with_all_records->to_string(false, true, false, false, false, true, false, false);
//
//    std::string expected = "##Dummy_header;\n0\t4\t.\t0\t0\t.\t.\tSVTYPE=INDEL;GRAPHTYPE=SIMPLE\tGT\n";
//
//    EXPECT_EQ(actual, expected);
//}
//
//TEST_F(VCFTest___to_string___Fixture, graph_type_is_simple_sv_is_ph_snps) {
//    std::string actual = vcf_with_all_records->to_string(false, true, false, false, false, false, true, false);
//
//    std::string expected = "##Dummy_header;\n0\t5\t.\t0\t0\t.\t.\tSVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE\tGT\n";
//
//    EXPECT_EQ(actual, expected);
//}
//
//TEST_F(VCFTest___to_string___Fixture, graph_type_is_simple_sv_is_complex) {
//    std::string actual = vcf_with_all_records->to_string(false, true, false, false, false, false, false, true);
//
//    std::string expected = "##Dummy_header;\n0\t6\t.\t0\t0\t.\t.\tSVTYPE=COMPLEX;GRAPHTYPE=SIMPLE\tGT\n";
//
//    EXPECT_EQ(actual, expected);
//}
//
//
//TEST_F(VCFTest___to_string___Fixture, record_with_dot_allele) {
//    std::string actual = vcf_with_all_records->to_string(true, false, false, false, false, false, false, false);
//
//    std::string expected = "##Dummy_header;\n0\t7\t.\t.\t.\t.\t.\t.;.\tGT\n";
//
//    EXPECT_EQ(actual, expected);
//}
//
//
//TEST_F(VCFTest___to_string___Fixture, all_records_filtered_out) {
//    std::string actual = vcf_with_all_records->to_string(false, false, false, false, false, false, false, false);
//
//    std::string expected = "##Dummy_header;\n";
//
//    EXPECT_EQ(actual, expected);
//}
//
//
//TEST_F(VCFTest___to_string___Fixture, no_records_filtered_out) {
//    std::string actual = vcf_with_all_records->to_string(true, true, true, true, true, true, true, true);
//
//    std::string expected = "##Dummy_header;\n";
//    expected += "0\t1\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE=SIMPLE\tGT\n";
//    expected += "0\t2\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE=NESTED\tGT\n";
//    expected += "0\t3\t.\t0\t0\t.\t.\tSVTYPE=SNP;GRAPHTYPE=TOO_MANY_ALTS\tGT\n";
//    expected += "0\t4\t.\t0\t0\t.\t.\tSVTYPE=INDEL;GRAPHTYPE=SIMPLE\tGT\n";
//    expected += "0\t5\t.\t0\t0\t.\t.\tSVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE\tGT\n";
//    expected += "0\t6\t.\t0\t0\t.\t.\tSVTYPE=COMPLEX;GRAPHTYPE=SIMPLE\tGT\n";
//    expected += "0\t7\t.\t.\t.\t.\t.\t.;.\tGT\n";
//
//    EXPECT_EQ(actual, expected);
//}
//
//
//TEST(VCFTest, filter) {
//    VCF vcf, vcf1, vcf2, vcf3, vcf4;
//    vcf.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
//    vcf.add_record("chrom1", 46, "T", "TA", "SVTYPE=INDEL;GRAPHTYPE=NESTED");
//    vcf.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE");
//    vcf.add_record("chrom1", 79, "CTT", "ATA", "SVTYPE=PH_SNPs;GRAPHTYPE=NESTED");
//    vcf.save("vcf_filter_test.vcf", false, true, false, false, true, false, true, false);
//
//    vcf1.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
//    vcf1.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE");
//    vcf2.load("vcf_filter_test.vcf");
//    EXPECT_EQ(vcf2, vcf1);
//
//    vcf.save("vcf_filter_test.vcf", false, true, true, false, false, false, true, false);
//    vcf3.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE");
//    vcf3.add_record("chrom1", 79, "CTT", "ATA", "SVTYPE=PH_SNPs;GRAPHTYPE=NESTED");
//    vcf4.load("vcf_filter_test.vcf");
//    EXPECT_EQ(vcf3, vcf4);
//}