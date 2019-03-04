#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "vcf.h"
#include "vcfrecord.h"
#include "interval.h"
#include "localnode.h"
#include <stdint.h>
#include <iostream>


using namespace std;

TEST(VCFTest, add_record_with_values) {

    VCF vcf;
    EXPECT_EQ((uint) 0, vcf.records.size());
    vcf.add_record("chrom1", 5, "A", "G");
    EXPECT_EQ((uint) 1, vcf.records.size());
}

TEST(VCFTest, add_record_twice_with_values) {

    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 5, "A", "G");
    EXPECT_EQ((uint) 1, vcf.records.size());
}

TEST(VCFTest, add_two_records_with_values) {

    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    EXPECT_EQ((uint) 2, vcf.records.size());
}

TEST(VCFTest, add_two_records_and_a_repeat_with_values) {

    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 5, "A", "G");
    EXPECT_EQ((uint) 2, vcf.records.size());
}

TEST(VCFTest, add_record_by_record) {
    VCF vcf;
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_record(vr, empty);
    EXPECT_EQ((uint) 1, vcf.records.size());
}

TEST(VCFTest, add_record_by_record_and_values) {
    VCF vcf;
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_record(vr, empty);
    vcf.add_record("chrom1", 79, "C", "G");
    EXPECT_EQ((uint) 1, vcf.records.size());
}

TEST(VCFTest, add_record_by_values_and_record) {
    VCF vcf;
    vcf.add_record("chrom1", 79, "C", "G");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_record(vr, empty);
    EXPECT_EQ((uint) 1, vcf.records.size());
}

TEST(VCFTest, add_record_by_record_returned_by_reference) {
    VCF vcf;
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    VCFRecord &ref_vr = vcf.add_record(vr, empty);
    EXPECT_EQ(ref_vr.chrom, "chrom1");
    EXPECT_EQ(ref_vr.pos, (uint) 79);
}

TEST(VCFTest, add_samples_empty) {
    VCF vcf;
    std::vector<std::string> samples;
    vcf.add_samples(samples);
    EXPECT_EQ(vcf.samples.size(), (uint)0);
    EXPECT_EQ(vcf.records.size(), (uint)0);
}

TEST(VCFTest, add_samples_simple) {
    VCF vcf;
    std::vector<std::string> samples = {"hello", "there", "people"};
    vcf.add_samples(samples);
    EXPECT_ITERABLE_EQ(std::vector<std::string>, samples, vcf.samples);
    EXPECT_EQ(vcf.records.size(), (uint)0);
}

TEST(VCFTest, add_samples_with_record) {
    VCF vcf;
    vcf.add_sample_gt("sample", "chrom1", 5, "A", "G");

    std::vector<std::string> samples = {"hello", "there", "people"};
    vcf.add_samples(samples);

    std::vector<std::string> exp_samples = {"sample", "hello", "there", "people"};

    EXPECT_ITERABLE_EQ(std::vector<std::string>, exp_samples, vcf.samples);
    EXPECT_EQ(vcf.records.size(), (uint)1);
    EXPECT_EQ(vcf.records[0].samples.size(), exp_samples.size());
}

TEST(VCFTest, add_sample_gt) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    vcf.add_sample_gt("sample", "chrom1", 46, "T", "TA");
    uint j = 1;
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.records[1].samples.size());
    EXPECT_EQ((uint16_t) 1, vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[0].samples.size());
    EXPECT_TRUE(vcf.records[0].samples[0].find("GT") == vcf.records[0].samples[0].end());
    EXPECT_EQ(j, vcf.records[2].samples.size());
    EXPECT_TRUE(vcf.records[2].samples[0].find("GT") == vcf.records[2].samples[0].end());
    EXPECT_EQ(j, vcf.records[3].samples.size());
    EXPECT_TRUE(vcf.records[3].samples[0].find("GT") == vcf.records[3].samples[0].end());

    vcf.add_sample_gt("sample", "chrom1", 79, "C", "C");
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.records[1].samples.size());
    EXPECT_EQ((uint16_t) 1, vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[0].samples.size());
    EXPECT_TRUE(vcf.records[0].samples[0].find("GT") == vcf.records[0].samples[0].end());
    EXPECT_EQ(j, vcf.records[2].samples.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[2].samples[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[3].samples.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[3].samples[0]["GT"][0]);
}

TEST(VCFTest, add_record_by_record_with_existing_sample) {
    VCF vcf;
    vcf.add_sample_gt("sample", "chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    VCFRecord &ref_vr = vcf.add_record(vr, empty);
    EXPECT_EQ(ref_vr.chrom, "chrom1");
    EXPECT_EQ(ref_vr.pos, (uint) 79);
    EXPECT_EQ(ref_vr.samples.size(), (uint) 1);
}

TEST(VCFTest, add_record_by_record_with_same_existing_sample) {
    VCF vcf;
    vcf.add_sample_gt("sample", "chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::unordered_map<std::string, std::vector<uint16_t>> empty_map;
    vr.samples.emplace_back(empty_map);
    vr.samples[0]["GT"] = {1};
    std::vector<std::string> samples = {"sample"};
    VCFRecord &ref_vr = vcf.add_record(vr, samples);
    EXPECT_EQ(ref_vr.chrom, "chrom1");
    EXPECT_EQ(ref_vr.pos, (uint) 79);
    EXPECT_EQ(ref_vr.samples.size(), (uint) 1);
    EXPECT_ITERABLE_EQ(vector<std::string>, samples, vcf.samples);
    EXPECT_FALSE(ref_vr.samples[0].find("GT") == ref_vr.samples[0].end());
    EXPECT_EQ(ref_vr.samples[0]["GT"][0], (uint16_t)1);
}

TEST(VCFTest, add_record_by_record_with_different_existing_sample) {
    VCF vcf;
    vcf.add_sample_gt("sample", "chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::unordered_map<std::string, std::vector<uint16_t>> empty_map;
    vr.samples.emplace_back(empty_map);
    vr.samples[0]["GT"] = {1};
    std::vector<std::string> samples = {"sample1"};
    VCFRecord &ref_vr = vcf.add_record(vr, samples);
    EXPECT_EQ(ref_vr.chrom, "chrom1");
    EXPECT_EQ(ref_vr.pos, (uint) 79);
    EXPECT_EQ(ref_vr.samples.size(), (uint) 2);
    std::vector<std::string> exp_samples = {"sample", "sample1"};
    EXPECT_ITERABLE_EQ(vector<std::string>, exp_samples, vcf.samples);
    EXPECT_TRUE(ref_vr.samples[0].find("GT") == ref_vr.samples[0].end());
    EXPECT_FALSE(ref_vr.samples[1].find("GT") == ref_vr.samples[0].end());
    EXPECT_EQ(ref_vr.samples[1]["GT"][0], (uint16_t)1);
}

TEST(VCFTest, add_sample_ref_alleles) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_record("chrom2", 30, "C", "A");

    vcf.add_sample_ref_alleles("sample", "chrom1", 15, 78);
    EXPECT_EQ((uint) 1, vcf.samples.size());
    EXPECT_EQ((uint) 5, vcf.records.size());
    EXPECT_EQ((uint) 1, vcf.records[0].samples.size());
    EXPECT_TRUE(vcf.records[0].samples[0].find("GT") == vcf.records[0].samples[0].end());
    EXPECT_EQ((uint) 1, vcf.records[1].samples.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ((uint) 1, vcf.records[2].samples.size());
    EXPECT_TRUE(vcf.records[2].samples[0].find("GT") == vcf.records[2].samples[0].end());
    EXPECT_EQ((uint) 1, vcf.records[3].samples.size());
    EXPECT_TRUE(vcf.records[3].samples[0].find("GT") == vcf.records[3].samples[0].end());
    EXPECT_EQ((uint) 1, vcf.records[4].samples.size());
    EXPECT_TRUE(vcf.records[4].samples[0].find("GT") == vcf.records[4].samples[0].end());

    vcf.add_sample_ref_alleles("sample2", "chrom1", 5, 46);
    EXPECT_EQ((uint) 2, vcf.samples.size());
    EXPECT_EQ((uint) 5, vcf.records.size());
    EXPECT_EQ((uint) 2, vcf.records[0].samples.size());
    EXPECT_EQ((uint16_t) 0, vcf.records[0].samples[1]["GT"][0]);
    EXPECT_EQ((uint) 2, vcf.records[1].samples.size());
    EXPECT_TRUE(vcf.records[1].samples[1].find("GT") == vcf.records[1].samples[1].end());
    EXPECT_EQ((uint) 2, vcf.records[2].samples.size());
    EXPECT_TRUE(vcf.records[2].samples[1].find("GT") == vcf.records[2].samples[1].end());
    EXPECT_EQ((uint) 2, vcf.records[3].samples.size());
    EXPECT_TRUE(vcf.records[3].samples[1].find("GT") == vcf.records[3].samples[1].end());
    EXPECT_EQ((uint) 2, vcf.records[4].samples.size());
    EXPECT_TRUE(vcf.records[4].samples[1].find("GT") == vcf.records[4].samples[1].end());
}

TEST(VCFTest, reorder_add_record_and_sample) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_sample_gt("sample1", "chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_sample_gt("sample2", "chrom1", 79, "C", "C");
    vcf.add_sample_gt("sample1", "chrom1", 79, "C", "A");

    vcf.sort_records();

    EXPECT_EQ((uint) 2, vcf.samples.size());
    EXPECT_EQ((uint) 4, vcf.records.size());
    EXPECT_EQ((uint) 2, vcf.records[0].samples.size());
    EXPECT_EQ((uint) 2, vcf.records[1].samples.size());
    EXPECT_EQ((uint) 2, vcf.records[2].samples.size());
    EXPECT_EQ((uint) 2, vcf.records[3].samples.size());
    EXPECT_TRUE(vcf.records[0].samples[0].find("GT") == vcf.records[0].samples[0].end());
    EXPECT_EQ((uint16_t) 1, vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 1, vcf.records[2].samples[0]["GT"][0]);
    EXPECT_TRUE(vcf.records[3].samples[0].find("GT") == vcf.records[3].samples[0].end());
    EXPECT_TRUE(vcf.records[0].samples[1].find("GT") == vcf.records[0].samples[1].end());
    EXPECT_TRUE(vcf.records[1].samples[1].find("GT") == vcf.records[1].samples[1].end());
    EXPECT_EQ((uint16_t) 0, vcf.records[2].samples[1]["GT"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[3].samples[1]["GT"][0]);

}


TEST(VCFTest, clear) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_record(vr, empty);
    uint j = 3;
    EXPECT_EQ(j, vcf.records.size());

    vcf.clear();
    j = 0;
    EXPECT_EQ(j, vcf.records.size());
}

TEST(VCFTest, append_vcf_simple_case) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    VCF new_vcf;
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom2", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_record("chrom2", 79, "C", "A");

    vcf.append_vcf(new_vcf);
    EXPECT_EQ((uint) 8, vcf.records.size());
    for (uint i = 0; i < 4; ++i) {
        EXPECT_EQ(vcf.records[i].chrom, "chrom1");
    }
    for (uint i = 4; i < 8; ++i) {
        EXPECT_EQ(vcf.records[i].chrom, "chrom2");
    }
    EXPECT_EQ((uint) 5, vcf.records[4].pos);
    EXPECT_EQ("TA", vcf.records[5].alt[0]);
    EXPECT_EQ((uint) 79, vcf.records[6].pos);
    EXPECT_EQ("A", vcf.records[7].alt[0]);
}

TEST(VCFTest, append_vcf_some_duplicate_records) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    VCF new_vcf;
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom1", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_record("chrom1", 79, "C", "A");

    vcf.append_vcf(new_vcf);
    EXPECT_EQ((uint) 6, vcf.records.size());
    for (uint i = 0; i < 4; ++i) {
        EXPECT_EQ(vcf.records[i].chrom, "chrom1");
    }
    for (uint i = 4; i < 6; ++i) {
        EXPECT_EQ(vcf.records[i].chrom, "chrom2");
    }
    EXPECT_EQ((uint) 5, vcf.records[4].pos);
    EXPECT_EQ((uint) 79, vcf.records[5].pos);
}

TEST(VCFTest, append_vcf_one_sample) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_sample_gt("sample", "chrom1", 79, "C", "G");

    VCF new_vcf;
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom1", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_record("chrom1", 79, "C", "A");

    vcf.append_vcf(new_vcf);
    EXPECT_EQ((uint) 1, vcf.samples.size());
    EXPECT_EQ("sample", vcf.samples[0]);
    EXPECT_EQ((uint) 1, vcf.records[0].samples.size());
    EXPECT_EQ((uint) 1, vcf.records[5].samples.size());
    bool found_gt = vcf.records[2].samples[0].find("GT") != vcf.records[2].samples[0].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint) 1, vcf.records[2].samples[0]["GT"][0]);
    found_gt = vcf.records[0].samples[0].find("GT") != vcf.records[0].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[4].samples[0].find("GT") != vcf.records[4].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[3].samples[0].find("GT") != vcf.records[3].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[5].samples[0].find("GT") != vcf.records[5].samples[0].end();
    EXPECT_FALSE(found_gt);
}

TEST(VCFTest, append_vcf_one_sample_in_new_vcf) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    VCF new_vcf;
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom1", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_record("chrom1", 79, "C", "A");
    new_vcf.add_sample_gt("sample", "chrom2", 5, "A", "G");

    vcf.append_vcf(new_vcf);
    EXPECT_EQ((uint) 1, vcf.samples.size());
    EXPECT_EQ("sample", vcf.samples[0]);
    EXPECT_EQ((uint) 1, vcf.records[0].samples.size());
    EXPECT_EQ((uint) 1, vcf.records[5].samples.size());
    bool found_gt = vcf.records[4].samples[0].find("GT") != vcf.records[4].samples[0].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint) 1, vcf.records[4].samples[0]["GT"][0]);
    found_gt = vcf.records[0].samples[0].find("GT") != vcf.records[0].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[2].samples[0].find("GT") != vcf.records[2].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[3].samples[0].find("GT") != vcf.records[3].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[5].samples[0].find("GT") != vcf.records[5].samples[0].end();
    EXPECT_FALSE(found_gt);
}

TEST(VCFTest, append_vcf_shared_sample) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_sample_gt("sample", "chrom1", 46, "T", "TA");

    VCF new_vcf;
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom1", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_record("chrom1", 79, "C", "A");
    new_vcf.add_sample_gt("sample", "chrom1", 46, "T", "TA");

    vcf.append_vcf(new_vcf);
    EXPECT_EQ((uint) 1, vcf.samples.size());
    EXPECT_EQ("sample", vcf.samples[0]);
    EXPECT_EQ((uint) 1, vcf.records[0].samples.size());
    EXPECT_EQ((uint) 1, vcf.records[5].samples.size());
    bool found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint) 1, vcf.records[1].samples[0]["GT"][0]);
    found_gt = vcf.records[0].samples[0].find("GT") != vcf.records[0].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[4].samples[0].find("GT") != vcf.records[4].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[2].samples[0].find("GT") != vcf.records[2].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[3].samples[0].find("GT") != vcf.records[3].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[5].samples[0].find("GT") != vcf.records[5].samples[0].end();
    EXPECT_FALSE(found_gt);
}

TEST(VCFTest, append_vcf_shared_samples_different_order) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    vcf.add_sample_gt("sample", "chrom1", 46, "T", "TA");

    VCF new_vcf;
    new_vcf.add_record("chrom1", 79, "C", "A");
    new_vcf.add_record("chrom2", 5, "A", "G");
    new_vcf.add_record("chrom1", 46, "T", "TA");
    new_vcf.add_record("chrom2", 79, "C", "G");
    new_vcf.add_sample_gt("sample1", "chrom1", 46, "T", "T");
    new_vcf.add_sample_gt("sample1", "chrom1", 79, "C", "A");

    vcf.append_vcf(new_vcf);

    EXPECT_EQ((uint) 2, vcf.samples.size());
    vector<string> v = {"sample", "sample1"};
    EXPECT_ITERABLE_EQ(vector<string>, v, vcf.samples);
    EXPECT_EQ((uint) 2, vcf.records[0].samples.size());
    EXPECT_EQ((uint) 2, vcf.records[1].samples.size());
    EXPECT_EQ((uint) 2, vcf.records[2].samples.size());
    EXPECT_EQ((uint) 2, vcf.records[3].samples.size());
    EXPECT_EQ((uint) 2, vcf.records[4].samples.size());
    EXPECT_EQ((uint) 2, vcf.records[5].samples.size());

    vector<uint16_t> alt_gt = {1};
    vector<uint16_t> ref_gt = {0};

    EXPECT_FALSE(vcf.records[0].samples[0].find("GT") != vcf.records[0].samples[0].end());
    EXPECT_FALSE(vcf.records[0].samples[1].find("GT") != vcf.records[0].samples[1].end());

    EXPECT_TRUE(vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end());
    EXPECT_TRUE(vcf.records[1].samples[1].find("GT") != vcf.records[1].samples[1].end());
    EXPECT_ITERABLE_EQ(vector<uint16_t>, vcf.records[1].samples[0]["GT"], alt_gt);
    EXPECT_ITERABLE_EQ(vector<uint16_t>, vcf.records[1].samples[1]["GT"], ref_gt);

    EXPECT_FALSE(vcf.records[2].samples[0].find("GT") != vcf.records[2].samples[0].end());
    EXPECT_FALSE(vcf.records[2].samples[1].find("GT") != vcf.records[2].samples[1].end());

    EXPECT_FALSE(vcf.records[3].samples[0].find("GT") != vcf.records[3].samples[0].end());
    EXPECT_TRUE(vcf.records[3].samples[1].find("GT") != vcf.records[3].samples[1].end());
    EXPECT_ITERABLE_EQ(vector<uint16_t>, vcf.records[3].samples[1]["GT"], alt_gt);

    EXPECT_FALSE(vcf.records[4].samples[0].find("GT") != vcf.records[4].samples[0].end());
    EXPECT_FALSE(vcf.records[4].samples[1].find("GT") != vcf.records[4].samples[1].end());

    EXPECT_FALSE(vcf.records[5].samples[0].find("GT") != vcf.records[5].samples[0].end());
    EXPECT_FALSE(vcf.records[5].samples[1].find("GT") != vcf.records[5].samples[1].end());
}

TEST(VCFTest, sort_records) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_sample_gt("sample", "chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_record("chrom2", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom2", 79, "C", "G");
    vcf.sort_records();

    EXPECT_EQ((uint) 6, vcf.records.size());
    for (uint i = 0; i < 4; ++i) {
        EXPECT_EQ("chrom1", vcf.records[i].chrom);
    }
    for (uint i = 4; i < 6; ++i) {
        EXPECT_EQ("chrom2", vcf.records[i].chrom);
    }
    EXPECT_EQ((uint) 5, vcf.records[0].pos);
    EXPECT_EQ((uint) 5, vcf.records[4].pos);
    EXPECT_EQ((uint) 46, vcf.records[1].pos);
    EXPECT_EQ((uint) 79, vcf.records[2].pos);
    EXPECT_EQ((uint) 79, vcf.records[3].pos);
    EXPECT_EQ((uint) 79, vcf.records[5].pos);
    EXPECT_EQ("G", vcf.records[3].alt[0]);
    EXPECT_EQ("G", vcf.records[5].alt[0]);
}

TEST(VCFTest, pos_in_range) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom2", 20, "A", "G");
    vcf.add_record("chrom2", 79, "C", "G");
    //vcf.sort();

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

TEST(VCFTest, genotype) {
    // not a snp site
    // missing count data
    // not confident
    // confident and have right gt
    // confident and have wrong gt
    // one sample needs regenotyping and the other doesn't
    VCF vcf;

    vcf.add_record("chrom2", 79, "C", "G");

    vcf.add_sample_gt("sample", "chrom1", 2, "T", "TA");
    vcf.add_sample_gt("sample", "chrom1", 5, "A", "G");
    vcf.add_sample_gt("sample", "chrom1", 79, "C", "A");
    vcf.add_sample_gt("sample", "chrom2", 20, "A", "G");
    vcf.add_sample_gt("sample", "chrom2", 79, "C", "C");
    vcf.add_sample_gt("sample", "chrom2", 80, "A", "C");

    vcf.add_sample_gt("asample", "chrom1", 2, "T", "TA");
    vcf.add_sample_gt("asample", "chrom1", 5, "A", "A");
    vcf.add_sample_gt("asample", "chrom1", 79, "C", "A");
    vcf.add_sample_gt("asample", "chrom2", 20, "A", "G");
    vcf.add_sample_gt("asample", "chrom2", 79, "C", "C");
    vcf.add_sample_gt("asample", "chrom2", 80, "A", "A");

    vcf.sort_records();
    std::vector<float> f = {0.0, 0.0};

    // record 0, not a snp site
    vcf.records[0].samples[0]["MEAN_FWD_COVG"] = {0, 10};
    vcf.records[0].samples[0]["MEAN_REV_COVG"] = {1, 20};
    vcf.records[0].samples[1]["MEAN_FWD_COVG"] = {1, 15};
    vcf.records[0].samples[1]["MEAN_REV_COVG"] = {2, 24};
    vcf.records[0].set_format(0,"GAPS", f);
    vcf.records[0].set_format(1,"GAPS", f);


    // record 1, different genotypes but both correct
    vcf.records[1].samples[0]["MEAN_FWD_COVG"] = {0, 10};
    vcf.records[1].samples[0]["MEAN_REV_COVG"] = {1, 20};
    vcf.records[1].samples[1]["MEAN_FWD_COVG"] = {10, 1};
    vcf.records[1].samples[1]["MEAN_REV_COVG"] = {21, 2};
    vcf.records[1].set_format(0,"GAPS", f);
    vcf.records[1].set_format(1,"GAPS", f);

    // record 2, same genotypes first correct
    vcf.records[2].samples[0]["MEAN_FWD_COVG"] = {0, 10};
    vcf.records[2].samples[0]["MEAN_REV_COVG"] = {1, 20};
    vcf.records[2].samples[1]["MEAN_FWD_COVG"] = {10, 1};
    vcf.records[2].samples[1]["MEAN_REV_COVG"] = {21, 2};
    vcf.records[2].set_format(0,"GAPS", f);
    vcf.records[2].set_format(1,"GAPS", f);

    // record 3, same genotypes both wrong
    vcf.records[3].samples[0]["MEAN_FWD_COVG"] = {20, 1};
    vcf.records[3].samples[0]["MEAN_REV_COVG"] = {21, 2};
    vcf.records[3].samples[1]["MEAN_FWD_COVG"] = {10, 1};
    vcf.records[3].samples[1]["MEAN_REV_COVG"] = {21, 2};
    vcf.records[3].set_format(0,"GAPS", f);
    vcf.records[3].set_format(1,"GAPS", f);

    // record 4, missing count data for first sample
    vcf.records[4].samples[0]["MEAN_FWD_COVG"] = {0, 10};
    vcf.records[4].samples[0]["MEAN_REV_COVG"] = {20};
    vcf.records[4].samples[1]["MEAN_FWD_COVG"] = {10, 1};
    vcf.records[4].samples[1]["MEAN_REV_COVG"] = {21, 2};
    vcf.records[4].set_format(0,"GAPS", f);
    vcf.records[4].set_format(1,"GAPS", f);

    // record 5, not confident for second sample
    vcf.records[5].samples[0]["MEAN_FWD_COVG"] = {0, 10};
    vcf.records[5].samples[0]["MEAN_REV_COVG"] = {1, 20};
    vcf.records[5].samples[1]["MEAN_FWD_COVG"] = {2, 1};
    vcf.records[5].samples[1]["MEAN_REV_COVG"] = {4, 2};
    vcf.records[5].set_format(0,"GAPS", f);
    vcf.records[5].set_format(1,"GAPS", f);

    vcf.genotype({30, 30}, 0.01, 30, 0, 1, 0, 0, true);

    cout << vcf << endl;

    // not genotyped first record
    EXPECT_EQ((uint16_t) 1, vcf.records[0].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 1, vcf.records[0].samples[1]["GT"][0]);
    bool found_confidence = vcf.records[0].regt_samples[0].find("GT_CONF")!=vcf.records[0].regt_samples[0].end();
    EXPECT_FALSE(found_confidence);
    found_confidence = vcf.records[0].regt_samples[1].find("GT_CONF")!=vcf.records[0].regt_samples[1].end();
    EXPECT_FALSE(found_confidence);

    // both correct
    EXPECT_EQ((uint) 2, vcf.records[1].samples.size());
    bool found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_TRUE(found_gt);
    found_gt = vcf.records[1].samples[1].find("GT") != vcf.records[1].samples[1].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint) 1, vcf.records[1].samples[0]["GT"].size());
    EXPECT_EQ((uint) 1, vcf.records[1].samples[1]["GT"].size());
    EXPECT_EQ((uint16_t) 1, vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1].samples[1]["GT"][0]);
    // first correct
    EXPECT_EQ((uint16_t) 1, vcf.records[2].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[2].samples[1]["GT"][0]);
    // both wrong
    EXPECT_EQ((uint16_t) 0, vcf.records[3].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[3].samples[1]["GT"][0]);
    // first missing data
    EXPECT_EQ((uint) 0, vcf.records[4].samples[0]["GT"].size());
    EXPECT_EQ((uint16_t) 0, vcf.records[4].samples[1]["GT"][0]);
    // second not confident
    EXPECT_EQ((uint16_t) 1, vcf.records[5].samples[0]["GT"][0]);
    EXPECT_EQ((uint) 0, vcf.records[5].samples[1]["GT"].size());
}

TEST(VCFTest, genotype_with_all_sites) {
    // not a snp site
    // missing count data
    // not confident
    // confident and have right gt
    // confident and have wrong gt
    // one sample needs regenotyping and the other doesn't
    VCF vcf;

    vcf.add_record("chrom2", 79, "CC", "GC");

    vcf.add_sample_gt("sample", "chrom1", 2, "T", "TA");
    vcf.add_sample_gt("sample", "chrom1", 5, "AC", "GC");
    vcf.add_sample_gt("sample", "chrom1", 79, "CC", "AC");
    vcf.add_sample_gt("sample", "chrom2", 20, "AC", "GC");
    vcf.add_sample_gt("sample", "chrom2", 79, "CC", "CC");
    vcf.add_sample_gt("sample", "chrom2", 80, "AC", "CC");

    vcf.add_sample_gt("asample", "chrom1", 2, "T", "TA");
    vcf.add_sample_gt("asample", "chrom1", 5, "AC", "AC");
    vcf.add_sample_gt("asample", "chrom1", 79, "CC", "AC");
    vcf.add_sample_gt("asample", "chrom2", 20, "AC", "GC");
    vcf.add_sample_gt("asample", "chrom2", 79, "CC", "CC");
    vcf.add_sample_gt("asample", "chrom2", 80, "AC", "AC");

    vcf.sort_records();
    std::vector<float> f = {0.0, 0.0};

    // record 0, not a snp site
    vcf.records[0].samples[0]["MEAN_FWD_COVG"].push_back(0);
    vcf.records[0].samples[0]["MEAN_REV_COVG"].push_back(1);
    vcf.records[0].samples[0]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[0].samples[0]["MEAN_REV_COVG"].push_back(20);
    vcf.records[0].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[0].samples[1]["MEAN_REV_COVG"].push_back(2);
    vcf.records[0].samples[1]["MEAN_FWD_COVG"].push_back(15);
    vcf.records[0].samples[1]["MEAN_REV_COVG"].push_back(24);
    vcf.records[0].set_format(0,"GAPS", f);
    vcf.records[0].set_format(1,"GAPS", f);

    // record 1, different genotypes but both correct
    vcf.records[1].samples[0]["MEAN_FWD_COVG"].push_back(0);
    vcf.records[1].samples[0]["MEAN_REV_COVG"].push_back(1);
    vcf.records[1].samples[0]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[1].samples[0]["MEAN_REV_COVG"].push_back(20);
    vcf.records[1].samples[1]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[1].samples[1]["MEAN_REV_COVG"].push_back(21);
    vcf.records[1].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[1].samples[1]["MEAN_REV_COVG"].push_back(2);
    vcf.records[1].set_format(0,"GAPS", f);
    vcf.records[1].set_format(1,"GAPS", f);

    // record 2, same genotypes first correct
    vcf.records[2].samples[0]["MEAN_FWD_COVG"].push_back(0);
    vcf.records[2].samples[0]["MEAN_REV_COVG"].push_back(1);
    vcf.records[2].samples[0]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[2].samples[0]["MEAN_REV_COVG"].push_back(20);
    vcf.records[2].samples[1]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[2].samples[1]["MEAN_REV_COVG"].push_back(21);
    vcf.records[2].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[2].samples[1]["MEAN_REV_COVG"].push_back(2);
    vcf.records[2].set_format(0,"GAPS", f);
    vcf.records[2].set_format(1,"GAPS", f);

    // record 3, same genotypes both wrong
    vcf.records[3].samples[0]["MEAN_FWD_COVG"].push_back(20);
    vcf.records[3].samples[0]["MEAN_REV_COVG"].push_back(21);
    vcf.records[3].samples[0]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[3].samples[0]["MEAN_REV_COVG"].push_back(2);
    vcf.records[3].samples[1]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[3].samples[1]["MEAN_REV_COVG"].push_back(21);
    vcf.records[3].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[3].samples[1]["MEAN_REV_COVG"].push_back(2);
    vcf.records[3].set_format(0,"GAPS", f);
    vcf.records[3].set_format(1,"GAPS", f);

    // record 4, missing count data for first sample
    vcf.records[4].samples[0]["MEAN_FWD_COVG"].push_back(0);
    vcf.records[4].samples[0]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[4].samples[0]["MEAN_REV_COVG"].push_back(20);
    vcf.records[4].samples[1]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[4].samples[1]["MEAN_REV_COVG"].push_back(21);
    vcf.records[4].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[4].samples[1]["MEAN_REV_COVG"].push_back(2);
    vcf.records[4].set_format(0,"GAPS", f);
    vcf.records[4].set_format(1,"GAPS", f);

    // record 5, not confident for second sample
    vcf.records[5].samples[0]["MEAN_FWD_COVG"].push_back(0);
    vcf.records[5].samples[0]["MEAN_REV_COVG"].push_back(1);
    vcf.records[5].samples[0]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[5].samples[0]["MEAN_REV_COVG"].push_back(20);
    vcf.records[5].samples[1]["MEAN_FWD_COVG"].push_back(2);
    vcf.records[5].samples[1]["MEAN_REV_COVG"].push_back(4);
    vcf.records[5].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[5].samples[1]["MEAN_REV_COVG"].push_back(2);
    vcf.records[5].set_format(0,"GAPS", f);
    vcf.records[5].set_format(1,"GAPS", f);

    bool snps_only = false;
    vcf.genotype({30, 30}, 0.01, 30, 0, 1, 0, 0, snps_only);

    cout << vcf << endl;

    // first record now genotyped
    EXPECT_EQ((uint16_t) 1, vcf.records[0].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 1, vcf.records[0].samples[1]["GT"][0]);
    bool found_confidence = vcf.records[0].regt_samples[0].find("GT_CONF")!=vcf.records[0].regt_samples[0].end();
    EXPECT_TRUE(found_confidence);
    found_confidence = vcf.records[0].regt_samples[1].find("GT_CONF")!=vcf.records[0].regt_samples[1].end();
    EXPECT_TRUE(found_confidence);
    // both correct
    EXPECT_EQ((uint) 2, vcf.records[1].samples.size());
    bool found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_TRUE(found_gt);
    found_gt = vcf.records[1].samples[1].find("GT") != vcf.records[1].samples[1].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint) 1, vcf.records[1].samples[0]["GT"].size());
    EXPECT_EQ((uint) 1, vcf.records[1].samples[1]["GT"].size());
    EXPECT_EQ((uint16_t) 1, vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[1].samples[1]["GT"][0]);
    // first correct
    EXPECT_EQ((uint16_t) 1, vcf.records[2].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[2].samples[1]["GT"][0]);
    // both wrong
    EXPECT_EQ((uint16_t) 0, vcf.records[3].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[3].samples[1]["GT"][0]);
    // first missing data
    EXPECT_EQ((uint) 0, vcf.records[4].samples[0]["GT"].size());
    EXPECT_EQ((uint16_t) 0, vcf.records[4].samples[1]["GT"][0]);
    // second not confident
    EXPECT_EQ((uint16_t) 1, vcf.records[5].samples[0]["GT"][0]);
    EXPECT_EQ((uint) 0, vcf.records[5].samples[1]["GT"].size());

}

TEST(VCFTest, clean) {
    VCF vcf;

    VCFRecord dummy;
    std::vector<std::string> empty = {};
    vcf.add_record(dummy, empty);
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_sample_gt("sample", "chrom1", 2, "T", "TA");
    vcf.add_sample_gt("sample", "chrom1", 5, "A", "G");
    vcf.add_sample_gt("sample", "chrom1", 79, "C", "A");
    vcf.records[2].clear();
    EXPECT_EQ((uint) 5, vcf.records.size());

    vcf.clean();
    EXPECT_EQ((uint) 3, vcf.records.size());
    EXPECT_EQ((uint) 79, vcf.records[0].pos);
    EXPECT_EQ((uint) 1, vcf.records[0].alt.size());
    EXPECT_EQ("G", vcf.records[0].alt[0]);
    EXPECT_EQ((uint) 5, vcf.records[1].pos);
    EXPECT_EQ((uint) 79, vcf.records[2].pos);
    EXPECT_EQ((uint) 1, vcf.records[2].alt.size());
    EXPECT_EQ("A", vcf.records[2].alt[0]);
}

TEST(VCFTest, add_formats) {
    VCF vcf;
    std::vector<std::string> formats = {"GT", "LIKELIHOOD", "GT_CONF", "MEAN_FWD_COVG", "MEAN_REV_COVG", "GAPS"};

    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_sample_gt("sample", "chrom1", 46, "CTT", "TA");

    vcf.add_formats(formats);
    cout << vcf << endl;

    for (const auto& record: vcf.records){
        for (const auto &f: formats){
            EXPECT_TRUE(std::find(record.format.begin(), record.format.end(), f)!=record.format.end());
        }
    }
}

TEST(VCFTest, merge_multi_allelic) {
    VCF vcf;
    // no gt
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 5, "A", "C");
    // gt
    vcf.add_record("chrom1", 46, "CTT", "A");
    vcf.add_record("chrom1", 46, "CTT", "TA");
    vcf.add_sample_gt("sample", "chrom1", 46, "CTT", "TA");
    vcf.add_sample_gt("sample", "chrom1", 46, "CTT", "A");
    // likelihoods too
    vcf.add_record("chrom1", 76, "CTT", "A");
    vcf.add_record("chrom1", 76, "CTT", "TA");
    vcf.add_sample_gt("sample", "chrom1", 76, "CTT", "TA");
    vcf.add_sample_gt("sample", "chrom1", 76, "CTT", "A");
    unordered_map<string, vector<float>> dummy;
    vcf.records[4].regt_samples.push_back(dummy);
    vcf.records[5].regt_samples.push_back(dummy);
    vcf.records[4].regt_samples[0]["LIKELIHOOD"] = {-50, -3};
    vcf.records[5].regt_samples[0]["LIKELIHOOD"] = {-50, -16};
    vcf.records[4].regt_samples[0]["GT_CONF"] = {47};
    vcf.records[5].regt_samples[0]["GT_CONF"] = {56};
    vcf.records[4].samples[0]["MEAN_FWD_COVG"] = {2, 30};
    vcf.records[5].samples[0]["MEAN_FWD_COVG"] = {2, 30};
    vcf.records[4].samples[0]["MEAN_REV_COVG"] = {2, 30};
    vcf.records[5].samples[0]["MEAN_REV_COVG"] = {2, 30};
    vcf.records[4].regt_samples[0]["GAPS"] = {4, 0};
    vcf.records[5].regt_samples[0]["GAPS"] = {4, 1};
    // incompatible
    vcf.add_record("chrom1", 85, "A", "G");
    vcf.add_record("chrom1", 85, "T", "C");

    vcf.merge_multi_allelic();
    std::vector<std::string> formats = {"GT", "LIKELIHOOD", "GT_CONF", "MEAN_FWD_COVG", "MEAN_REV_COVG", "GAPS"};
    vcf.add_formats(formats);
    cout << vcf << endl;

    EXPECT_EQ((uint) 5, vcf.records.size());
    EXPECT_EQ((uint) 5, vcf.records[0].pos);
    EXPECT_EQ((uint) 2, vcf.records[0].alt.size());
    EXPECT_EQ((uint) 1, vcf.records[0].samples.size());
    EXPECT_EQ((uint) 0, vcf.records[0].samples[0].size());

    EXPECT_EQ((uint) 46, vcf.records[1].pos);
    EXPECT_EQ((uint) 2, vcf.records[1].alt.size());
    EXPECT_EQ((uint) 1, vcf.records[1].samples.size());
    bool found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint) 0, vcf.records[1].samples[0]["GT"].size());

    EXPECT_EQ((uint) 76, vcf.records[2].pos);
    EXPECT_EQ((uint) 2, vcf.records[2].alt.size());
    EXPECT_EQ((uint) 1, vcf.records[2].samples.size());
    found_gt = vcf.records[2].samples[0].find("GT") != vcf.records[2].samples[0].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint) 1, vcf.records[2].samples[0]["GT"][0]);
    EXPECT_EQ((uint) 3, vcf.records[2].regt_samples[0].size());
    bool found = vcf.records[2].regt_samples[0].find("LIKELIHOOD") != vcf.records[2].regt_samples[0].end();
    EXPECT_TRUE(found);
    EXPECT_EQ((uint) 3, vcf.records[2].regt_samples[0]["LIKELIHOOD"].size());
    EXPECT_EQ(-50.0, vcf.records[2].regt_samples[0]["LIKELIHOOD"][0]);
    EXPECT_EQ(-3.0, vcf.records[2].regt_samples[0]["LIKELIHOOD"][1]);
    EXPECT_EQ(-16.0, vcf.records[2].regt_samples[0]["LIKELIHOOD"][2]);
    EXPECT_EQ((uint) 3, vcf.records[2].regt_samples[0]["GAPS"].size());
    EXPECT_EQ(4, vcf.records[2].regt_samples[0]["GAPS"][0]);
    EXPECT_EQ(0, vcf.records[2].regt_samples[0]["GAPS"][1]);
    EXPECT_EQ(1, vcf.records[2].regt_samples[0]["GAPS"][2]);
    found = vcf.records[2].regt_samples[0].find("GT_CONF") != vcf.records[2].regt_samples[0].end();
    EXPECT_TRUE(found);
    EXPECT_EQ((uint) 1, vcf.records[2].regt_samples[0]["GT_CONF"].size());
    EXPECT_EQ(13.0, vcf.records[2].regt_samples[0]["GT_CONF"][0]);

    EXPECT_EQ((uint) 85, vcf.records[3].pos);
    EXPECT_EQ((uint) 1, vcf.records[3].alt.size());
    EXPECT_EQ((uint) 85, vcf.records[4].pos);
    EXPECT_EQ((uint) 1, vcf.records[4].alt.size());
}

TEST(VCFTest, correct_dot_alleles) {
    VCF vcf;
    // at start
    vcf.add_sample_gt("sample", "chrom1", 0, ".", "TA");
    vcf.add_sample_gt("sample", "chrom2", 0, "T", ".");
    // in middle
    vcf.add_sample_gt("sample", "chrom1", 35, ".", "A");
    vcf.add_sample_gt("sample", "chrom2", 35, "TA", ".");
    // multiple alts
    vcf.add_sample_gt("sample", "chrom1", 44, "TA", "T");
    vcf.add_sample_gt("sample", "chrom1", 44, "TA", ".");
    vcf.add_sample_gt("sample", "chrom2", 44, ".", "T");
    vcf.add_sample_gt("sample", "chrom2", 44, ".", "TA");

    string vcf_ref = "TATATGTGTC"
            "GCGACACTGC"
            "ATGCATGCAT"
            "AGTCCTAAAG"
            "TCCTTAAACG"
            "TTTATAGTCG";

    vcf.correct_dot_alleles(vcf_ref, "chrom1");
    vcf.correct_dot_alleles(vcf_ref, "chrom2");

    EXPECT_EQ(vcf.records[0].ref, "T");
    EXPECT_EQ(vcf.records[1].ref, "C");
    EXPECT_EQ(vcf.records[2].ref, "TTA");
    EXPECT_EQ(vcf.records[3].ref, "TA");
    EXPECT_EQ(vcf.records[4].ref, "TA");
    EXPECT_EQ(vcf.records[5].ref, "CTA");
    EXPECT_EQ(vcf.records[6].ref, "T");
    EXPECT_EQ(vcf.records[7].ref, "T");

    EXPECT_EQ(vcf.records[0].alt.size(), 1);
    EXPECT_EQ(vcf.records[0].alt[0], "TAT");
    EXPECT_EQ(vcf.records[1].alt.size(), 1);
    EXPECT_EQ(vcf.records[1].alt[0], "CA");
    EXPECT_EQ(vcf.records[2].alt.size(), 1);
    EXPECT_EQ(vcf.records[2].alt[0], "T");
    EXPECT_EQ(vcf.records[3].alt.size(), 1);
    EXPECT_EQ(vcf.records[3].alt[0], "T");
    EXPECT_EQ(vcf.records[4].alt.size(), 1);
    EXPECT_EQ(vcf.records[4].alt[0], "A");
    EXPECT_EQ(vcf.records[5].alt.size(), 1);
    EXPECT_EQ(vcf.records[5].alt[0], "C");
    EXPECT_EQ(vcf.records[6].alt.size(), 1);
    EXPECT_EQ(vcf.records[6].alt[0], "TT");
    EXPECT_EQ(vcf.records[7].alt.size(), 1);
    EXPECT_EQ(vcf.records[7].alt[0], "TTA");
}

TEST(VCFTest, make_gt_compatible) {
    VCF vcf;
    // no gt
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 5, "A", "C");
    // gt incompatible no likelihoods
    vcf.add_record("chrom1", 46, "CTT", "A");
    vcf.add_record("chrom1", 46, "CTT", "TA");
    vcf.add_sample_gt("sample", "chrom1", 46, "CTT", "TA");
    vcf.add_sample_gt("sample", "chrom1", 46, "CTT", "A");
    // gt incompatible, likelihoods too both alts
    vcf.add_record("chrom1", 76, "CTT", "A");
    vcf.add_record("chrom1", 76, "CTT", "TA");
    vcf.add_sample_gt("sample", "chrom1", 76, "CTT", "TA");
    vcf.add_sample_gt("sample", "chrom1", 76, "CTT", "A");
    unordered_map<string, vector<float>> dummy;
    vcf.records[4].regt_samples.push_back(dummy);
    vcf.records[5].regt_samples.push_back(dummy);
    vcf.records[4].regt_samples[0]["LIKELIHOOD"] = {-50, -3};
    vcf.records[5].regt_samples[0]["LIKELIHOOD"] = {-50, -16};
    vcf.records[4].regt_samples[0]["GT_CONF"] = {47};
    vcf.records[5].regt_samples[0]["GT_CONF"] = {56};
    // gt incompatible one ref, ref correct
    vcf.add_record("chrom1", 85, "A", "G");
    vcf.add_record("chrom1", 85, "A", "C");
    vcf.add_sample_gt("sample", "chrom1", 85, "A", "A");
    vcf.records[6].samples[0]["GT"] = {1};
    vcf.records[6].regt_samples.push_back(dummy);
    vcf.records[7].regt_samples.push_back(dummy);
    vcf.records[6].regt_samples[0]["LIKELIHOOD"] = {-5, -30};
    vcf.records[7].regt_samples[0]["LIKELIHOOD"] = {-5, -16};
    vcf.records[6].regt_samples[0]["GT_CONF"] = {47};
    vcf.records[7].regt_samples[0]["GT_CONF"] = {56};
    // gt incompatible one ref, ref wrong
    vcf.add_record("chrom1", 95, "A", "G");
    vcf.add_record("chrom1", 95, "A", "C");
    vcf.add_sample_gt("sample", "chrom1", 95, "A", "A");
    vcf.records[8].samples[0]["GT"] = {1};
    vcf.records[8].regt_samples.push_back(dummy);
    vcf.records[9].regt_samples.push_back(dummy);
    vcf.records[8].regt_samples[0]["LIKELIHOOD"] = {-50, -3};
    vcf.records[9].regt_samples[0]["LIKELIHOOD"] = {-50, -60};
    vcf.records[8].regt_samples[0]["GT_CONF"] = {47};
    vcf.records[9].regt_samples[0]["GT_CONF"] = {10};

    vcf.make_gt_compatible();
    cout << vcf << endl;

    bool found_gt = vcf.records[0].samples[0].find("GT") != vcf.records[0].samples[0].end();
    EXPECT_FALSE(found_gt);
    found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_FALSE(found_gt);
    EXPECT_EQ((uint) 0, vcf.records[2].samples[0]["GT"].size());
    EXPECT_EQ((uint) 0, vcf.records[3].samples[0]["GT"].size());
    EXPECT_EQ((uint16_t) 1, vcf.records[4].samples[0]["GT"][0]);
    EXPECT_EQ((uint) 0, vcf.records[5].samples[0]["GT"].size());
    EXPECT_EQ((uint16_t) 0, vcf.records[6].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[7].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 1, vcf.records[8].samples[0]["GT"][0]);
    EXPECT_EQ((uint16_t) 0, vcf.records[9].samples[0]["GT"].size());
}

TEST(VCFTest, equals) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_record(vr, empty);
    EXPECT_EQ(vcf, vcf);

    // different order
    VCF vcf1;
    vcf1.add_record("chrom1", 5, "A", "G");
    vcf1.add_record(vr, empty);
    vcf1.add_record("chrom1", 46, "T", "TA");
    EXPECT_EQ(vcf1, vcf1);
    EXPECT_EQ(vcf, vcf1);
    EXPECT_EQ(vcf1, vcf);

    // same length, one different
    VCF vcf2;
    vcf2.add_record("chrom1", 10, "A", "G");
    vcf2.add_record(vr, empty);
    vcf2.add_record("chrom1", 46, "T", "TA");
    EXPECT_EQ(vcf2, vcf2);
    EXPECT_EQ((vcf == vcf2), false);
    EXPECT_EQ((vcf2 == vcf), false);

    // different length
    VCF vcf3;
    vcf3.add_record("chrom1", 5, "A", "G");
    vcf3.add_record(vr, empty);
    vcf3.add_record("chrom1", 46, "T", "TA");
    vcf3.add_record("chrom1", 30, "G", "CC");
    EXPECT_EQ(vcf3, vcf3);
    EXPECT_EQ((vcf == vcf3), false);
    EXPECT_EQ((vcf3 == vcf), false);
}

TEST(VCFTest, save) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_record(vr, empty);
    uint j = 3;
    EXPECT_EQ(j, vcf.records.size());

    vcf.save("vcf_test.vcf");
}

TEST(VCFTest, load) {
    VCF vcf, vcf1;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    std::vector<std::string> empty = {};
    vcf.add_record(vr, empty);
    uint j = 3;
    EXPECT_EQ(j, vcf.records.size());

    cout << vcf << endl;
    vcf1.load("vcf_test.vcf");
    cout << vcf1 << endl;

    EXPECT_EQ(vcf == vcf1, true);
}

TEST(VCFTest, filter) {
    VCF vcf, vcf1, vcf2, vcf3, vcf4;
    vcf.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
    vcf.add_record("chrom1", 46, "T", "TA", "SVTYPE=INDEL;GRAPHTYPE=NESTED");
    vcf.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE");
    vcf.add_record("chrom1", 79, "CTT", "ATA", "SVTYPE=PH_SNPs;GRAPHTYPE=NESTED");
    vcf.save("vcf_filter_test.vcf", true, false, false, false, false, false, false);

    vcf1.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
    vcf1.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE");
    vcf2.load("vcf_filter_test.vcf");
    EXPECT_EQ(vcf2 == vcf1, true);

    vcf.save("vcf_filter_test.vcf", false, false, false, false, false, true, false);
    vcf3.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
    vcf3.add_record("chrom1", 79, "CTT", "ATA", "SVTYPE=SNP;GRAPHTYPE=NESTED");
    vcf4.load("vcf_filter_test.vcf");
    EXPECT_EQ(vcf3 == vcf4, true);
}