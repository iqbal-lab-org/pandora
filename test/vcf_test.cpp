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
    EXPECT_EQ((uint)0, vcf.records.size());
    vcf.add_record("chrom1", 5, "A", "G");
    EXPECT_EQ((uint)1, vcf.records.size());
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
    vcf.add_record(vr);
    EXPECT_EQ((uint) 1, vcf.records.size());
}

TEST(VCFTest, add_record_by_record_and_values) {
    VCF vcf;
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    vcf.add_record(vr);
    vcf.add_record("chrom1", 79, "C", "G");
    EXPECT_EQ((uint) 1, vcf.records.size());
}

TEST(VCFTest, add_record_by_values_and_record) {
    VCF vcf;
    vcf.add_record("chrom1", 79, "C", "G");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    vcf.add_record(vr);
    EXPECT_EQ((uint) 1, vcf.records.size());
}

TEST(VCFTest, add_record_by_record_returned_by_reference) {
    VCF vcf;
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    VCFRecord& ref_vr = vcf.add_record(vr);
    EXPECT_EQ(ref_vr.chrom, "chrom1");
    EXPECT_EQ(ref_vr.pos, (uint)79);
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
    EXPECT_EQ((uint8_t) 1, vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[0].samples.size());
    EXPECT_TRUE(vcf.records[0].samples[0].find("GT") == vcf.records[0].samples[0].end());
    EXPECT_EQ(j, vcf.records[2].samples.size());
    EXPECT_TRUE(vcf.records[2].samples[0].find("GT") == vcf.records[2].samples[0].end());
    EXPECT_EQ(j, vcf.records[3].samples.size());
    EXPECT_TRUE(vcf.records[3].samples[0].find("GT") == vcf.records[3].samples[0].end());

    vcf.add_sample_gt("sample", "chrom1", 79, "C", "C");
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.records[1].samples.size());
    EXPECT_EQ((uint8_t) 1, vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[0].samples.size());
    EXPECT_TRUE(vcf.records[0].samples[0].find("GT") == vcf.records[0].samples[0].end());
    EXPECT_EQ(j, vcf.records[2].samples.size());
    EXPECT_EQ((uint8_t) 0, vcf.records[2].samples[0]["GT"][0]);
    EXPECT_EQ(j, vcf.records[3].samples.size());
    EXPECT_EQ((uint8_t) 0, vcf.records[3].samples[0]["GT"][0]);
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
    EXPECT_EQ((uint8_t) 0, vcf.records[1].samples[0]["GT"][0]);
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
    EXPECT_EQ((uint8_t) 0, vcf.records[0].samples[1]["GT"][0]);
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
    EXPECT_EQ((uint8_t) 1, vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ((uint8_t) 1, vcf.records[2].samples[0]["GT"][0]);
    EXPECT_TRUE(vcf.records[3].samples[0].find("GT") == vcf.records[3].samples[0].end());
    EXPECT_TRUE(vcf.records[0].samples[1].find("GT") == vcf.records[0].samples[1].end());
    EXPECT_TRUE(vcf.records[1].samples[1].find("GT") == vcf.records[1].samples[1].end());
    EXPECT_EQ((uint8_t) 0, vcf.records[2].samples[1]["GT"][0]);
    EXPECT_EQ((uint8_t) 0, vcf.records[3].samples[1]["GT"][0]);

}


TEST(VCFTest, clear) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    vcf.add_record(vr);
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
    EXPECT_EQ((uint)8, vcf.records.size());
    for (uint i=0; i<4; ++i) {
        EXPECT_EQ(vcf.records[i].chrom, "chrom1");
    }
    for (uint i=4; i<8; ++i) {
        EXPECT_EQ(vcf.records[i].chrom, "chrom2");
    }
    EXPECT_EQ((uint)5,vcf.records[4].pos);
    EXPECT_EQ("TA", vcf.records[5].alt[0]);
    EXPECT_EQ((uint)79,vcf.records[6].pos);
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
    EXPECT_EQ((uint)6, vcf.records.size());
    for (uint i=0; i<4; ++i) {
        EXPECT_EQ(vcf.records[i].chrom, "chrom1");
    }
    for (uint i=4; i<6; ++i) {
        EXPECT_EQ(vcf.records[i].chrom, "chrom2");
    }
    EXPECT_EQ((uint)5,vcf.records[4].pos);
    EXPECT_EQ((uint)79,vcf.records[5].pos);
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
    EXPECT_EQ((uint)1, vcf.samples.size());
    EXPECT_EQ("sample", vcf.samples[0]);
    EXPECT_EQ((uint)1, vcf.records[0].samples.size());
    EXPECT_EQ((uint)1, vcf.records[5].samples.size());
    bool found_gt = vcf.records[2].samples[0].find("GT") != vcf.records[2].samples[0].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint)1, vcf.records[2].samples[0]["GT"][0]);
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
    EXPECT_EQ((uint)1, vcf.samples.size());
    EXPECT_EQ("sample", vcf.samples[0]);
    EXPECT_EQ((uint)1, vcf.records[0].samples.size());
    EXPECT_EQ((uint)1, vcf.records[5].samples.size());
    bool found_gt = vcf.records[4].samples[0].find("GT") != vcf.records[4].samples[0].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint)1, vcf.records[4].samples[0]["GT"][0]);
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
    EXPECT_EQ((uint)1, vcf.samples.size());
    EXPECT_EQ("sample", vcf.samples[0]);
    EXPECT_EQ((uint)1, vcf.records[0].samples.size());
    EXPECT_EQ((uint)1, vcf.records[5].samples.size());
    bool found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint)1, vcf.records[1].samples[0]["GT"][0]);
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

    EXPECT_EQ((uint)2, vcf.samples.size());
    vector<string> v = {"sample", "sample1"};
    EXPECT_ITERABLE_EQ(vector<string>, v, vcf.samples);
    EXPECT_EQ((uint)2, vcf.records[0].samples.size());
    EXPECT_EQ((uint)2, vcf.records[5].samples.size());


    bool found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint)1, vcf.records[1].samples[0]["GT"][0]);
    found_gt = vcf.records[1].samples[1].find("GT") != vcf.records[1].samples[1].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint)0, vcf.records[1].samples[1]["GT"][0]);
    found_gt = vcf.records[3].samples[1].find("GT") != vcf.records[3].samples[1].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint)1, vcf.records[3].samples[1]["GT"][0]);

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
    for (uint i = 0; i < 4; ++i){
        EXPECT_EQ("chrom1", vcf.records[i].chrom);
    }
    for (uint i = 4; i < 6; ++i){
        EXPECT_EQ("chrom2", vcf.records[i].chrom);
    }
    EXPECT_EQ((uint)5, vcf.records[0].pos);
    EXPECT_EQ((uint)5, vcf.records[4].pos);
    EXPECT_EQ((uint)46, vcf.records[1].pos);
    EXPECT_EQ((uint)79, vcf.records[2].pos);
    EXPECT_EQ((uint)79, vcf.records[3].pos);
    EXPECT_EQ((uint)79, vcf.records[5].pos);
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

    EXPECT_TRUE(vcf.pos_in_range(4,6,"chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(5,6,"chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(4,5,"chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(4,6,"chrom2"));

    EXPECT_TRUE(vcf.pos_in_range(45,47,"chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(46,47,"chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(45,46,"chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(45,47,"chrom2"));

    EXPECT_TRUE(vcf.pos_in_range(78,80,"chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(79,80,"chrom1"));
    EXPECT_FALSE(vcf.pos_in_range(78,79,"chrom1"));
    EXPECT_TRUE(vcf.pos_in_range(78,80,"chrom2"));

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

    vcf.add_sample_gt("sample","chrom1", 2, "T", "TA");
    vcf.add_sample_gt("sample","chrom1", 5, "A", "G");
    vcf.add_sample_gt("sample","chrom1", 79, "C", "A");
    vcf.add_sample_gt("sample","chrom2", 20, "A", "G");
    vcf.add_sample_gt("sample","chrom2", 79, "C", "C");
    vcf.add_sample_gt("sample","chrom2", 80, "A", "C");

    vcf.add_sample_gt("asample","chrom1", 2, "T", "TA");
    vcf.add_sample_gt("asample","chrom1", 5, "A", "A");
    vcf.add_sample_gt("asample","chrom1", 79, "C", "A");
    vcf.add_sample_gt("asample","chrom2", 20, "A", "G");
    vcf.add_sample_gt("asample","chrom2", 79, "C", "C");
    vcf.add_sample_gt("asample","chrom2", 80, "A", "A");

    vcf.sort_records();

    // record 0, not a snp site
    vcf.records[0].samples[0]["MEAN_FWD_COVG"] = {0,10};
    vcf.records[0].samples[0]["MEAN_REV_COVG"] = {1,20};
    vcf.records[0].samples[1]["MEAN_FWD_COVG"] = {1,15};
    vcf.records[0].samples[1]["MEAN_REV_COVG"] = {2,24};

    // record 1, different genotypes but both correct
    vcf.records[1].samples[0]["MEAN_FWD_COVG"] = {0,10};
    vcf.records[1].samples[0]["MEAN_REV_COVG"] = {1,20};
    vcf.records[1].samples[1]["MEAN_FWD_COVG"] = {10,1};
    vcf.records[1].samples[1]["MEAN_REV_COVG"] = {21,2};

    // record 2, same genotypes first correct
    vcf.records[2].samples[0]["MEAN_FWD_COVG"] = {0,10};
    vcf.records[2].samples[0]["MEAN_REV_COVG"] = {1,20};
    vcf.records[2].samples[1]["MEAN_FWD_COVG"] = {10,1};
    vcf.records[2].samples[1]["MEAN_REV_COVG"] = {21,2};

    // record 3, same genotypes both wrong
    vcf.records[3].samples[0]["MEAN_FWD_COVG"] = {20,1};
    vcf.records[3].samples[0]["MEAN_REV_COVG"] = {21,2};
    vcf.records[3].samples[1]["MEAN_FWD_COVG"] = {10,1};
    vcf.records[3].samples[1]["MEAN_REV_COVG"] = {21,2};

    // record 4, missing count data for first sample
    vcf.records[4].samples[0]["MEAN_FWD_COVG"] = {0,10};
    vcf.records[4].samples[0]["MEAN_REV_COVG"] = {20};
    vcf.records[4].samples[1]["MEAN_FWD_COVG"] = {10,1};
    vcf.records[4].samples[1]["MEAN_REV_COVG"] = {21,2};

    // record 5, not confident for second sample
    vcf.records[5].samples[0]["MEAN_FWD_COVG"] = {0,10};
    vcf.records[5].samples[0]["MEAN_REV_COVG"] = {1,20};
    vcf.records[5].samples[1]["MEAN_FWD_COVG"] = {2,1};
    vcf.records[5].samples[1]["MEAN_REV_COVG"] = {4,2};

    cout << vcf << endl;

    vcf.genotype(30,0.01,30);

    cout << vcf << endl;

    // not genotyped first record
    EXPECT_EQ((uint8_t)1,vcf.records[0].samples[0]["GT"][0]);
    EXPECT_EQ((uint8_t)1,vcf.records[0].samples[1]["GT"][0]);
    bool found_confidence = !vcf.records[0].regt_samples.empty();
    EXPECT_FALSE(found_confidence);
    // both correct
    EXPECT_EQ((uint)2, vcf.records[1].samples.size());
    bool found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_TRUE(found_gt);
    found_gt = vcf.records[1].samples[1].find("GT") != vcf.records[1].samples[1].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint)1, vcf.records[1].samples[0]["GT"].size());
    EXPECT_EQ((uint)1, vcf.records[1].samples[1]["GT"].size());
    EXPECT_EQ((uint8_t)1,vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ((uint8_t)0,vcf.records[1].samples[1]["GT"][0]);
    // first correct
    EXPECT_EQ((uint8_t)1,vcf.records[2].samples[0]["GT"][0]);
    EXPECT_EQ((uint8_t)0,vcf.records[2].samples[1]["GT"][0]);
    // both wrong
    EXPECT_EQ((uint8_t)0,vcf.records[3].samples[0]["GT"][0]);
    EXPECT_EQ((uint8_t)0,vcf.records[3].samples[1]["GT"][0]);
    // first missing data
    EXPECT_EQ((uint)0, vcf.records[4].samples[0]["GT"].size());
    EXPECT_EQ((uint8_t)0,vcf.records[4].samples[1]["GT"][0]);
    // second not confident
    EXPECT_EQ((uint8_t)1,vcf.records[5].samples[0]["GT"][0]);
    EXPECT_EQ((uint)0, vcf.records[5].samples[1]["GT"].size());
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

    vcf.add_sample_gt("sample","chrom1", 2, "T", "TA");
    vcf.add_sample_gt("sample","chrom1", 5, "AC", "GC");
    vcf.add_sample_gt("sample","chrom1", 79, "CC", "AC");
    vcf.add_sample_gt("sample","chrom2", 20, "AC", "GC");
    vcf.add_sample_gt("sample","chrom2", 79, "CC", "CC");
    vcf.add_sample_gt("sample","chrom2", 80, "AC", "CC");

    vcf.add_sample_gt("asample","chrom1", 2, "T", "TA");
    vcf.add_sample_gt("asample","chrom1", 5, "AC", "AC");
    vcf.add_sample_gt("asample","chrom1", 79, "CC", "AC");
    vcf.add_sample_gt("asample","chrom2", 20, "AC", "GC");
    vcf.add_sample_gt("asample","chrom2", 79, "CC", "CC");
    vcf.add_sample_gt("asample","chrom2", 80, "AC", "AC");

    vcf.sort_records();

    // record 0, not a snp site
    vcf.records[0].samples[0]["MEAN_FWD_COVG"].push_back(0);
    vcf.records[0].samples[0]["MEAN_REV_COVG"].push_back(1);
    vcf.records[0].samples[0]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[0].samples[0]["MEAN_REV_COVG"].push_back(20);
    vcf.records[0].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[0].samples[1]["MEAN_REV_COVG"].push_back(2);
    vcf.records[0].samples[1]["MEAN_FWD_COVG"].push_back(15);
    vcf.records[0].samples[1]["MEAN_REV_COVG"].push_back(24);

    // record 1, different genotypes but both correct
    vcf.records[1].samples[0]["MEAN_FWD_COVG"].push_back(0);
    vcf.records[1].samples[0]["MEAN_REV_COVG"].push_back(1);
    vcf.records[1].samples[0]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[1].samples[0]["MEAN_REV_COVG"].push_back(20);
    vcf.records[1].samples[1]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[1].samples[1]["MEAN_REV_COVG"].push_back(21);
    vcf.records[1].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[1].samples[1]["MEAN_REV_COVG"].push_back(2);

    // record 2, same genotypes first correct
    vcf.records[2].samples[0]["MEAN_FWD_COVG"].push_back(0);
    vcf.records[2].samples[0]["MEAN_REV_COVG"].push_back(1);
    vcf.records[2].samples[0]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[2].samples[0]["MEAN_REV_COVG"].push_back(20);
    vcf.records[2].samples[1]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[2].samples[1]["MEAN_REV_COVG"].push_back(21);
    vcf.records[2].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[2].samples[1]["MEAN_REV_COVG"].push_back(2);

    // record 3, same genotypes both wrong
    vcf.records[3].samples[0]["MEAN_FWD_COVG"].push_back(20);
    vcf.records[3].samples[0]["MEAN_REV_COVG"].push_back(21);
    vcf.records[3].samples[0]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[3].samples[0]["MEAN_REV_COVG"].push_back(2);
    vcf.records[3].samples[1]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[3].samples[1]["MEAN_REV_COVG"].push_back(21);
    vcf.records[3].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[3].samples[1]["MEAN_REV_COVG"].push_back(2);

    // record 4, missing count data for first sample
    vcf.records[4].samples[0]["MEAN_FWD_COVG"].push_back(0);
    vcf.records[4].samples[0]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[4].samples[0]["MEAN_REV_COVG"].push_back(20);
    vcf.records[4].samples[1]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[4].samples[1]["MEAN_REV_COVG"].push_back(21);
    vcf.records[4].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[4].samples[1]["MEAN_REV_COVG"].push_back(2);

    // record 5, not confident for second sample
    vcf.records[5].samples[0]["MEAN_FWD_COVG"].push_back(0);
    vcf.records[5].samples[0]["MEAN_REV_COVG"].push_back(1);
    vcf.records[5].samples[0]["MEAN_FWD_COVG"].push_back(10);
    vcf.records[5].samples[0]["MEAN_REV_COVG"].push_back(20);
    vcf.records[5].samples[1]["MEAN_FWD_COVG"].push_back(2);
    vcf.records[5].samples[1]["MEAN_REV_COVG"].push_back(4);
    vcf.records[5].samples[1]["MEAN_FWD_COVG"].push_back(1);
    vcf.records[5].samples[1]["MEAN_REV_COVG"].push_back(2);

    cout << vcf << endl;

    bool snps_only = false;
    vcf.genotype(30,0.01,30,snps_only);

    cout << vcf << endl;

    // first record now genotyped
    EXPECT_EQ((uint8_t)1,vcf.records[0].samples[0]["GT"][0]);
    EXPECT_EQ((uint8_t)1,vcf.records[0].samples[1]["GT"][0]);
    bool found_confidence = !vcf.records[0].regt_samples.empty();
    EXPECT_TRUE(found_confidence);
    // both correct
    EXPECT_EQ((uint)2, vcf.records[1].samples.size());
    bool found_gt = vcf.records[1].samples[0].find("GT") != vcf.records[1].samples[0].end();
    EXPECT_TRUE(found_gt);
    found_gt = vcf.records[1].samples[1].find("GT") != vcf.records[1].samples[1].end();
    EXPECT_TRUE(found_gt);
    EXPECT_EQ((uint)1, vcf.records[1].samples[0]["GT"].size());
    EXPECT_EQ((uint)1, vcf.records[1].samples[1]["GT"].size());
    EXPECT_EQ((uint8_t)1,vcf.records[1].samples[0]["GT"][0]);
    EXPECT_EQ((uint8_t)0,vcf.records[1].samples[1]["GT"][0]);
    // first correct
    EXPECT_EQ((uint8_t)1,vcf.records[2].samples[0]["GT"][0]);
    EXPECT_EQ((uint8_t)0,vcf.records[2].samples[1]["GT"][0]);
    // both wrong
    EXPECT_EQ((uint8_t)0,vcf.records[3].samples[0]["GT"][0]);
    EXPECT_EQ((uint8_t)0,vcf.records[3].samples[1]["GT"][0]);
    // first missing data
    EXPECT_EQ((uint)0, vcf.records[4].samples[0]["GT"].size());
    EXPECT_EQ((uint8_t)0,vcf.records[4].samples[1]["GT"][0]);
    // second not confident
    EXPECT_EQ((uint8_t)1,vcf.records[5].samples[0]["GT"][0]);
    EXPECT_EQ((uint)0, vcf.records[5].samples[1]["GT"].size());

}


TEST(VCFTest, equals) {
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    vcf.add_record(vr);
    EXPECT_EQ(vcf, vcf);

    // different order
    VCF vcf1;
    vcf1.add_record("chrom1", 5, "A", "G");
    vcf1.add_record(vr);
    vcf1.add_record("chrom1", 46, "T", "TA");
    EXPECT_EQ(vcf1, vcf1);
    EXPECT_EQ(vcf, vcf1);
    EXPECT_EQ(vcf1, vcf);

    // same length, one different
    VCF vcf2;
    vcf2.add_record("chrom1", 10, "A", "G");
    vcf2.add_record(vr);
    vcf2.add_record("chrom1", 46, "T", "TA");
    EXPECT_EQ(vcf2, vcf2);
    EXPECT_EQ((vcf == vcf2), false);
    EXPECT_EQ((vcf2 == vcf), false);

    // different length
    VCF vcf3;
    vcf3.add_record("chrom1", 5, "A", "G");
    vcf3.add_record(vr);
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
    vcf.add_record(vr);
    uint j = 3;
    EXPECT_EQ(j, vcf.records.size());

    vcf.save("vcf_test.vcf");
}

TEST(VCFTest, load) {
    VCF vcf, vcf1;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    vcf.add_record(vr);
    uint j = 3;
    EXPECT_EQ(j, vcf.records.size());

    cout << vcf << endl;
    vcf1.load("vcf_test.vcf");
    cout << vcf1 << endl;

    /*for(uint i=0; i!=vcf1.records.size(); ++i)
    {
        cout << vcf1.records[i];
    }*/
    EXPECT_EQ(vcf == vcf1, true);
}

TEST(VCFTest, filter) {
    VCF vcf, vcf1, vcf2, vcf3, vcf4;
    vcf.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
    vcf.add_record("chrom1", 46, "T", "TA", "SVTYPE=INDEL;GRAPHTYPE=NESTED");
    vcf.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE");
    vcf.add_record("chrom1", 79, "CTT", "ATA", "SVTYPE=PH_SNPs;GRAPHTYPE=NESTED");
    vcf.samples.push_back("dummy");

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