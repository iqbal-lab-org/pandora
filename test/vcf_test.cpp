#include "gtest/gtest.h"
#include "vcf.h"
#include "vcfrecord.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class VCFTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(VCFTest,addRecord){

    VCF vcf;
    uint j = 0;
    EXPECT_EQ(j, vcf.records.size());

    vcf.add_record("chrom1", 5, "A", "G");
    j = 1;
    EXPECT_EQ(j, vcf.records.size());

    // add the same one again
    vcf.add_record("chrom1", 5, "A", "G");
    EXPECT_EQ(j, vcf.records.size());

    // add a different one
    vcf.add_record("chrom1", 46, "T", "TA");
    j = 2;
    EXPECT_EQ(j, vcf.records.size());

    // add the first one again
    vcf.add_record("chrom1", 5, "A", "G");
    EXPECT_EQ(j, vcf.records.size());

    // use the other addition method
    VCFRecord vr("chrom1", 5, "A", "G");
    vcf.add_record(vr);
    EXPECT_EQ(j, vcf.records.size());

    // use the other addition method for a new one
    vr = VCFRecord("chrom1", 79, "C", "G");
    vcf.add_record(vr);
    j = 3;
    EXPECT_EQ(j, vcf.records.size());
}

TEST_F(VCFTest,clear){
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

TEST_F(VCFTest,equals){
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
    EXPECT_EQ((vcf==vcf2), false);
    EXPECT_EQ((vcf2==vcf), false);

    // different length
    VCF vcf3;
    vcf3.add_record("chrom1", 5, "A", "G");
    vcf3.add_record(vr);
    vcf3.add_record("chrom1", 46, "T", "TA");
    vcf3.add_record("chrom1", 30, "G", "CC");
    EXPECT_EQ(vcf3, vcf3);
    EXPECT_EQ((vcf==vcf3), false);
    EXPECT_EQ((vcf3==vcf), false);
}

TEST_F(VCFTest,save){
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    vcf.add_record(vr);
    uint j = 3;
    EXPECT_EQ(j, vcf.records.size());

    vcf.save("../test/test_cases/vcf_test.vcf");
}

TEST_F(VCFTest,load){
    VCF vcf, vcf1;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 79, "C", "G");
    vcf.add_record(vr);
    uint j = 3;
    EXPECT_EQ(j, vcf.records.size());

    vcf1.load("../test/test_cases/vcf_test.vcf");

    /*for(uint i=0; i!=vcf1.records.size(); ++i)
    {
        cout << vcf1.records[i];
    }*/
    EXPECT_EQ(vcf == vcf1, true);
}
