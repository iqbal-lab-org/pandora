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

TEST_F(VCFTest,add_record){

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

TEST_F(VCFTest, add_sample_gt)
{
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");

    vcf.add_sample_gt("sample", "chrom1", 46, "T", "TA");
    uint j = 1;
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.records[1].samples.size());
    EXPECT_EQ("1", vcf.records[1].samples[0]);
    EXPECT_EQ(j, vcf.records[0].samples.size());
    EXPECT_EQ(".", vcf.records[0].samples[0]);
    EXPECT_EQ(j, vcf.records[2].samples.size());
    EXPECT_EQ(".", vcf.records[2].samples[0]);

    vcf.add_sample_gt("sample", "chrom1", 79, "C", "C");
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.records[1].samples.size());
    EXPECT_EQ("1", vcf.records[1].samples[0]);
    EXPECT_EQ(j, vcf.records[0].samples.size());
    EXPECT_EQ(".", vcf.records[0].samples[0]);
    EXPECT_EQ(j, vcf.records[2].samples.size());
    EXPECT_EQ("0", vcf.records[2].samples[0]);
    EXPECT_EQ(j, vcf.records[3].samples.size());
    EXPECT_EQ("0", vcf.records[3].samples[0]);
}

TEST_F(VCFTest, add_sample_ref_alleles)
{
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_record("chrom1", 79, "C", "A");
    vcf.add_record("chrom2", 30, "C", "A");

    vcf.add_sample_ref_alleles("sample", "chrom1", 15, 78);
    EXPECT_EQ((uint)1, vcf.samples.size());
    EXPECT_EQ((uint)5, vcf.records.size());
    EXPECT_EQ((uint)1, vcf.records[0].samples.size());
    EXPECT_EQ(".", vcf.records[0].samples[0]);
    EXPECT_EQ((uint)1, vcf.records[1].samples.size());
    EXPECT_EQ("0", vcf.records[1].samples[0]);
    EXPECT_EQ((uint)1, vcf.records[2].samples.size());
    EXPECT_EQ(".", vcf.records[2].samples[0]);
    EXPECT_EQ((uint)1, vcf.records[3].samples.size());
    EXPECT_EQ(".", vcf.records[3].samples[0]);
    EXPECT_EQ((uint)1, vcf.records[4].samples.size());
    EXPECT_EQ(".", vcf.records[4].samples[0]);

    vcf.add_sample_ref_alleles("sample2", "chrom1", 5, 46);
    EXPECT_EQ((uint)2, vcf.samples.size());
    EXPECT_EQ((uint)5, vcf.records.size());
    EXPECT_EQ((uint)2, vcf.records[0].samples.size());
    EXPECT_EQ("0", vcf.records[0].samples[1]);
    EXPECT_EQ((uint)2, vcf.records[1].samples.size());
    EXPECT_EQ(".", vcf.records[1].samples[1]);
    EXPECT_EQ((uint)2, vcf.records[2].samples.size());
    EXPECT_EQ(".", vcf.records[2].samples[1]);
    EXPECT_EQ((uint)2, vcf.records[3].samples.size());
    EXPECT_EQ(".", vcf.records[3].samples[1]);
    EXPECT_EQ((uint)2, vcf.records[4].samples.size());
    EXPECT_EQ(".", vcf.records[4].samples[1]);
}

TEST_F(VCFTest, reorder_add_record_and_sample)
{
    VCF vcf;
    vcf.add_record("chrom1", 5, "A", "G");
    vcf.add_record("chrom1", 46, "T", "TA");
    vcf.add_sample_gt("sample1", "chrom1", 46, "T", "TA");
    vcf.add_record("chrom1", 79, "C", "G");
    vcf.add_sample_gt("sample2", "chrom1", 79, "C", "C");
    vcf.add_sample_gt("sample1", "chrom1", 79, "C", "A");
    vcf.sort_records();

    EXPECT_EQ((uint)2, vcf.samples.size());
    EXPECT_EQ((uint)4, vcf.records.size());
    EXPECT_EQ((uint)2, vcf.records[0].samples.size());
    EXPECT_EQ((uint)2, vcf.records[1].samples.size());
    EXPECT_EQ((uint)2, vcf.records[2].samples.size());
    EXPECT_EQ((uint)2, vcf.records[3].samples.size());
    EXPECT_EQ(".", vcf.records[0].samples[0]);
    EXPECT_EQ("1", vcf.records[1].samples[0]);
    EXPECT_EQ("1", vcf.records[2].samples[0]);
    EXPECT_EQ(".", vcf.records[3].samples[0]);
    EXPECT_EQ(".", vcf.records[0].samples[1]);
    EXPECT_EQ(".", vcf.records[1].samples[1]);
    EXPECT_EQ("0", vcf.records[2].samples[1]);
    EXPECT_EQ("0", vcf.records[3].samples[1]);

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

TEST_F(VCFTest, filter)
{
    VCF vcf, vcf1, vcf2, vcf3, vcf4;
    vcf.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
    vcf.add_record("chrom1", 46, "T", "TA", "SVTYPE=INDEL;GRAPHTYPE=NESTED");
    vcf.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE");
    vcf.add_record("chrom1", 79, "CTT", "ATA", "SVTYPE=PH_SNPs;GRAPHTYPE=NESTED");
    vcf.samples.push_back("dummy");

    vcf.save("../test/test_cases/vcf_filter_test.vcf", true, false, false, false, false, false, false);
    vcf1.add_record("chrom1", 5, "A", "G", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
    vcf1.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=PH_SNPs;GRAPHTYPE=SIMPLE");
    vcf2.load("../test/test_cases/vcf_filter_test.vcf");
    EXPECT_EQ(vcf2 == vcf1, true);

    vcf.save("../test/test_cases/vcf_filter_test.vcf", false, false, false, false, false, true, false);
    vcf3.add_record("chrom1", 79, "CTT", "GTA", "SVTYPE=SNP;GRAPHTYPE=SIMPLE");
    vcf3.add_record("chrom1", 79, "CTT", "ATA", "SVTYPE=SNP;GRAPHTYPE=NESTED");
    vcf4.load("../test/test_cases/vcf_filter_test.vcf");
    EXPECT_EQ(vcf3 == vcf4, true);

}
