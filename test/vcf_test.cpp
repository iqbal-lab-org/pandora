#include "gtest/gtest.h"
#include "vcf.h"
#include "vcfrecord.h"
#include "interval.h"
#include "localnode.h"
#include <stdint.h>
#include <iostream>

using namespace std;

TEST(VCFTest, add_record) {

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
    EXPECT_EQ((uint8_t) 1, vcf.records[1].samples[0]["GT"]);
    EXPECT_EQ(j, vcf.records[0].samples.size());
    EXPECT_TRUE(vcf.records[0].samples[0].find("GT") == vcf.records[0].samples[0].end());
    EXPECT_EQ(j, vcf.records[2].samples.size());
    EXPECT_TRUE(vcf.records[2].samples[0].find("GT") == vcf.records[2].samples[0].end());
    EXPECT_EQ(j, vcf.records[3].samples.size());
    EXPECT_TRUE(vcf.records[3].samples[0].find("GT") == vcf.records[3].samples[0].end());

    vcf.add_sample_gt("sample", "chrom1", 79, "C", "C");
    EXPECT_EQ(j, vcf.samples.size());
    EXPECT_EQ(j, vcf.records[1].samples.size());
    EXPECT_EQ((uint8_t) 1, vcf.records[1].samples[0]["GT"]);
    EXPECT_EQ(j, vcf.records[0].samples.size());
    EXPECT_TRUE(vcf.records[0].samples[0].find("GT") == vcf.records[0].samples[0].end());
    EXPECT_EQ(j, vcf.records[2].samples.size());
    EXPECT_EQ((uint8_t) 0, vcf.records[2].samples[0]["GT"]);
    EXPECT_EQ(j, vcf.records[3].samples.size());
    EXPECT_EQ((uint8_t) 0, vcf.records[3].samples[0]["GT"]);
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
    EXPECT_EQ((uint8_t) 0, vcf.records[1].samples[0]["GT"]);
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
    EXPECT_EQ((uint8_t) 0, vcf.records[0].samples[1]["GT"]);
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
    EXPECT_EQ((uint8_t) 1, vcf.records[1].samples[0]["GT"]);
    EXPECT_EQ((uint8_t) 1, vcf.records[2].samples[0]["GT"]);
    EXPECT_TRUE(vcf.records[3].samples[0].find("GT") == vcf.records[3].samples[0].end());
    EXPECT_TRUE(vcf.records[0].samples[1].find("GT") == vcf.records[0].samples[1].end());
    EXPECT_TRUE(vcf.records[1].samples[1].find("GT") == vcf.records[1].samples[1].end());
    EXPECT_EQ((uint8_t) 0, vcf.records[2].samples[1]["GT"]);
    EXPECT_EQ((uint8_t) 0, vcf.records[3].samples[1]["GT"]);

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

    vcf1.load("vcf_test.vcf");

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

TEST(VCFTest, write_aligned_fasta) {
    VCF vcf;
    vcf.add_record("chrom1", 1, "A", "G");
    vcf.add_record("chrom1", 3, "T", "TA");
    VCFRecord vr = VCFRecord("chrom1", 5, "C", "G");
    vcf.add_record(vr);
    uint j = 3;
    EXPECT_EQ(j, vcf.records.size());

    vector<LocalNodePtr> lmp;
    vcf.write_aligned_fasta("vcf1.multisample.fa", lmp);

    // add just the ref
    LocalNodePtr ln0(make_shared<LocalNode>("A", Interval(0, 1), 1));
    lmp.push_back(ln0);
    LocalNodePtr ln1(make_shared<LocalNode>("A", Interval(5, 6), 2));
    lmp.push_back(ln1);
    LocalNodePtr ln4(make_shared<LocalNode>("A", Interval(7, 8), 3));
    lmp.push_back(ln4);
    LocalNodePtr ln2(make_shared<LocalNode>("T", Interval(46, 47), 4));
    lmp.push_back(ln2);
    LocalNodePtr ln5(make_shared<LocalNode>("A", Interval(50, 51), 5));
    lmp.push_back(ln5);
    LocalNodePtr ln3(make_shared<LocalNode>("C", Interval(79, 80), 6));
    lmp.push_back(ln3);
    vcf.write_aligned_fasta("vcf2.multisample.fa", lmp);

    // now add a sample
    vcf.add_sample_gt("sample1", "chrom1", 46, "T", "TA");
    vcf.write_aligned_fasta("vcf3.multisample.fa", lmp);

}
