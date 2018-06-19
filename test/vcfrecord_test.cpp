#include "gtest/gtest.h"
#include "vcfrecord.h"
#include <stdint.h>
#include <iostream>

using namespace std;

TEST(VCFRecordTest, create) {

    VCFRecord vr;
    EXPECT_EQ(".", vr.chrom);
    EXPECT_EQ((uint) 0, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ(".", vr.ref);
    EXPECT_EQ(".", vr.alt);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);

    vr = VCFRecord("chrom1", 3, "A", "T");
    EXPECT_EQ("chrom1", vr.chrom);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.ref);
    EXPECT_EQ("T", vr.alt);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ("SVTYPE=SNP", vr.info);
}

TEST(VCFRecordLikelihoodTest, does_not_crash_with_no_samples) {
    VCFRecord vr("chrom1", 3, "A", "T");
    EXPECT_NO_FATAL_FAILURE(vr.likelihood(1,0.01));
}

TEST(VCFRecordLikelihoodTest, does_not_run_if_info_missing) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    m["nothing"] = 0;
    vr.samples.push_back(m);
    vr.likelihood(1,0.01);
    bool found_likelihood = vr.samples[0].find("REF_LIKELIHOOD") != vr.samples[0].end();
    EXPECT_FALSE(found_likelihood);

    vr.samples[0]["GT"] = 1;
    vr.likelihood(1,0.01);
    found_likelihood = vr.samples[0].find("REF_LIKELIHOOD") != vr.samples[0].end();
    EXPECT_FALSE(found_likelihood);

    vr.samples[0]["REF_MEAN_FWD_COVG"] = 1;
    vr.samples[0]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[0]["ALT_MEAN_FWD_COVG"] = 1;
    vr.likelihood(1,0.01);
    found_likelihood = vr.samples[0].find("REF_LIKELIHOOD") != vr.samples[0].end();
    EXPECT_FALSE(found_likelihood);

    vr.samples[0].erase("ALT_MEAN_FWD_COVG");
    vr.samples[0]["ALT_MEAN_REV_COVG"] = 1;
    vr.likelihood(1,0.01);
    found_likelihood = vr.samples[0].find("REF_LIKELIHOOD") != vr.samples[0].end();
    EXPECT_FALSE(found_likelihood);

    vr.samples[0].erase("REF_MEAN_FWD_COVG");
    vr.samples[0]["ALT_MEAN_FWD_COVG"] = 1;
    vr.likelihood(1,0.01);
    found_likelihood = vr.samples[0].find("REF_LIKELIHOOD") != vr.samples[0].end();
    EXPECT_FALSE(found_likelihood);

    vr.samples[0].erase("REF_MEAN_REV_COVG");
    vr.samples[0]["REF_MEAN_FWD_COVG"] = 1;
    vr.likelihood(1,0.01);
    found_likelihood = vr.samples[0].find("REF_LIKELIHOOD") != vr.samples[0].end();
    EXPECT_FALSE(found_likelihood);
}

TEST(VCFRecordLikelihoodTest, adds_likelihood_with_info) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.samples[0]["REF_MEAN_FWD_COVG"] = 1;
    vr.samples[0]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[0]["ALT_MEAN_FWD_COVG"] = 2;
    vr.samples[0]["ALT_MEAN_REV_COVG"] = 2;
    vr.likelihood(1, 0.01);
    bool found_likelihood = vr.samples[0].find("REF_LIKELIHOOD") != vr.samples[0].end();
    EXPECT_TRUE(found_likelihood);
    found_likelihood = vr.samples[0].find("ALT_LIKELIHOOD") != vr.samples[0].end();
    EXPECT_TRUE(found_likelihood);
}

TEST(VCFRecordLikelihoodTest, gets_correct_likelihood_simple_case) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.samples[0]["REF_MEAN_FWD_COVG"] = 1;
    vr.samples[0]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[0]["ALT_MEAN_FWD_COVG"] = 2;
    vr.samples[0]["ALT_MEAN_REV_COVG"] = 2;
    vr.likelihood(1, 0.01);
    float exp_likelihood = -1-log(2) + 4*log(0.01);
    EXPECT_FLOAT_EQ(exp_likelihood,vr.samples[0]["REF_LIKELIHOOD"]);
    exp_likelihood = -1-log(4)-log(3)-log(2)+2*log(0.01);
    EXPECT_FLOAT_EQ(exp_likelihood,vr.samples[0]["ALT_LIKELIHOOD"]);
}

TEST(VCFRecordLikelihoodTest, handles_ref_covg_0) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.samples[0]["REF_MEAN_FWD_COVG"] = 0;
    vr.samples[0]["REF_MEAN_REV_COVG"] = 0;
    vr.samples[0]["ALT_MEAN_FWD_COVG"] = 2;
    vr.samples[0]["ALT_MEAN_REV_COVG"] = 2;
    vr.likelihood(1, 0.01);
    float exp_likelihood = numeric_limits<float>::lowest();
    EXPECT_FLOAT_EQ(exp_likelihood,vr.samples[0]["REF_LIKELIHOOD"]);
    exp_likelihood = -1-log(4)-log(3)-log(2);
    EXPECT_FLOAT_EQ(exp_likelihood,vr.samples[0]["ALT_LIKELIHOOD"]);
}

TEST(VCFRecordLikelihoodTest, handles_alt_covg_0) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.samples[0]["REF_MEAN_FWD_COVG"] = 1;
    vr.samples[0]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[0]["ALT_MEAN_FWD_COVG"] = 0;
    vr.samples[0]["ALT_MEAN_REV_COVG"] = 0;
    vr.likelihood(1, 0.01);
    float exp_likelihood = numeric_limits<float>::lowest();
    EXPECT_FLOAT_EQ(exp_likelihood,vr.samples[0]["ALT_LIKELIHOOD"]);
    exp_likelihood = -1-log(2);
    EXPECT_FLOAT_EQ(exp_likelihood,vr.samples[0]["REF_LIKELIHOOD"]);
}

TEST(VCFRecordConfidenceTest, does_not_run_if_info_missing) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.confidence();
    bool found_confidence = vr.samples[0].find("CONFIDENCE") != vr.samples[0].end();
    EXPECT_FALSE(found_confidence);

    vr.samples[0]["REF_LIKELIHOOD"] = -1.0;
    vr.confidence();
    found_confidence = vr.samples[0].find("CONFIDENCE") != vr.samples[0].end();
    EXPECT_FALSE(found_confidence);

    vr.samples[0]["ALT_LIKELIHOOD"] = -1.0;
    vr.samples[0].erase("REF_LIKELIHOOD");
    vr.confidence();
    found_confidence = vr.samples[0].find("CONFIDENCE") != vr.samples[0].end();
    EXPECT_FALSE(found_confidence);
}

TEST(VCFRecordConfidenceTest, adds_confidence_with_info) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.samples[0]["REF_LIKELIHOOD"] = -1.0;
    vr.samples[0]["ALT_LIKELIHOOD"] = 0.0;
    vr.confidence();
    bool found_confidence = vr.samples[0].find("CONFIDENCE") != vr.samples[0].end();
    EXPECT_TRUE(found_confidence);
}

TEST(VCFRecordConfidenceTest, gets_correct_confidence_simple_case) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.samples[0]["REF_LIKELIHOOD"] = -1.0;
    vr.samples[0]["ALT_LIKELIHOOD"] = 0.0;
    vr.confidence();
    float exp_confidence = 1;
    EXPECT_FLOAT_EQ(exp_confidence,vr.samples[0]["CONFIDENCE"]);
}

TEST(VCFRecordConfidenceTest, handles_ref_covg_0) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.samples[0]["REF_LIKELIHOOD"] = numeric_limits<float>::lowest();
    vr.samples[0]["ALT_LIKELIHOOD"] = -1.5;
    vr.confidence();
    float exp_confidence = -numeric_limits<float>::lowest()-1.5;
    EXPECT_FLOAT_EQ(exp_confidence,vr.samples[0]["CONFIDENCE"]);
}

TEST(VCFRecordConfidenceTest, handles_alt_covg_0) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.samples[0]["REF_LIKELIHOOD"] = -1.5;
    vr.samples[0]["ALT_LIKELIHOOD"] = numeric_limits<float>::lowest();
    vr.confidence();
    float exp_confidence = -numeric_limits<float>::lowest()-1.5;
    EXPECT_FLOAT_EQ(exp_confidence,vr.samples[0]["CONFIDENCE"]);
}

TEST(VCFRecordSwapTest, swaps_correctly_all_values_present) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.samples[0]["REF_MEAN_FWD_COVG"] = 0;
    vr.samples[0]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[0]["ALT_MEAN_FWD_COVG"] = 2;
    vr.samples[0]["ALT_MEAN_REV_COVG"] = 3;
    vr.samples[0]["REF_MED_FWD_COVG"] = 4;
    vr.samples[0]["REF_MED_REV_COVG"] = 5;
    vr.samples[0]["ALT_MED_FWD_COVG"] = 6;
    vr.samples[0]["ALT_MED_REV_COVG"] = 7;
    vr.samples[0]["REF_SUM_FWD_COVG"] = 8;
    vr.samples[0]["REF_SUM_REV_COVG"] = 9;
    vr.samples[0]["ALT_SUM_FWD_COVG"] = 10;
    vr.samples[0]["ALT_SUM_REV_COVG"] = 11;
    vr.samples[0]["REF_LIKELIHOOD"] = 12;
    vr.samples[0]["ALT_LIKELIHOOD"] = 13;
    vr.samples[0]["GT"] = 14;
    vr.samples[0]["CONFIDENCE"] = 15;
    vr.swap_ref_and_alt_properties(vr.samples[0]);
    EXPECT_EQ(vr.samples[0]["REF_MEAN_FWD_COVG"],2);
    EXPECT_EQ(vr.samples[0]["REF_MEAN_REV_COVG"],3);
    EXPECT_EQ(vr.samples[0]["ALT_MEAN_FWD_COVG"],0);
    EXPECT_EQ(vr.samples[0]["ALT_MEAN_REV_COVG"],1);
    EXPECT_EQ(vr.samples[0]["REF_MED_FWD_COVG"],6);
    EXPECT_EQ(vr.samples[0]["REF_MED_REV_COVG"],7);
    EXPECT_EQ(vr.samples[0]["ALT_MED_FWD_COVG"],4);
    EXPECT_EQ(vr.samples[0]["ALT_MED_REV_COVG"],5);
    EXPECT_EQ(vr.samples[0]["REF_SUM_FWD_COVG"],10);
    EXPECT_EQ(vr.samples[0]["REF_SUM_REV_COVG"],11);
    EXPECT_EQ(vr.samples[0]["ALT_SUM_FWD_COVG"],8);
    EXPECT_EQ(vr.samples[0]["ALT_SUM_REV_COVG"],9);
    EXPECT_EQ(vr.samples[0]["REF_LIKELIHOOD"],13);
    EXPECT_EQ(vr.samples[0]["ALT_LIKELIHOOD"],12);
    EXPECT_EQ(vr.samples[0]["GT"],14);
    EXPECT_EQ(vr.samples[0]["CONFIDENCE"],15);
}

TEST(VCFRecordSwapTest, swaps_correctly_some_values_present) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    vr.samples.push_back(m);
    vr.samples[0]["REF_MEAN_FWD_COVG"] = 0;
    vr.samples[0]["REF_MEAN_REV_COVG"] = 1;

    vr.samples[0]["ALT_MEAN_REV_COVG"] = 3;
    vr.samples[0]["REF_MED_FWD_COVG"] = 4;
    vr.samples[0]["REF_MED_REV_COVG"] = 5;
    vr.samples[0]["ALT_MED_FWD_COVG"] = 6;
    vr.samples[0]["ALT_MED_REV_COVG"] = 7;

    vr.samples[0]["REF_LIKELIHOOD"] = 12;
    vr.samples[0]["ALT_LIKELIHOOD"] = 13;
    vr.samples[0]["GT"] = 14;
    vr.samples[0]["CONFIDENCE"] = 15;
    vr.swap_ref_and_alt_properties(vr.samples[0]);
    EXPECT_EQ(vr.samples[0]["REF_MEAN_FWD_COVG"],'\0');
    EXPECT_EQ(vr.samples[0]["REF_MEAN_REV_COVG"],3);
    EXPECT_EQ(vr.samples[0]["ALT_MEAN_FWD_COVG"],0);
    EXPECT_EQ(vr.samples[0]["ALT_MEAN_REV_COVG"],1);
    EXPECT_EQ(vr.samples[0]["REF_MED_FWD_COVG"],6);
    EXPECT_EQ(vr.samples[0]["REF_MED_REV_COVG"],7);
    EXPECT_EQ(vr.samples[0]["ALT_MED_FWD_COVG"],4);
    EXPECT_EQ(vr.samples[0]["ALT_MED_REV_COVG"],5);
    EXPECT_EQ(vr.samples[0]["REF_SUM_FWD_COVG"],'\0');
    EXPECT_EQ(vr.samples[0]["REF_SUM_REV_COVG"],'\0');
    EXPECT_EQ(vr.samples[0]["ALT_SUM_FWD_COVG"],'\0');
    EXPECT_EQ(vr.samples[0]["ALT_SUM_REV_COVG"],'\0');
    EXPECT_EQ(vr.samples[0]["REF_LIKELIHOOD"],13);
    EXPECT_EQ(vr.samples[0]["ALT_LIKELIHOOD"],12);
    EXPECT_EQ(vr.samples[0]["GT"],14);
    EXPECT_EQ(vr.samples[0]["CONFIDENCE"],15);
}


TEST(VCFRecordTest, equals) {
    VCFRecord vr;
    EXPECT_EQ(vr, vr);

    VCFRecord vr1("chrom1", 3, "A", "T");
    EXPECT_EQ(vr1, vr1);
    EXPECT_EQ((vr == vr1), false);
    EXPECT_EQ((vr1 == vr), false);

    VCFRecord vr2("chrom2", 3, "A", "T");
    EXPECT_EQ(vr2, vr2);
    EXPECT_EQ((vr2 == vr1), false);
    EXPECT_EQ((vr1 == vr2), false);

    VCFRecord vr3("chrom1", 6, "A", "T");
    EXPECT_EQ(vr3, vr3);
    EXPECT_EQ((vr3 == vr1), false);
    EXPECT_EQ((vr1 == vr3), false);

    VCFRecord vr4("chrom1", 3, "G", "T");
    EXPECT_EQ(vr4, vr4);
    EXPECT_EQ((vr4 == vr1), false);
    EXPECT_EQ((vr1 == vr4), false);

    VCFRecord vr5("chrom1", 3, "A", "G");
    EXPECT_EQ(vr5, vr5);
    EXPECT_EQ((vr5 == vr1), false);
    EXPECT_EQ((vr1 == vr5), false);

}

TEST(VCFRecordTest, less_than) {
    VCFRecord vr1("chrom1", 3, "A", "T");
    VCFRecord vr2("chrom2", 3, "A", "T");
    EXPECT_EQ((vr2 < vr1), false);
    EXPECT_EQ((vr1 < vr2), true);

    VCFRecord vr3("chrom1", 6, "A", "T");
    EXPECT_EQ((vr3 < vr1), false);
    EXPECT_EQ((vr1 < vr3), true);

    VCFRecord vr4("chrom1", 3, "G", "T");
    EXPECT_EQ((vr4 < vr1), false);
    EXPECT_EQ((vr1 < vr4), true);

    VCFRecord vr5("chrom1", 3, "A", "G");
    EXPECT_EQ((vr5 < vr1), true);
    EXPECT_EQ((vr1 < vr5), false);
}
