#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "vcfrecord.h"
#include <stdint.h>
#include <iostream>
#include <cmath>

using namespace std;

TEST(VCFRecordTest, create_empty) {

    VCFRecord vr;
    EXPECT_EQ(".", vr.chrom);
    EXPECT_EQ((uint) 0, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ(".", vr.ref);
    EXPECT_EQ(".", vr.alt);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);
    EXPECT_EQ((uint)0, vr.format.size());
    EXPECT_EQ((uint)0, vr.samples.size());
    EXPECT_EQ((uint)0, vr.regt_samples.size());
}

TEST(VCFRecordTest, create_with_values) {

    VCFRecord vr("chrom1", 3, "A", "T");
    EXPECT_EQ("chrom1", vr.chrom);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.ref);
    EXPECT_EQ("T", vr.alt);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ("SVTYPE=SNP", vr.info);
    EXPECT_EQ((uint)1, vr.format.size());
    EXPECT_EQ((uint)0, vr.samples.size());
    EXPECT_EQ((uint)0, vr.regt_samples.size());
}

TEST(VCFRecordTest, create_from_record) {

    VCFRecord template_vr("chrom1", 3, "A", "T");
    VCFRecord vr(template_vr);
    EXPECT_EQ("chrom1", vr.chrom);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.ref);
    EXPECT_EQ("T", vr.alt);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ("SVTYPE=SNP", vr.info);
    EXPECT_EQ((uint)1, vr.format.size());
    EXPECT_EQ((uint)0, vr.samples.size());
    EXPECT_EQ((uint)0, vr.regt_samples.size());
}

TEST(VCFRecordTest, add_formats_none) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vector<string> new_formats = {};
    vr.add_formats(new_formats);
    vector<string> expected_formats = {"GT"};
    EXPECT_ITERABLE_EQ(vector<string>,expected_formats, vr.format);
}

TEST(VCFRecordTest, add_formats_some) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vector<string> new_formats = {"hi", "there"};
    vr.add_formats(new_formats);
    vector<string> expected_formats = {"GT", "hi", "there"};
    EXPECT_ITERABLE_EQ(vector<string>,expected_formats, vr.format);
}

TEST(VCFRecordTest, add_formats_some_repeat) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vector<string> new_formats = {"hi", "there"};
    vr.add_formats(new_formats);
    vr.add_formats(new_formats);
    vector<string> expected_formats = {"GT", "hi", "there"};
    EXPECT_ITERABLE_EQ(vector<string>,expected_formats, vr.format);
}

TEST(VCFRecordTest, add_formats_some_overlapping) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vector<string> new_formats = {"hi", "there"};
    vr.add_formats(new_formats);
    new_formats = {"hi", "again"};
    vr.add_formats(new_formats);
    vector<string> expected_formats = {"GT", "hi", "there", "again"};
    EXPECT_ITERABLE_EQ(vector<string>,expected_formats, vr.format);
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
    bool found_likelihood = vr.regt_samples[0].find("REF_LIKELIHOOD") != vr.regt_samples[0].end();
    EXPECT_FALSE(found_likelihood);

    vr.samples[0]["GT"] = 1;
    vr.likelihood(1,0.01);
    found_likelihood = vr.regt_samples[0].find("REF_LIKELIHOOD") != vr.regt_samples[0].end();
    EXPECT_FALSE(found_likelihood);

    vr.samples[0]["REF_MEAN_FWD_COVG"] = 1;
    vr.samples[0]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[0]["ALT_MEAN_FWD_COVG"] = 1;
    vr.likelihood(1,0.01);
    found_likelihood = vr.regt_samples[0].find("REF_LIKELIHOOD") != vr.regt_samples[0].end();
    EXPECT_FALSE(found_likelihood);

    vr.samples[0].erase("ALT_MEAN_FWD_COVG");
    vr.samples[0]["ALT_MEAN_REV_COVG"] = 1;
    vr.likelihood(1,0.01);
    found_likelihood = vr.regt_samples[0].find("REF_LIKELIHOOD") != vr.regt_samples[0].end();
    EXPECT_FALSE(found_likelihood);

    vr.samples[0].erase("REF_MEAN_FWD_COVG");
    vr.samples[0]["ALT_MEAN_FWD_COVG"] = 1;
    vr.likelihood(1,0.01);
    found_likelihood = vr.regt_samples[0].find("REF_LIKELIHOOD") != vr.regt_samples[0].end();
    EXPECT_FALSE(found_likelihood);

    vr.samples[0].erase("REF_MEAN_REV_COVG");
    vr.samples[0]["REF_MEAN_FWD_COVG"] = 1;
    vr.likelihood(1,0.01);
    found_likelihood = vr.regt_samples[0].find("REF_LIKELIHOOD") != vr.regt_samples[0].end();
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
    bool found_likelihood = vr.regt_samples[0].find("REF_LIKELIHOOD") != vr.regt_samples[0].end();
    EXPECT_TRUE(found_likelihood);
    found_likelihood = vr.regt_samples[0].find("ALT_LIKELIHOOD") != vr.regt_samples[0].end();
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
    EXPECT_FLOAT_EQ(exp_likelihood,vr.regt_samples[0]["REF_LIKELIHOOD"]);
    exp_likelihood = -1-log(4)-log(3)-log(2)+2*log(0.01);
    EXPECT_FLOAT_EQ(exp_likelihood,vr.regt_samples[0]["ALT_LIKELIHOOD"]);
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
    EXPECT_FLOAT_EQ(exp_likelihood,vr.regt_samples[0]["REF_LIKELIHOOD"]);
    exp_likelihood = -1-log(4)-log(3)-log(2);
    EXPECT_FLOAT_EQ(exp_likelihood,vr.regt_samples[0]["ALT_LIKELIHOOD"]);
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
    EXPECT_FLOAT_EQ(exp_likelihood,vr.regt_samples[0]["ALT_LIKELIHOOD"]);
    exp_likelihood = -1-log(2);
    EXPECT_FLOAT_EQ(exp_likelihood,vr.regt_samples[0]["REF_LIKELIHOOD"]);
}

TEST(VCFRecordConfidenceTest, does_not_run_if_info_missing) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, float> m;
    vr.regt_samples.push_back(m);
    vr.confidence();
    bool found_confidence = vr.regt_samples[0].find("CONFIDENCE") != vr.regt_samples[0].end();
    EXPECT_FALSE(found_confidence);

    vr.regt_samples[0]["REF_LIKELIHOOD"] = -1.0;
    vr.confidence();
    found_confidence = vr.regt_samples[0].find("CONFIDENCE") != vr.regt_samples[0].end();
    EXPECT_FALSE(found_confidence);

    vr.regt_samples[0]["ALT_LIKELIHOOD"] = -1.0;
    vr.regt_samples[0].erase("REF_LIKELIHOOD");
    vr.confidence();
    found_confidence = vr.regt_samples[0].find("CONFIDENCE") != vr.regt_samples[0].end();
    EXPECT_FALSE(found_confidence);
}

TEST(VCFRecordConfidenceTest, adds_confidence_with_info) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, float> m;
    m.reserve(3);
    vr.regt_samples.push_back(m);
    vr.regt_samples[0]["REF_LIKELIHOOD"] = -1.0;
    vr.regt_samples[0]["ALT_LIKELIHOOD"] = -2.5;
    vr.confidence();
    bool found_confidence = vr.regt_samples[0].find("CONFIDENCE") != vr.regt_samples[0].end();
    EXPECT_TRUE(found_confidence);
}

TEST(VCFRecordConfidenceTest, gets_correct_confidence_simple_case) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, float> m;
    vr.regt_samples.push_back(m);
    vr.regt_samples[0]["REF_LIKELIHOOD"] = -1.0;
    vr.regt_samples[0]["ALT_LIKELIHOOD"] = 0.0;
    vr.confidence();
    float exp_confidence = 1;
    EXPECT_FLOAT_EQ(exp_confidence,vr.regt_samples[0]["CONFIDENCE"]);
}

TEST(VCFRecordConfidenceTest, handles_ref_covg_0) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, float> m;
    vr.regt_samples.push_back(m);
    vr.regt_samples[0]["REF_LIKELIHOOD"] = numeric_limits<float>::lowest();
    vr.regt_samples[0]["ALT_LIKELIHOOD"] = -1.5;
    vr.confidence();
    float exp_confidence = -numeric_limits<float>::lowest()-1.5;
    EXPECT_FLOAT_EQ(exp_confidence,vr.regt_samples[0]["CONFIDENCE"]);
}

TEST(VCFRecordConfidenceTest, handles_alt_covg_0) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, float> m;
    vr.regt_samples.push_back(m);
    vr.regt_samples[0]["REF_LIKELIHOOD"] = -1.5;
    vr.regt_samples[0]["ALT_LIKELIHOOD"] = numeric_limits<float>::lowest();
    vr.confidence();
    float exp_confidence = -numeric_limits<float>::lowest()-1.5;
    EXPECT_FLOAT_EQ(exp_confidence,vr.regt_samples[0]["CONFIDENCE"]);
}

TEST(VCFRecordRegenotypeTest, correctly_regenotypes) {
    // sample 0 missing confidence
    // sample 1 confidence below threshold
    // sample 2 confidence above threshold, but has correct GT 0 already
    // sample 3 confidence above threshold, but has correct GT 1 already
    // sample 4 confidence above threshold, has incorrect GT 0
    // sample 5 confidence above threshold, has incorrect GT 1

    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    unordered_map<string, float> rm;
    for (uint i = 0; i < 6; ++i) {
        vr.regt_samples.push_back(rm);
        vr.samples.push_back(m);
    }

    vr.samples[0]["REF_MEAN_FWD_COVG"] = 0;
    vr.samples[0]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[0]["ALT_MEAN_FWD_COVG"] = 2;
    vr.samples[0]["ALT_MEAN_REV_COVG"] = 3;
    vr.regt_samples[0]["REF_LIKELIHOOD"] = 4;
    vr.regt_samples[0]["ALT_LIKELIHOOD"] = 5;
    vr.samples[0]["GT"] = 1;

    vr.samples[1]["REF_MEAN_FWD_COVG"] = 0;
    vr.samples[1]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[1]["ALT_MEAN_FWD_COVG"] = 2;
    vr.samples[1]["ALT_MEAN_REV_COVG"] = 3;
    vr.regt_samples[1]["REF_LIKELIHOOD"] = 4;
    vr.regt_samples[1]["ALT_LIKELIHOOD"] = 5;
    vr.samples[1]["GT"] = 1;
    vr.regt_samples[1]["CONFIDENCE"] = 1;

    vr.samples[2]["REF_MEAN_FWD_COVG"] = 0;
    vr.samples[2]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[2]["ALT_MEAN_FWD_COVG"] = 2;
    vr.samples[2]["ALT_MEAN_REV_COVG"] = 3;
    vr.regt_samples[2]["REF_LIKELIHOOD"] = 6;
    vr.regt_samples[2]["ALT_LIKELIHOOD"] = 4;
    vr.samples[2]["GT"] = 0;
    vr.regt_samples[2]["CONFIDENCE"] = 2;

    vr.samples[3]["REF_MEAN_FWD_COVG"] = 0;
    vr.samples[3]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[3]["ALT_MEAN_FWD_COVG"] = 2;
    vr.samples[3]["ALT_MEAN_REV_COVG"] = 3;
    vr.regt_samples[3]["REF_LIKELIHOOD"] = 4;
    vr.regt_samples[3]["ALT_LIKELIHOOD"] = 6;
    vr.samples[3]["GT"] = 1;
    vr.regt_samples[3]["CONFIDENCE"] = 2;

    vr.samples[4]["REF_MEAN_FWD_COVG"] = 0;
    vr.samples[4]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[4]["ALT_MEAN_FWD_COVG"] = 2;
    vr.samples[4]["ALT_MEAN_REV_COVG"] = 3;
    vr.regt_samples[4]["REF_LIKELIHOOD"] = 6;
    vr.regt_samples[4]["ALT_LIKELIHOOD"] = 4;
    vr.samples[4]["GT"] = 1;
    vr.regt_samples[4]["CONFIDENCE"] = 2;

    vr.samples[5]["REF_MEAN_FWD_COVG"] = 0;
    vr.samples[5]["REF_MEAN_REV_COVG"] = 1;
    vr.samples[5]["ALT_MEAN_FWD_COVG"] = 2;
    vr.samples[5]["ALT_MEAN_REV_COVG"] = 3;
    vr.regt_samples[5]["REF_LIKELIHOOD"] = 4;
    vr.regt_samples[5]["ALT_LIKELIHOOD"] = 6;
    vr.samples[5]["GT"] = 0;
    vr.regt_samples[5]["CONFIDENCE"] = 2;

    vr.regenotype(1);
    EXPECT_EQ(vr.samples[0]["REF_MEAN_FWD_COVG"],0);
    EXPECT_EQ(vr.samples[0]["REF_MEAN_REV_COVG"],1);
    EXPECT_EQ(vr.samples[0]["ALT_MEAN_FWD_COVG"],2);
    EXPECT_EQ(vr.samples[0]["ALT_MEAN_REV_COVG"],3);
    EXPECT_EQ(vr.regt_samples[0]["REF_LIKELIHOOD"],4);
    EXPECT_EQ(vr.regt_samples[0]["ALT_LIKELIHOOD"],5);
    EXPECT_EQ(vr.samples[0]["GT"],'\0');

    EXPECT_EQ(vr.samples[1]["REF_MEAN_FWD_COVG"],0);
    EXPECT_EQ(vr.samples[1]["REF_MEAN_REV_COVG"],1);
    EXPECT_EQ(vr.samples[1]["ALT_MEAN_FWD_COVG"],2);
    EXPECT_EQ(vr.samples[1]["ALT_MEAN_REV_COVG"],3);
    EXPECT_EQ(vr.regt_samples[1]["REF_LIKELIHOOD"],4);
    EXPECT_EQ(vr.regt_samples[1]["ALT_LIKELIHOOD"],5);
    EXPECT_EQ(vr.samples[1]["GT"],'\0');

    EXPECT_EQ(vr.samples[2]["REF_MEAN_FWD_COVG"],0);
    EXPECT_EQ(vr.samples[2]["REF_MEAN_REV_COVG"],1);
    EXPECT_EQ(vr.samples[2]["ALT_MEAN_FWD_COVG"],2);
    EXPECT_EQ(vr.samples[2]["ALT_MEAN_REV_COVG"],3);
    EXPECT_EQ(vr.regt_samples[2]["REF_LIKELIHOOD"],6);
    EXPECT_EQ(vr.regt_samples[2]["ALT_LIKELIHOOD"],4);
    EXPECT_EQ(vr.samples[2]["GT"],0);

    EXPECT_EQ(vr.samples[3]["REF_MEAN_FWD_COVG"],0);
    EXPECT_EQ(vr.samples[3]["REF_MEAN_REV_COVG"],1);
    EXPECT_EQ(vr.samples[3]["ALT_MEAN_FWD_COVG"],2);
    EXPECT_EQ(vr.samples[3]["ALT_MEAN_REV_COVG"],3);
    EXPECT_EQ(vr.regt_samples[3]["REF_LIKELIHOOD"],4);
    EXPECT_EQ(vr.regt_samples[3]["ALT_LIKELIHOOD"],6);
    EXPECT_EQ(vr.samples[3]["GT"],1);

    EXPECT_EQ(vr.samples[4]["REF_MEAN_FWD_COVG"],0);
    EXPECT_EQ(vr.samples[4]["REF_MEAN_REV_COVG"],1);
    EXPECT_EQ(vr.samples[4]["ALT_MEAN_FWD_COVG"],2);
    EXPECT_EQ(vr.samples[4]["ALT_MEAN_REV_COVG"],3);
    EXPECT_EQ(vr.regt_samples[4]["REF_LIKELIHOOD"],6);
    EXPECT_EQ(vr.regt_samples[4]["ALT_LIKELIHOOD"],4);
    EXPECT_EQ(vr.samples[4]["GT"],0);

    EXPECT_EQ(vr.samples[5]["REF_MEAN_FWD_COVG"],0);
    EXPECT_EQ(vr.samples[5]["REF_MEAN_REV_COVG"],1);
    EXPECT_EQ(vr.samples[5]["ALT_MEAN_FWD_COVG"],2);
    EXPECT_EQ(vr.samples[5]["ALT_MEAN_REV_COVG"],3);
    EXPECT_EQ(vr.regt_samples[5]["REF_LIKELIHOOD"],4);
    EXPECT_EQ(vr.regt_samples[5]["ALT_LIKELIHOOD"],6);
    EXPECT_EQ(vr.samples[5]["GT"],1);
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

TEST(VCFRecordTest, ostream) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vector<string> v = {"chrom1","3",".","A","T",".",".","SVTYPE=SNP","GT"};
    stringstream out;
    out << vr;
    string rr;
    for (auto  s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
}

TEST(VCFRecordTest, ostream_with_sample_not_all_info_in_formats) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    m["GT"] = 1;
    m["pringle"] = 2;
    vr.samples.push_back(m);
    vector<string> v = {"chrom1","3",".","A","T",".",".","SVTYPE=SNP","GT"};
    stringstream out;
    out << vr;
    string rr;
    for (auto  s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
    uint8_t u=1;
    uint ru;
    out >> ru;
    EXPECT_EQ(u, ru);
}

TEST(VCFRecordTest, ostream_with_sample_including_all_formats) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    m["GT"] = 0;
    m["pringle"] = 2;
    vr.samples.push_back(m);
    vr.add_formats({"pringle"});
    vector<string> v = {"chrom1","3",".","A","T",".",".","SVTYPE=SNP","GT:pringle"};
    vector<uint8_t> vu = {0,2};
    stringstream out;
    out << vr;
    string rr;
    for (auto  s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
    uint ru;
    for (auto  s : vu) {
        out >> ru;
        EXPECT_EQ(s, ru);
        out.ignore(1, ':');
    }

}

TEST(VCFRecordTest, ostream_with_sample_more_formats_than_info) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    m["GT"] = 0;
    vr.samples.push_back(m);
    vr.add_formats({"pringle"});
    vector<string> v = {"chrom1","3",".","A","T",".",".","SVTYPE=SNP","GT:pringle"};
    vector<uint8_t> vu = {0};
    stringstream out;
    out << vr;
    string rr;
    for (auto  s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
    uint ru;
    uint8_t u=0;
    out >> ru;
    EXPECT_EQ(u, ru);
    out >> rr;
    EXPECT_EQ(":.", rr);
}

TEST(VCFRecordTest, ostream_with_sample_more_formats_than_info_regt) {
    VCFRecord vr("chrom1", 3, "A", "T");
    unordered_map<string, uint8_t> m;
    m["GT"] = 0;
    vr.samples.push_back(m);
    unordered_map<string, float> n;
    n["pringle"] = 0.1;
    vr.regt_samples.push_back(n);
    vr.add_formats({"pringle"});
    vector<string> v = {"chrom1","3",".","A","T",".",".","SVTYPE=SNP","GT:pringle"};
    stringstream out;
    out << vr;
    string rr;
    for (auto  s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
    uint ru;
    uint8_t u = 0;
    out >> ru;
    EXPECT_EQ(u, ru);
    out.ignore(1,':');
    float rf=0.0,f=0.1;
    out >> rf;
    EXPECT_EQ(f,rf);
}


