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
    EXPECT_EQ((uint) 0, vr.alts.size());
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);
    EXPECT_EQ((uint) 0, vr.format.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());
}

TEST(VCFRecordTest, create_with_values) {

    VCFRecord vr("chrom1", 3, "A", "T");
    EXPECT_EQ("chrom1", vr.chrom);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.ref);
    EXPECT_EQ("T", vr.alts[0]);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ("SVTYPE=SNP", vr.info);
    EXPECT_EQ((uint) 1, vr.format.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());
}

TEST(VCFRecordTest, create_from_record) {

    VCFRecord template_vr("chrom1", 3, "A", "T");
    VCFRecord vr(template_vr);
    EXPECT_EQ("chrom1", vr.chrom);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.ref);
    EXPECT_EQ("T", vr.alts[0]);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ("SVTYPE=SNP", vr.info);
    EXPECT_EQ((uint) 1, vr.format.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());
}

TEST(VCFRecordTest, create_from_record_with_samples) {

    VCFRecord template_vr("chrom1", 3, "A", "T");
    template_vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(2);
    template_vr.sampleIndex_to_format_to_sampleInfo[0]["GT"] = {0};
    template_vr.sampleIndex_to_format_to_sampleInfo[1]["GT"] = {1};
    template_vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(2);
    template_vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"] = {4.0};
    template_vr.sampleIndex_to_format_to_sampleGenotypedInfo[1]["GT_CONF"] = {6.3};
    VCFRecord vr(template_vr);
    EXPECT_EQ("chrom1", vr.chrom);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.ref);
    EXPECT_EQ("T", vr.alts[0]);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ("SVTYPE=SNP", vr.info);
    EXPECT_EQ((uint) 1, vr.format.size());
    EXPECT_EQ((uint) 2, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 2, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());
}

TEST(VCFRecordTest, clear_simple) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.clear();
    EXPECT_EQ(".", vr.chrom);
    EXPECT_EQ((uint) 0, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ(".", vr.ref);
    EXPECT_EQ((uint) 0, vr.alts.size());
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);
    EXPECT_EQ((uint) 0, vr.format.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());
}

TEST(VCFRecordTest, clear_with_samples) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(2);
    vr.sampleIndex_to_format_to_sampleInfo[0]["GT"] = {0};
    vr.sampleIndex_to_format_to_sampleInfo[1]["GT"] = {1};
    vr.clear();
    EXPECT_EQ(".", vr.chrom);
    EXPECT_EQ((uint) 0, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ(".", vr.ref);
    EXPECT_EQ((uint) 0, vr.alts.size());
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);
    EXPECT_EQ((uint) 0, vr.format.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());
}

TEST(VCFRecordTest, clear_with_sample_confs) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(2);
    vr.sampleIndex_to_format_to_sampleInfo[0]["GT"] = {0};
    vr.sampleIndex_to_format_to_sampleInfo[1]["GT"] = {1};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(2);
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {0, 4, 5};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[1]["LIKELIHOOD"] = {1, 2};
    vr.add_formats({"GT", "LIKELIHOOD"});
    vr.clear();
    EXPECT_EQ(".", vr.chrom);
    EXPECT_EQ((uint) 0, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ(".", vr.ref);
    EXPECT_EQ((uint) 0, vr.alts.size());
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);
    EXPECT_EQ((uint) 0, vr.format.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());
}

TEST(VCFRecordTest, clear_sample) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(2);
    vr.sampleIndex_to_format_to_sampleInfo[0]["GT"] = {0};
    vr.sampleIndex_to_format_to_sampleInfo[1]["GT"] = {1};

    // index out of range, nothin happens
    vr.clear_sample(3);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ("A", vr.ref);
    EXPECT_EQ((uint) 1, vr.alts.size());
    EXPECT_EQ("T", vr.alts[0]);
    EXPECT_EQ((uint) 2, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 1, vr.sampleIndex_to_format_to_sampleInfo[0].size());
    EXPECT_EQ((uint) 1, vr.sampleIndex_to_format_to_sampleInfo[1].size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());

    // index in range
    vr.clear_sample(1);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ("A", vr.ref);
    EXPECT_EQ((uint) 1, vr.alts.size());
    EXPECT_EQ("T", vr.alts[0]);
    EXPECT_EQ((uint) 2, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 1, vr.sampleIndex_to_format_to_sampleInfo[0].size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleInfo[1].size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());

    // index in range, nothing to do
    vr.clear_sample(1);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ("A", vr.ref);
    EXPECT_EQ((uint) 1, vr.alts.size());
    EXPECT_EQ("T", vr.alts[0]);
    EXPECT_EQ((uint) 2, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 1, vr.sampleIndex_to_format_to_sampleInfo[0].size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleInfo[1].size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());

    // index in range, last sample
    vr.clear_sample(0);
    EXPECT_EQ((uint) 0, vr.pos);
    EXPECT_EQ(".", vr.ref);
    EXPECT_EQ((uint) 0, vr.alts.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleInfo.size());
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_format_to_sampleGenotypedInfo.size());
}

TEST(VCFRecordTest, add_formats_none) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vector<string> new_formats = {};
    vr.add_formats(new_formats);
    vector<string> expected_formats = {"GT"};
    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
}

TEST(VCFRecordTest, add_formats_some) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vector<string> new_formats = {"hi", "there"};
    vr.add_formats(new_formats);
    vector<string> expected_formats = {"GT", "hi", "there"};
    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
}

TEST(VCFRecordTest, add_formats_some_repeat) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vector<string> new_formats = {"hi", "there"};
    vr.add_formats(new_formats);
    vr.add_formats(new_formats);
    vector<string> expected_formats = {"GT", "hi", "there"};
    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
}

TEST(VCFRecordTest, add_formats_some_overlapping) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vector<string> new_formats = {"hi", "there"};
    vr.add_formats(new_formats);
    new_formats = {"hi", "again"};
    vr.add_formats(new_formats);
    vector<string> expected_formats = {"GT", "hi", "there", "again"};
    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
}

TEST(VCFRecordTest, add_format_death_no_samples) {
    VCFRecord vr("chrom1", 3, "A", "T");
    uint16_t v = 20;
    EXPECT_DEATH(vr.set_format(0, "hello", v), "");
    float w = 20.0;
    EXPECT_DEATH(vr.set_format(0, "hello", w), "");
}

TEST(VCFRecordTest, add_format_cap_too_big) {
    VCFRecord vr("chrom1", 3, "A", "T");
    uint32_t v = 60000000;
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.set_format(0, "hello", v);
    EXPECT_EQ(vr.get_format_u(0, "hello")[0], std::numeric_limits<uint16_t>::max() - 1);
}

TEST(VCFRecordTest, add_format_new_uint) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    uint16_t v = 20;
    vr.set_format(0, "hello", v);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo.size(), 1);
    EXPECT_TRUE(vr.sampleIndex_to_format_to_sampleInfo[0].find("hello") != vr.sampleIndex_to_format_to_sampleInfo[0].end());
    std::vector<uint16_t> exp_v = {v};
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, vr.sampleIndex_to_format_to_sampleInfo[0]["hello"], exp_v);
    std::vector<std::string> exp_f = {"GT", "hello"};
    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
}

TEST(VCFRecordTest, add_format_old_uint_overwritten) {
    VCFRecord vr("chrom1", 3, "A", "T");
    SampleInfo<uint16_t> sampleInfo;
    sampleInfo["hello"] = {10};
    vr.sampleIndex_to_format_to_sampleInfo.push_back(sampleInfo);
    uint16_t v = 20;
    vr.set_format(0, "hello", v);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo.size(), 1);
    EXPECT_TRUE(vr.sampleIndex_to_format_to_sampleInfo[0].find("hello") != vr.sampleIndex_to_format_to_sampleInfo[0].end());
    std::vector<uint16_t> exp_v = {v};
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, vr.sampleIndex_to_format_to_sampleInfo[0]["hello"], exp_v);
    std::vector<std::string> exp_f = {"GT", "hello"};
    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
}

TEST(VCFRecordTest, add_format_new_float) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    float v = 20.0;
    vr.set_format(0, "hello", v);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo.size(), 1);
    EXPECT_TRUE(vr.sampleIndex_to_format_to_sampleGenotypedInfo[0].find("hello") != vr.sampleIndex_to_format_to_sampleGenotypedInfo[0].end());
    std::vector<float> exp_v = {v};
    EXPECT_ITERABLE_EQ(std::vector<float>, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["hello"], exp_v);
    std::vector<std::string> exp_f = {"GT", "hello"};
    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
}

TEST(VCFRecordTest, add_format_old_float_overwritten) {
    VCFRecord vr("chrom1", 3, "A", "T");
    SampleInfo<float> sample_info;
    sample_info["hello"].push_back(10.0);
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back(sample_info);
    float v = 20.0;
    vr.set_format(0, "hello", v);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo.size(), 1);
    EXPECT_TRUE(vr.sampleIndex_to_format_to_sampleGenotypedInfo[0].find("hello") != vr.sampleIndex_to_format_to_sampleGenotypedInfo[0].end());
    std::vector<float> exp_v = {v};
    EXPECT_ITERABLE_EQ(std::vector<float>, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["hello"], exp_v);
    std::vector<std::string> exp_f = {"GT", "hello"};
    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
}

TEST(VCFRecordTest, append_format_old_uint) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    uint16_t v = 10;
    vr.set_format(0, "hello", v);
    v = 20;
    vr.append_format(0, "hello", v);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo.size(), 1);
    EXPECT_TRUE(vr.sampleIndex_to_format_to_sampleInfo[0].find("hello") != vr.sampleIndex_to_format_to_sampleInfo[0].end());
    std::vector<uint16_t> exp_v = {10, 20};
    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, vr.sampleIndex_to_format_to_sampleInfo[0]["hello"], exp_v);
    std::vector<std::string> exp_f = {"GT", "hello"};
    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
}

TEST(VCFRecordTest, append_format_old_float) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
    float v = 10.0;
    vr.set_format(0, "hello", v);
    v = 20.0;
    vr.append_format(0, "hello", v);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo.size(), 1);
    EXPECT_TRUE(vr.sampleIndex_to_format_to_sampleGenotypedInfo[0].find("hello") != vr.sampleIndex_to_format_to_sampleGenotypedInfo[0].end());
    std::vector<float> exp_v = {10.0, 20.0};
    EXPECT_ITERABLE_EQ(std::vector<float>, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["hello"], exp_v);
    std::vector<std::string> exp_f = {"GT", "hello"};
    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
}

TEST(VCFRecordTest, get_format_float) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
    float v = 10.0;
    vr.set_format(0, "hello", v);

    auto res = vr.get_format_f(1,"hello");
    EXPECT_EQ(res.size(), 0);
    EXPECT_TRUE(res.empty());

    res = vr.get_format_f(0,"help");
    EXPECT_EQ(res.size(), 0);
    EXPECT_TRUE(res.empty());

    res = vr.get_format_f(0,"hello");
    EXPECT_EQ(res.size(), 1);
    EXPECT_FALSE(res.empty());
}

TEST(VCFRecordTest, get_format_uint) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    uint16_t v = 10;
    vr.set_format(0, "hello", v);

    auto res = vr.get_format_u(1,"hello");
    EXPECT_EQ(res.size(), 0);
    EXPECT_TRUE(res.empty());

    res = vr.get_format_u(0,"help");
    EXPECT_EQ(res.size(), 0);
    EXPECT_TRUE(res.empty());

    res = vr.get_format_u(0,"hello");
    EXPECT_EQ(res.size(), 1);
    EXPECT_FALSE(res.empty());
}

TEST(VCFRecordLikelihoodTest, does_not_crash_with_no_samples) {
    VCFRecord vr("chrom1", 3, "A", "T");
    EXPECT_NO_FATAL_FAILURE(vr.likelihood({}, 0.01, 0));
}

TEST(VCFRecordLikelihoodTest, does_not_run_if_info_missing) {
    VCFRecord vr("chrom1", 3, "A", "T");
    SampleInfo<uint16_t> sample_info;
    sample_info["nothing"] = {0};
    vr.sampleIndex_to_format_to_sampleInfo.push_back(sample_info);
    assert(vr.sampleIndex_to_format_to_sampleInfo.size() > 0);
    std::vector<float> f = {0.0, 0.0};
    vr.set_format(0,"GAPS", f);
    vr.likelihood({1}, 0.01, 0);
    bool found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
    EXPECT_FALSE(found_likelihood);

    vr.sampleIndex_to_format_to_sampleInfo[0]["GT"] = {1};
    vr.likelihood({1}, 0.01, 0);
    found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
    EXPECT_FALSE(found_likelihood);

    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 1};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1};
    vr.likelihood({1}, 0.01, 0);
    found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
    EXPECT_FALSE(found_likelihood);

    vr.sampleIndex_to_format_to_sampleInfo[0].erase("MEAN_FWD_COVG");
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 1};
    vr.likelihood({1}, 0.01, 0);
    found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
    EXPECT_FALSE(found_likelihood);

    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1};
    vr.likelihood({1}, 0.01, 0);
    found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
    EXPECT_FALSE(found_likelihood);

    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 1};
    vr.sampleIndex_to_format_to_sampleInfo[0].erase("MEAN_REV_COVG");
    vr.likelihood({1}, 0.01, 0);
    found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
    EXPECT_FALSE(found_likelihood);
}

TEST(VCFRecordLikelihoodTest, adds_likelihood_with_info) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 2};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 2};
    std::vector<float> f = {0.0, 0.0};
    vr.set_format(0,"GAPS", f);
    vr.likelihood({1}, 0.01, 0);
    bool found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
    EXPECT_TRUE(found_likelihood);
}

TEST(VCFRecordLikelihoodTest, gets_correct_likelihood_simple_case) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 2};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 2};
    std::vector<float> f = {0.0, 0.0};
    vr.set_format(0,"GAPS", f);
    vr.likelihood({1}, 0.01, 0);
    float exp_likelihood = -1 - log(2) + 4 * log(0.01) + log(1-exp(-(float(1))));
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][0]);
    exp_likelihood = -1 - log(4) - log(3) - log(2) + 2 * log(0.01) + log(1-exp(-(float(1))));
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][1]);
}

TEST(VCFRecordLikelihoodTest, gets_correct_likelihood_with_min_covg_threshold) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 2};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 2};
    std::vector<float> f = {0.0, 0.0};
    vr.set_format(0,"GAPS", f);
    vr.likelihood({1}, 0.01, 3);

    float exp_likelihood = 4 * log(0.01) - 1 + log(1-exp(-(float(1))));
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][0]);
    exp_likelihood = - 1 - log(4) - log(3) - log(2) + log(1-exp(-(float(1))));
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][1]);
}

TEST(VCFRecordLikelihoodTest, handles_ref_covg_0) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 2};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {0, 2};
    std::vector<float> f = {0.0, 0.0};
    vr.set_format(0,"GAPS", f);
    vr.likelihood({1}, 0.01, 0);
    float exp_likelihood = -1 + 4 * log(0.01) + log(1-exp(-(float(1))));
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][0]);
    exp_likelihood = -1 - log(4) - log(3) - log(2) + log(1-exp(-(float(1))));
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][1]);
}

TEST(VCFRecordLikelihoodTest, handles_alt_covg_0) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 0};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 0};
    std::vector<float> f = {0.0, 0.0};
    vr.set_format(0,"GAPS", f);
    vr.likelihood({1}, 0.01, 0);
    float exp_likelihood = -1 + 2 * log(0.01) + log(1-exp(-(float(1))));
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][1]);
    exp_likelihood = -1 - log(2) + log(1-exp(-(float(1))));
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][0]);
}

TEST(VCFRecordLikelihoodTest, gets_correct_likelihood_gaps) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 2};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 2};
    std::vector<float> f = {0.5, 0.8};
    vr.set_format(0,"GAPS", f);
    vr.likelihood({1}, 0.01, 0);
    float exp_likelihood = -1 - log(2) + 4 * log(0.01) + 0.5*log(1-exp(-(float(1)))) - 0.5;
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][0]);
    exp_likelihood = -1 - log(4) - log(3) - log(2) + 2 * log(0.01) + 0.2*log(1-exp(-(float(1)))) - 0.8;
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][1]);
}

TEST(VCFRecordLikelihoodTest, death_not_enough_covgs) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(2);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 2};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 2};
    vr.sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"] = {1, 2};
    vr.sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"] = {1, 2};
    std::vector<float> f = {0.5, 0.8};
    vr.set_format(0,"GAPS", f);
    vr.set_format(1,"GAPS", f);
    EXPECT_DEATH(vr.likelihood({1}, 0.01, 0), "");
}

TEST(VCFRecordLikelihoodTest, samples_with_different_depths) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(2);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 2};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 2};
    vr.sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"] = {1, 2};
    vr.sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"] = {1, 2};
    std::vector<float> f = {0.5, 0.8};
    vr.set_format(0,"GAPS", f);
    vr.set_format(1,"GAPS", f);
    vr.likelihood({1,2}, 0.01, 0);

    float exp_likelihood = -1 - log(2) + 4 * log(0.01) + 0.5*log(1-exp(-(float(1)))) - 0.5;
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][0]);
    exp_likelihood = -1 - log(4) - log(3) - log(2) + 2 * log(0.01) + 0.2*log(1-exp(-(float(1)))) - 0.8;
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][1]);
    exp_likelihood = 2*log(2) -2 - log(2) + 4 * log(0.01) + 0.5*log(1-exp(-(float(2)))) - 2*0.5;
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[1]["LIKELIHOOD"][0]);
    exp_likelihood = 4*log(2) -2 - log(4) - log(3) - log(2) + 2 * log(0.01) + 0.2*log(1-exp(-(float(2)))) - 2*0.8;
    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_format_to_sampleGenotypedInfo[1]["LIKELIHOOD"][1]);
}

TEST(VCFRecordConfidenceTest, does_not_run_if_info_missing) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
    vr.confidence();
    bool found_confidence = !vr.get_format_f(0,"GT_CONF").empty();
    EXPECT_FALSE(found_confidence);
    std::vector<float> f = {-1.0};
    vr.set_format(0,"LIKELIHOOD", f);
    EXPECT_DEATH(vr.confidence(), "");
}

TEST(VCFRecordConfidenceTest, adds_confidence_with_info) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-1.0, -2.5};
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 0};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {0, 0};
    vr.confidence();
    bool found_confidence = !vr.get_format_f(0,"GT_CONF").empty();
    EXPECT_TRUE(found_confidence);
}

TEST(VCFRecordConfidenceTest, gets_correct_confidence_simple_case) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-1.0, 0.0};
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 0};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {0, 0};
    vr.confidence();
    float exp_confidence = 1;
    EXPECT_FLOAT_EQ(exp_confidence, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"][0]);
}

TEST(VCFRecordConfidenceTest, gets_correct_confidence_two_alts) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.alts.push_back("C");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-14.0, -6.0, -3.0};
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 0, 0};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {0, 0, 0};
    vr.confidence();
    float exp_confidence = 3;
    EXPECT_FLOAT_EQ(exp_confidence, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"][0]);
}

TEST(VCFRecordConfidenceTest, gets_correct_confidence_min_total) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.alts.push_back("C");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-14.0, -6.0, -3.0};
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 0, 1};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {0, 0, 1};
    vr.confidence(3,0);
    float exp_confidence = 0;
    EXPECT_FLOAT_EQ(exp_confidence, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"][0]);
    vr.confidence(2,0);
    exp_confidence = 3;
    EXPECT_FLOAT_EQ(exp_confidence, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"][0]);
}

TEST(VCFRecordConfidenceTest, gets_correct_confidence_min_diff) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.alts.push_back("C");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-14.0, -6.0, -3.0};
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 2, 4};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {0, 0, 1};
    vr.confidence(0,4);
    float exp_confidence = 0;
    EXPECT_FLOAT_EQ(exp_confidence, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"][0]);
    vr.confidence(0,3);
    exp_confidence = 3;
    EXPECT_FLOAT_EQ(exp_confidence, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"][0]);
}

TEST(VCFRecordConfidenceTest, handles_ref_covg_0) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {std::numeric_limits<float>::lowest(), -1.5};
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 0};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {0, 0};
    vr.confidence();
    float exp_confidence = -std::numeric_limits<float>::lowest() - 1.5;
    EXPECT_FLOAT_EQ(exp_confidence, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"][0]);
}

TEST(VCFRecordConfidenceTest, handles_alt_covg_0) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {-1.5, std::numeric_limits<float>::lowest()};
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(1);
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 0};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {0, 0};
    vr.confidence();
    float exp_confidence = -std::numeric_limits<float>::lowest() - 1.5;
    EXPECT_FLOAT_EQ(exp_confidence, vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["GT_CONF"][0]);
}

TEST(VCFRecordRegenotypeTest, correctly_genotypes) {
    // sample 0 missing confidence
    // sample 1 confidence below threshold
    // sample 2 confidence above threshold, but has correct GT 0 already
    // sample 3 confidence above threshold, but has correct GT 1 already
    // sample 4 confidence above threshold, has incorrect GT 0
    // sample 5 confidence above threshold, has incorrect GT 1

    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back_several_empty_sample_infos(6);
    vr.sampleIndex_to_format_to_sampleInfo.push_back_several_empty_sample_infos(6);

    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 2};
    vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 3};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"] = {4, 5};
    vr.sampleIndex_to_format_to_sampleInfo[0]["GT"] = {1};

    vr.sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"] = {0, 2};
    vr.sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"] = {1, 3};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[1]["LIKELIHOOD"] = {4, 5};
    vr.sampleIndex_to_format_to_sampleInfo[1]["GT"] = {1};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[1]["GT_CONF"] = {1};

    vr.sampleIndex_to_format_to_sampleInfo[2]["MEAN_FWD_COVG"] = {0, 2};
    vr.sampleIndex_to_format_to_sampleInfo[2]["MEAN_REV_COVG"] = {1, 3};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[2]["LIKELIHOOD"] = {6, 4};
    vr.sampleIndex_to_format_to_sampleInfo[2]["GT"] = {0};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[2]["GT_CONF"] = {2};

    vr.sampleIndex_to_format_to_sampleInfo[3]["MEAN_FWD_COVG"] = {0, 2};
    vr.sampleIndex_to_format_to_sampleInfo[3]["MEAN_REV_COVG"] = {1, 3};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[3]["LIKELIHOOD"] = {4, 6};
    vr.sampleIndex_to_format_to_sampleInfo[3]["GT"] = {1};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[3]["GT_CONF"] = {2};

    vr.sampleIndex_to_format_to_sampleInfo[4]["MEAN_FWD_COVG"] = {0, 2};
    vr.sampleIndex_to_format_to_sampleInfo[4]["MEAN_REV_COVG"] = {1, 3};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[4]["LIKELIHOOD"] = {6, 4};
    vr.sampleIndex_to_format_to_sampleInfo[4]["GT"] = {1};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[4]["GT_CONF"] = {2};

    vr.sampleIndex_to_format_to_sampleInfo[5]["MEAN_FWD_COVG"] = {0, 2};
    vr.sampleIndex_to_format_to_sampleInfo[5]["MEAN_REV_COVG"] = {1, 3};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[5]["LIKELIHOOD"] = {4, 6};
    vr.sampleIndex_to_format_to_sampleInfo[5]["GT"] = {0};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo[5]["GT_CONF"] = {2};

    vr.genotype(1);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"][0], 0);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"][0], 1);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_FWD_COVG"][1], 2);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[0]["MEAN_REV_COVG"][1], 3);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][0], 4);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[0]["LIKELIHOOD"][1], 5);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[0]["GT"].size(), (uint) 0);

    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"][0], 0);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"][0], 1);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[1]["MEAN_FWD_COVG"][1], 2);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[1]["MEAN_REV_COVG"][1], 3);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[1]["LIKELIHOOD"][0], 4);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[1]["LIKELIHOOD"][1], 5);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[1]["GT"].size(), (uint) 0);

    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[2]["MEAN_FWD_COVG"][0], 0);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[2]["MEAN_REV_COVG"][0], 1);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[2]["MEAN_FWD_COVG"][1], 2);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[2]["MEAN_REV_COVG"][1], 3);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[2]["LIKELIHOOD"][0], 6);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[2]["LIKELIHOOD"][1], 4);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[2]["GT"][0], 0);

    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[3]["MEAN_FWD_COVG"][0], 0);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[3]["MEAN_REV_COVG"][0], 1);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[3]["MEAN_FWD_COVG"][1], 2);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[3]["MEAN_REV_COVG"][1], 3);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[3]["LIKELIHOOD"][0], 4);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[3]["LIKELIHOOD"][1], 6);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[3]["GT"][0], 1);

    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[4]["MEAN_FWD_COVG"][0], 0);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[4]["MEAN_REV_COVG"][0], 1);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[4]["MEAN_FWD_COVG"][1], 2);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[4]["MEAN_REV_COVG"][1], 3);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[4]["LIKELIHOOD"][0], 6);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[4]["LIKELIHOOD"][1], 4);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[4]["GT"][0], 0);

    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[5]["MEAN_FWD_COVG"][0], 0);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[5]["MEAN_REV_COVG"][0], 1);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[5]["MEAN_FWD_COVG"][1], 2);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[5]["MEAN_REV_COVG"][1], 3);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[5]["LIKELIHOOD"][0], 4);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleGenotypedInfo[5]["LIKELIHOOD"][1], 6);
    EXPECT_EQ(vr.sampleIndex_to_format_to_sampleInfo[5]["GT"][0], 1);
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
    vector<string> v = {"chrom1", "4", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT"};
    stringstream out;
    out << vr;
    string rr;
    for (const auto &s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
}

TEST(VCFRecordTest, ostream_with_sample_not_all_info_in_formats) {
    VCFRecord vr("chrom1", 3, "A", "T");
    SampleInfo<uint16_t> sample_info;
    sample_info["GT"] = {1};
    sample_info["pringle"] = {2};
    vr.sampleIndex_to_format_to_sampleInfo.push_back(sample_info);
    vector<string> v = {"chrom1", "4", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT"};
    stringstream out;
    out << vr;
    string rr;
    for (const auto &s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
    uint16_t u = 1;
    uint ru;
    out >> ru;
    EXPECT_EQ(u, ru);
}

TEST(VCFRecordTest, ostream_with_sample_including_all_formats) {
    VCFRecord vr("chrom1", 3, "A", "T");
    SampleInfo<uint16_t> sample_info;
    sample_info["GT"] = {0};
    sample_info["pringle"] = {2};
    vr.sampleIndex_to_format_to_sampleInfo.push_back(sample_info);
    vr.add_formats({"pringle"});
    vector<string> v = {"chrom1", "4", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT:pringle"};
    vector<uint16_t> vu = {0, 2};
    stringstream out;
    out << vr;
    string rr;
    for (const auto &s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
    uint ru;
    for (const auto &s : vu) {
        out >> ru;
        EXPECT_EQ(s, ru);
        out.ignore(1, ':');
    }

}

TEST(VCFRecordTest, ostream_with_sample_more_formats_than_info) {
    VCFRecord vr("chrom1", 3, "A", "T");
    SampleInfo<uint16_t> sample_info;
    sample_info["GT"] = {0};
    vr.sampleIndex_to_format_to_sampleInfo.push_back(sample_info);
    vr.add_formats({"pringle"});
    vector<string> v = {"chrom1", "4", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT:pringle"};
    vector<uint16_t> vu = {0};
    stringstream out;
    out << vr;
    string rr;
    for (const auto &s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
    uint ru;
    uint16_t u = 0;
    out >> ru;
    EXPECT_EQ(u, ru);
    out >> rr;
    EXPECT_EQ(":.", rr);
}

TEST(VCFRecordTest, ostream_with_sample_more_formats_than_info_regt) {
    VCFRecord vr("chrom1", 3, "A", "T");
    SampleInfo<uint16_t> sample_info_int;
    sample_info_int["GT"] = {0};
    vr.sampleIndex_to_format_to_sampleInfo.push_back(sample_info_int);
    SampleInfo<float> sample_info_float;
    sample_info_float["pringle"] = {0.1};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back(sample_info_float);
    vr.add_formats({"pringle"});
    vector<string> v = {"chrom1", "4", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT:pringle"};
    stringstream out;
    out << vr;
    string rr;
    for (const auto &s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
    uint ru;
    uint16_t u = 0;
    out >> ru;
    EXPECT_EQ(u, ru);
    out.ignore(1, ':');
    float rf = 0.0, f = 0.1;
    out >> rf;
    EXPECT_EQ(f, rf);
}

TEST(VCFRecordTest, ostream_with_zero_pos) {
    VCFRecord vr("chrom1", 0, "A", "T");
    SampleInfo<uint16_t> sample_info_int;
    sample_info_int["GT"] = {0};
    vr.sampleIndex_to_format_to_sampleInfo.push_back(sample_info_int);
    SampleInfo<float> sample_info_float;
    sample_info_float["pringle"] = {0.1};
    vr.sampleIndex_to_format_to_sampleGenotypedInfo.push_back(sample_info_float);
    vr.add_formats({"pringle"});
    vector<string> v = {"chrom1", "1", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT:pringle"};
    stringstream out;
    out << vr;
    string rr;
    for (const auto &s : v) {
        out >> rr;
        EXPECT_EQ(s, rr);
    }
    uint ru;
    uint16_t u = 0;
    out >> ru;
    EXPECT_EQ(u, ru);
    out.ignore(1, ':');
    float rf = 0.0, f = 0.1;
    out >> rf;
    EXPECT_EQ(f, rf);
}


TEST(VCFRecordTest, get_longest_allele_length___longest_allele_is_ref) {
    VCFRecord vr("chrom1", 1, "ACGT", "A", ".", ".");

    size_t actual = vr.get_longest_allele_length();
    size_t expected = 4;

    EXPECT_EQ(actual, expected);
}

TEST(VCFRecordTest, get_longest_allele_length___longest_allele_is_alt) {
    VCFRecord vr("chrom1", 1, "ACGT", "A", ".", ".");
    vr.alts.push_back("ACGTTTTAC");
    vr.alts.push_back("C");

    size_t actual = vr.get_longest_allele_length();
    size_t expected = 9;

    EXPECT_EQ(actual, expected);
}