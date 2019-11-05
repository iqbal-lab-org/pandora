#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "test_macro.cpp"
#include "vcfrecord.h"
#include <stdint.h>
#include <iostream>
#include <cmath>
#include "test_helpers.h"


using namespace std;
using ::testing::Return;

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
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_sampleInfo.size());
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
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_sampleInfo.size());
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
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_sampleInfo.size());
}

TEST(VCFRecordTest, create_from_record_with_samples) {

    VCFRecord template_vr("chrom1", 3, "A", "T");
    template_vr.sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(2, &default_genotyping_options);
    VCFRecord vr(template_vr);
    EXPECT_EQ("chrom1", vr.chrom);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.ref);
    EXPECT_EQ("T", vr.alts[0]);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ("SVTYPE=SNP", vr.info);
    EXPECT_EQ((uint) 2, vr.sampleIndex_to_sampleInfo.size());
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
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_sampleInfo.size());
}

TEST(VCFRecordTest, clear_with_samples) {
    VCFRecord vr("chrom1", 3, "A", "T");
    vr.sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(2, &default_genotyping_options);
    vr.clear();
    EXPECT_EQ(".", vr.chrom);
    EXPECT_EQ((uint) 0, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ(".", vr.ref);
    EXPECT_EQ((uint) 0, vr.alts.size());
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);
    EXPECT_EQ((uint) 0, vr.sampleIndex_to_sampleInfo.size());
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


TEST(VCFRecordTest, alts_to_string___no_alts) {
    VCFRecord vcf_record;

    std::string actual = vcf_record.alts_to_string();

    EXPECT_EQ(".", actual);
}

TEST(VCFRecordTest, alts_to_string___three_alts) {
    VCFRecord vcf_record;
    vcf_record.alts.push_back("A1");
    vcf_record.alts.push_back("A2");
    vcf_record.alts.push_back("A3");

    std::string actual = vcf_record.alts_to_string();

    EXPECT_EQ("A1,A2,A3", actual);
}


TEST(VCFRecordTest, get_format___no_flags_set___expects_death) {
    VCFRecord vcf_record;
    EXPECT_DEATH(vcf_record.get_format(false, false), "");
}

TEST(VCFRecordTest, get_format___both_flags_set___expects_death) {
    VCFRecord vcf_record;
    EXPECT_DEATH(vcf_record.get_format(true, true), "");
}

TEST(VCFRecordTest, get_format___genotyping_from_maximum_likelihood) {
    VCFRecord vcf_record;

    std::string actual = vcf_record.get_format(true, false);

    std::string expected = "GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS";
    EXPECT_EQ(actual, expected);
}

TEST(VCFRecordTest, get_format___genotyping_from_coverage) {
    VCFRecord vcf_record;

    std::string actual = vcf_record.get_format(false, true);

    std::string expected = "GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF";
    EXPECT_EQ(actual, expected);
}


class VCFRecordTest___to_string___Fixture : public ::testing::Test {
private:
    class VCFRecordMock : public VCFRecord {
    public:
        using VCFRecord::VCFRecord;
        MOCK_METHOD(std::string, alts_to_string, (), (const));
        MOCK_METHOD(std::string, get_format, (bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage), (const));
        MOCK_METHOD(std::string, sample_infos_to_string, (bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage), (const));
    };

public:
    VCFRecordTest___to_string___Fixture() : vcf_record("chrom", 10, "ref", "", "info", "graph_type_info") {}

    void SetUp() override {
    }

    void TearDown() override {
    }

    VCFRecordMock vcf_record;
};

TEST_F(VCFRecordTest___to_string___Fixture, to_string) {
    EXPECT_CALL(vcf_record, alts_to_string)
    .Times(1)
    .WillOnce(Return("alts"));
    EXPECT_CALL(vcf_record, get_format)
    .Times(1)
    .WillOnce(Return("format"));
    EXPECT_CALL(vcf_record, sample_infos_to_string)
    .Times(1)
    .WillOnce(Return("sample_infos"));

    std::string actual = vcf_record.to_string(true, false);
    std::string expected("chrom\t11\t.\tref\talts\t.\t.\tinfo;graph_type_info\tformat\tsample_infos");

    EXPECT_EQ(actual, expected);
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






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PREVIOUS TESTS FROM VCF_RECORD FOLLOW
// COMMENTED OUT == NO NEED ANYMORE AND I PUT THE REASON WHY IT IS NOT NEEDED
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// format-related methods
// REASON WHY THESE ARE COMMENTED OUT: NO NEED ANYMORE, FORMATS ARE NOT VARIABLE ANYMORE, WHICH RENDERS THE CODE SIMPLER
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TEST(VCFRecordTest, add_formats_none)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vector<string> new_formats = {};
//    vr.add_formats(new_formats);
//    vector<string> expected_formats = { "GT" };
//    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
//}
//
//TEST(VCFRecordTest, add_formats_some)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vector<string> new_formats = { "hi", "there" };
//    vr.add_formats(new_formats);
//    vector<string> expected_formats = { "GT", "hi", "there" };
//    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
//}
//
//TEST(VCFRecordTest, add_formats_some_repeat)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vector<string> new_formats = { "hi", "there" };
//    vr.add_formats(new_formats);
//    vr.add_formats(new_formats);
//    vector<string> expected_formats = { "GT", "hi", "there" };
//    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
//}
//
//TEST(VCFRecordTest, add_formats_some_overlapping)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vector<string> new_formats = { "hi", "there" };
//    vr.add_formats(new_formats);
//    new_formats = { "hi", "again" };
//    vr.add_formats(new_formats);
//    vector<string> expected_formats = { "GT", "hi", "there", "again" };
//    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
//}
//
//TEST(VCFRecordTest, add_format_death_no_samples)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    uint16_t v = 20;
//    EXPECT_DEATH(vr.set_format(0, "hello", v), "");
//    float w = 20.0;
//    EXPECT_DEATH(vr.set_format(0, "hello", w), "");
//}
//
//TEST(VCFRecordTest, add_format_cap_too_big)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    uint32_t v = 60000000;
//    unordered_map<string, vector<uint16_t>> m;
//    vr.samples.push_back(m);
//    vr.set_format(0, "hello", v);
//    EXPECT_EQ(vr.get_format_u(0, "hello")[0], std::numeric_limits<uint16_t>::max() - 1);
//}
//
//TEST(VCFRecordTest, add_format_new_uint)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    unordered_map<string, vector<uint16_t>> m;
//    vr.samples.push_back(m);
//    uint16_t v = 20;
//    vr.set_format(0, "hello", v);
//    EXPECT_EQ(vr.samples.size(), 1);
//    EXPECT_TRUE(vr.samples[0].find("hello") != vr.samples[0].end());
//    std::vector<uint16_t> exp_v = { v };
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, vr.samples[0]["hello"], exp_v);
//    std::vector<std::string> exp_f = { "GT", "hello" };
//    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
//}
//
//TEST(VCFRecordTest, add_format_old_uint_overwritten)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    unordered_map<string, vector<uint16_t>> m;
//    m["hello"] = { 10 };
//    vr.samples.push_back(m);
//    uint16_t v = 20;
//    vr.set_format(0, "hello", v);
//    EXPECT_EQ(vr.samples.size(), 1);
//    EXPECT_TRUE(vr.samples[0].find("hello") != vr.samples[0].end());
//    std::vector<uint16_t> exp_v = { v };
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, vr.samples[0]["hello"], exp_v);
//    std::vector<std::string> exp_f = { "GT", "hello" };
//    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
//}
//
//TEST(VCFRecordTest, add_format_new_float)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    unordered_map<string, vector<uint16_t>> m;
//    vr.samples.push_back(m);
//    float v = 20.0;
//    vr.set_format(0, "hello", v);
//    EXPECT_EQ(vr.regt_samples.size(), 1);
//    EXPECT_TRUE(vr.regt_samples[0].find("hello") != vr.regt_samples[0].end());
//    std::vector<float> exp_v = { v };
//    EXPECT_ITERABLE_EQ(std::vector<float>, vr.regt_samples[0]["hello"], exp_v);
//    std::vector<std::string> exp_f = { "GT", "hello" };
//    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
//}
//
//TEST(VCFRecordTest, add_format_old_float_overwritten)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    unordered_map<string, vector<float>> m;
//    m["hello"] = {};
//    m["hello"].push_back(10.0);
//    vr.regt_samples.push_back(m);
//    float v = 20.0;
//    vr.set_format(0, "hello", v);
//    EXPECT_EQ(vr.regt_samples.size(), 1);
//    EXPECT_TRUE(vr.regt_samples[0].find("hello") != vr.regt_samples[0].end());
//    std::vector<float> exp_v = { v };
//    EXPECT_ITERABLE_EQ(std::vector<float>, vr.regt_samples[0]["hello"], exp_v);
//    std::vector<std::string> exp_f = { "GT", "hello" };
//    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
//}
//
//TEST(VCFRecordTest, append_format_old_uint)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    unordered_map<string, vector<uint16_t>> m;
//    vr.samples.push_back(m);
//    uint16_t v = 10;
//    vr.set_format(0, "hello", v);
//    v = 20;
//    vr.append_format(0, "hello", v);
//    EXPECT_EQ(vr.samples.size(), 1);
//    EXPECT_TRUE(vr.samples[0].find("hello") != vr.samples[0].end());
//    std::vector<uint16_t> exp_v = { 10, 20 };
//    EXPECT_ITERABLE_EQ(std::vector<uint16_t>, vr.samples[0]["hello"], exp_v);
//    std::vector<std::string> exp_f = { "GT", "hello" };
//    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
//}
//
//TEST(VCFRecordTest, append_format_old_float)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    unordered_map<string, vector<float>> m;
//    vr.regt_samples.push_back(m);
//    float v = 10.0;
//    vr.set_format(0, "hello", v);
//    v = 20.0;
//    vr.append_format(0, "hello", v);
//    EXPECT_EQ(vr.regt_samples.size(), 1);
//    EXPECT_TRUE(vr.regt_samples[0].find("hello") != vr.regt_samples[0].end());
//    std::vector<float> exp_v = { 10.0, 20.0 };
//    EXPECT_ITERABLE_EQ(std::vector<float>, vr.regt_samples[0]["hello"], exp_v);
//    std::vector<std::string> exp_f = { "GT", "hello" };
//    EXPECT_ITERABLE_EQ(std::vector<std::string>, vr.format, exp_f);
//}
//
//TEST(VCFRecordTest, get_format_float)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    unordered_map<string, vector<float>> m;
//    m.reserve(3);
//    vr.regt_samples.push_back(m);
//    float v = 10.0;
//    vr.set_format(0, "hello", v);
//
//    auto res = vr.get_format_f(1, "hello");
//    EXPECT_EQ(res.size(), 0);
//    EXPECT_TRUE(res.empty());
//
//    res = vr.get_format_f(0, "help");
//    EXPECT_EQ(res.size(), 0);
//    EXPECT_TRUE(res.empty());
//
//    res = vr.get_format_f(0, "hello");
//    EXPECT_EQ(res.size(), 1);
//    EXPECT_FALSE(res.empty());
//}
//
//TEST(VCFRecordTest, get_format_uint)
//{
//    VCFRecord vr("chrom1", 3, "A", "T");
//    unordered_map<string, vector<uint16_t>> m;
//    m.reserve(3);
//    vr.samples.push_back(m);
//    uint16_t v = 10;
//    vr.set_format(0, "hello", v);
//
//    auto res = vr.get_format_u(1, "hello");
//    EXPECT_EQ(res.size(), 0);
//    EXPECT_TRUE(res.empty());
//
//    res = vr.get_format_u(0, "help");
//    EXPECT_EQ(res.size(), 0);
//    EXPECT_TRUE(res.empty());
//
//    res = vr.get_format_u(0, "hello");
//    EXPECT_EQ(res.size(), 1);
//    EXPECT_FALSE(res.empty());
//}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ostream methods
//ALL ostream TESTS WERE REPLACED BY to_string TESTS - MUCH OF THE COMPLEXITY HERE WAS MOVED TO OTHER METHODS/CONCEPTS
//e.g. VCFRecord::alts_to_string(), VCFRecord::get_format(), SampleIndexToSampleInfo::to_string(), SampleInfo::to_string(), etc...
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TEST(VCFRecordTest, ostream) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vector<string> v = {"chrom1", "4", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT"};
//    stringstream out;
//    out << vr;
//    string rr;
//    for (const auto &s : v) {
//        out >> rr;
//        EXPECT_EQ(s, rr);
//    }
//}
//
//TEST(VCFRecordTest, ostream_with_sample_not_all_info_in_formats) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    SampleInfo<uint16_t> sample_info;
//    sample_info["GT"] = {1};
//    sample_info["pringle"] = {2};
//    vr.sampleIndex_to_sampleInfo.push_back(sample_info);
//    vector<string> v = {"chrom1", "4", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT"};
//    stringstream out;
//    out << vr;
//    string rr;
//    for (const auto &s : v) {
//        out >> rr;
//        EXPECT_EQ(s, rr);
//    }
//    uint16_t u = 1;
//    uint ru;
//    out >> ru;
//    EXPECT_EQ(u, ru);
//}
//
//TEST(VCFRecordTest, ostream_with_sample_including_all_formats) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    SampleInfo<uint16_t> sample_info;
//    sample_info["GT"] = {0};
//    sample_info["pringle"] = {2};
//    vr.sampleIndex_to_sampleInfo.push_back(sample_info);
//    vr.add_formats({"pringle"});
//    vector<string> v = {"chrom1", "4", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT:pringle"};
//    vector<uint16_t> vu = {0, 2};
//    stringstream out;
//    out << vr;
//    string rr;
//    for (const auto &s : v) {
//        out >> rr;
//        EXPECT_EQ(s, rr);
//    }
//    uint ru;
//    for (const auto &s : vu) {
//        out >> ru;
//        EXPECT_EQ(s, ru);
//        out.ignore(1, ':');
//    }
//
//}
//
//TEST(VCFRecordTest, ostream_with_sample_more_formats_than_info) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    SampleInfo<uint16_t> sample_info;
//    sample_info["GT"] = {0};
//    vr.sampleIndex_to_sampleInfo.push_back(sample_info);
//    vr.add_formats({"pringle"});
//    vector<string> v = {"chrom1", "4", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT:pringle"};
//    vector<uint16_t> vu = {0};
//    stringstream out;
//    out << vr;
//    string rr;
//    for (const auto &s : v) {
//        out >> rr;
//        EXPECT_EQ(s, rr);
//    }
//    uint ru;
//    uint16_t u = 0;
//    out >> ru;
//    EXPECT_EQ(u, ru);
//    out >> rr;
//    EXPECT_EQ(":.", rr);
//}
//
//TEST(VCFRecordTest, ostream_with_sample_more_formats_than_info_regt) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    SampleInfo<uint16_t> sample_info_int;
//    sample_info_int["GT"] = {0};
//    vr.sampleIndex_to_sampleInfo.push_back(sample_info_int);
//    SampleInfo<float> sample_info_float;
//    sample_info_float["pringle"] = {0.1};
//    vr.sampleIndex_to_sampleInfo.push_back(sample_info_float);
//    vr.add_formats({"pringle"});
//    vector<string> v = {"chrom1", "4", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT:pringle"};
//    stringstream out;
//    out << vr;
//    string rr;
//    for (const auto &s : v) {
//        out >> rr;
//        EXPECT_EQ(s, rr);
//    }
//    uint ru;
//    uint16_t u = 0;
//    out >> ru;
//    EXPECT_EQ(u, ru);
//    out.ignore(1, ':');
//    float rf = 0.0, f = 0.1;
//    out >> rf;
//    EXPECT_EQ(f, rf);
//}
//
//TEST(VCFRecordTest, ostream_with_zero_pos) {
//    VCFRecord vr("chrom1", 0, "A", "T");
//    SampleInfo<uint16_t> sample_info_int;
//    sample_info_int["GT"] = {0};
//    vr.sampleIndex_to_sampleInfo.push_back(sample_info_int);
//    SampleInfo<float> sample_info_float;
//    sample_info_float["pringle"] = {0.1};
//    vr.sampleIndex_to_sampleInfo.push_back(sample_info_float);
//    vr.add_formats({"pringle"});
//    vector<string> v = {"chrom1", "1", ".", "A", "T", ".", ".", "SVTYPE=SNP", "GT:pringle"};
//    stringstream out;
//    out << vr;
//    string rr;
//    for (const auto &s : v) {
//        out >> rr;
//        EXPECT_EQ(s, rr);
//    }
//    uint ru;
//    uint16_t u = 0;
//    out >> ru;
//    EXPECT_EQ(u, ru);
//    out.ignore(1, ':');
//    float rf = 0.0, f = 0.1;
//    out >> rf;
//    EXPECT_EQ(f, rf);
//}