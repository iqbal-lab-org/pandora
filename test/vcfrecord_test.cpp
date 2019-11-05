#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "vcfrecord.h"
#include <stdint.h>
#include <iostream>
#include <cmath>
#include "test_helpers.h"


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

//
//
//TEST(VCFRecordTest, equals) {
//    VCFRecord vr;
//    EXPECT_EQ(vr, vr);
//
//    VCFRecord vr1("chrom1", 3, "A", "T");
//    EXPECT_EQ(vr1, vr1);
//    EXPECT_EQ((vr == vr1), false);
//    EXPECT_EQ((vr1 == vr), false);
//
//    VCFRecord vr2("chrom2", 3, "A", "T");
//    EXPECT_EQ(vr2, vr2);
//    EXPECT_EQ((vr2 == vr1), false);
//    EXPECT_EQ((vr1 == vr2), false);
//
//    VCFRecord vr3("chrom1", 6, "A", "T");
//    EXPECT_EQ(vr3, vr3);
//    EXPECT_EQ((vr3 == vr1), false);
//    EXPECT_EQ((vr1 == vr3), false);
//
//    VCFRecord vr4("chrom1", 3, "G", "T");
//    EXPECT_EQ(vr4, vr4);
//    EXPECT_EQ((vr4 == vr1), false);
//    EXPECT_EQ((vr1 == vr4), false);
//
//    VCFRecord vr5("chrom1", 3, "A", "G");
//    EXPECT_EQ(vr5, vr5);
//    EXPECT_EQ((vr5 == vr1), false);
//    EXPECT_EQ((vr1 == vr5), false);
//
//}
//
//TEST(VCFRecordTest, less_than) {
//    VCFRecord vr1("chrom1", 3, "A", "T");
//    VCFRecord vr2("chrom2", 3, "A", "T");
//    EXPECT_EQ((vr2 < vr1), false);
//    EXPECT_EQ((vr1 < vr2), true);
//
//    VCFRecord vr3("chrom1", 6, "A", "T");
//    EXPECT_EQ((vr3 < vr1), false);
//    EXPECT_EQ((vr1 < vr3), true);
//
//    VCFRecord vr4("chrom1", 3, "G", "T");
//    EXPECT_EQ((vr4 < vr1), false);
//    EXPECT_EQ((vr1 < vr4), true);
//
//    VCFRecord vr5("chrom1", 3, "A", "G");
//    EXPECT_EQ((vr5 < vr1), true);
//    EXPECT_EQ((vr1 < vr5), false);
//}
//
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
//
//
//TEST(VCFRecordTest, get_longest_allele_length___longest_allele_is_ref) {
//    VCFRecord vr("chrom1", 1, "ACGT", "A", ".", ".");
//
//    size_t actual = vr.get_longest_allele_length();
//    size_t expected = 4;
//
//    EXPECT_EQ(actual, expected);
//}
//
//TEST(VCFRecordTest, get_longest_allele_length___longest_allele_is_alt) {
//    VCFRecord vr("chrom1", 1, "ACGT", "A", ".", ".");
//    vr.alts.push_back("ACGTTTTAC");
//    vr.alts.push_back("C");
//
//    size_t actual = vr.get_longest_allele_length();
//    size_t expected = 9;
//
//    EXPECT_EQ(actual, expected);
//}