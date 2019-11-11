#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "test_macro.cpp"
#include "vcfrecord.h"
#include <stdint.h>
#include <iostream>
#include <cmath>
#include "test_helpers.h"
#include "vcf.h"


using namespace std;
using ::testing::Return;

TEST(VCFRecordTest, create_empty) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf);
    EXPECT_EQ(&vcf, vr.parent_vcf);
    EXPECT_EQ(".", vr.chrom);
    EXPECT_EQ((uint) 0, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ(".", vr.get_ref());
    EXPECT_EQ((uint) 0, vr.get_alts().size());
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);
    EXPECT_EQ(1, vr.sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint) 1, vr.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST(VCFRecordTest, create_with_values) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf,"chrom1", 3, "A", "T");
    EXPECT_EQ(&vcf, vr.parent_vcf);
    EXPECT_EQ("chrom1", vr.chrom);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.get_ref());
    EXPECT_EQ(1, vr.get_alts().size());
    EXPECT_EQ("T", vr.get_alts()[0]);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ("SVTYPE=SNP", vr.info);
    EXPECT_EQ(1, vr.sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint) 2, vr.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST(VCFRecordTest, create_from_record) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord template_vr(&vcf,"chrom1", 3, "A", "T");
    VCFRecord vr(template_vr);
    EXPECT_EQ(&vcf, vr.parent_vcf);
    EXPECT_EQ("chrom1", vr.chrom);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.get_ref());
    EXPECT_EQ(1, vr.get_alts().size());
    EXPECT_EQ("T", vr.get_alts()[0]);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ("SVTYPE=SNP", vr.info);
    EXPECT_EQ(1, vr.sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint) 2, vr.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST(VCFRecordTest, create_from_record_with_samples) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord template_vr(&vcf,"chrom1", 3, "A", "T");
    VCFRecord vr(template_vr);
    EXPECT_EQ(&vcf, vr.parent_vcf);
    EXPECT_EQ("chrom1", vr.chrom);
    EXPECT_EQ((uint) 3, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.get_ref());
    EXPECT_EQ(1, vr.get_alts().size());
    EXPECT_EQ("T", vr.get_alts()[0]);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ("SVTYPE=SNP", vr.info);
    EXPECT_EQ(1, vr.sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint) 2, vr.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}



class VCFRecordTest___default_VCF_Record___Fixture : public ::testing::Test {
public:
    VCFRecordTest___default_VCF_Record___Fixture() : vcf(create_VCF_with_default_parameters()), vcf_record(&vcf) {}
    VCF vcf;
    VCFRecord vcf_record;
};


TEST_F(VCFRecordTest___default_VCF_Record___Fixture, set_ref___ref_is_empty) {
    vcf_record.set_ref("");

    EXPECT_EQ(".", vcf_record.get_ref());
    EXPECT_EQ(1, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, set_ref___ref_is_dot) {
    vcf_record.set_ref(".");

    EXPECT_EQ(".", vcf_record.get_ref());
    EXPECT_EQ(1, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, set_ref___ref_is_valid) {
    vcf_record.set_ref("ACGT");

    EXPECT_EQ("ACGT", vcf_record.get_ref());
    EXPECT_EQ(1, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, add_new_alt___alt_is_empty___not_added) {
    vcf_record.add_new_alt("");

    EXPECT_EQ(0, vcf_record.get_alts().size());
    EXPECT_EQ(1, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, add_new_alt___alt_is_dot___not_added) {
    vcf_record.add_new_alt(".");

    EXPECT_EQ(0, vcf_record.get_alts().size());
    EXPECT_EQ(1, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, add_new_alt___alt_is_valid) {
    vcf_record.add_new_alt("AC");

    EXPECT_EQ(1, vcf_record.get_alts().size());
    EXPECT_EQ("AC", vcf_record.get_alts()[0]);
    EXPECT_EQ(2, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, add_new_alt___add_two_valid_alts) {
    vcf_record.add_new_alt("AC");
    vcf_record.add_new_alt("AG");

    EXPECT_EQ(2, vcf_record.get_alts().size());
    EXPECT_EQ("AC", vcf_record.get_alts()[0]);
    EXPECT_EQ("AG", vcf_record.get_alts()[1]);
    EXPECT_EQ(3, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, add_new_alt___add_two_valid_alts_and_several_repeated_alts___repeated_do_not_get_added) {
    vcf_record.add_new_alt("AC");
    vcf_record.add_new_alt("AG");
    vcf_record.add_new_alt("AC");
    vcf_record.add_new_alt("AG");

    EXPECT_EQ(2, vcf_record.get_alts().size());
    EXPECT_EQ("AC", vcf_record.get_alts()[0]);
    EXPECT_EQ("AG", vcf_record.get_alts()[1]);
    EXPECT_EQ(3, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, add_new_alts___add_two_valid_alts_and_several_repeated_alts___repeated_do_not_get_added) {
    std::vector<std::string> alts{"AC", "AG", ".", "", "AC", "AG"};

    vcf_record.add_new_alts(alts.begin(), alts.end());

    EXPECT_EQ(2, vcf_record.get_alts().size());
    EXPECT_EQ("AC", vcf_record.get_alts()[0]);
    EXPECT_EQ("AG", vcf_record.get_alts()[1]);
    EXPECT_EQ(3, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST(VCFRecordTest, clear_simple) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf,"chrom1", 3, "A", "T");
    vr.clear();
    EXPECT_EQ(&vcf, vr.parent_vcf);
    EXPECT_EQ(".", vr.chrom);
    EXPECT_EQ((uint) 0, vr.pos);
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ(".", vr.get_ref());
    EXPECT_EQ((uint) 0, vr.get_alts().size());
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);
    EXPECT_EQ((uint) 1, vr.sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint) 1, vr.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST(VCFRecordTest, equals) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf);
    EXPECT_EQ(vr, vr);

    VCFRecord vr1(&vcf, "chrom1", 3, "A", "T");
    EXPECT_EQ(vr1, vr1);
    EXPECT_EQ((vr == vr1), false);
    EXPECT_EQ((vr1 == vr), false);

    VCFRecord vr2(&vcf, "chrom2", 3, "A", "T");
    EXPECT_EQ(vr2, vr2);
    EXPECT_EQ((vr2 == vr1), false);
    EXPECT_EQ((vr1 == vr2), false);

    VCFRecord vr3(&vcf, "chrom1", 6, "A", "T");
    EXPECT_EQ(vr3, vr3);
    EXPECT_EQ((vr3 == vr1), false);
    EXPECT_EQ((vr1 == vr3), false);

    VCFRecord vr4(&vcf, "chrom1", 3, "G", "T");
    EXPECT_EQ(vr4, vr4);
    EXPECT_EQ((vr4 == vr1), false);
    EXPECT_EQ((vr1 == vr4), false);

    VCFRecord vr5(&vcf, "chrom1", 3, "A", "G");
    EXPECT_EQ(vr5, vr5);
    EXPECT_EQ((vr5 == vr1), false);
    EXPECT_EQ((vr1 == vr5), false);

}

TEST(VCFRecordTest, less_than) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr1(&vcf, "chrom1", 3, "A", "T");
    VCFRecord vr2(&vcf, "chrom2", 3, "A", "T");
    EXPECT_EQ((vr2 < vr1), false);
    EXPECT_EQ((vr1 < vr2), true);

    VCFRecord vr3(&vcf, "chrom1", 6, "A", "T");
    EXPECT_EQ((vr3 < vr1), false);
    EXPECT_EQ((vr1 < vr3), true);

    VCFRecord vr4(&vcf, "chrom1", 3, "G", "T");
    EXPECT_EQ((vr4 < vr1), false);
    EXPECT_EQ((vr1 < vr4), true);

    VCFRecord vr5(&vcf, "chrom1", 3, "A", "G");
    EXPECT_EQ((vr5 < vr1), true);
    EXPECT_EQ((vr1 < vr5), false);
}


TEST(VCFRecordTest, alts_to_string___no_alts) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);

    std::string actual = vcf_record.alts_to_string();

    EXPECT_EQ(".", actual);
}

TEST(VCFRecordTest, alts_to_string___three_alts) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);
    vcf_record.add_new_alt("A1");
    vcf_record.add_new_alt("A2");
    vcf_record.add_new_alt("A3");

    std::string actual = vcf_record.alts_to_string();

    EXPECT_EQ("A1,A2,A3", actual);
}


TEST(VCFRecordTest, get_format___no_flags_set___expects_death) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);
    EXPECT_DEATH(vcf_record.get_format(false, false), "");
}

TEST(VCFRecordTest, get_format___both_flags_set___expects_death) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);
    EXPECT_DEATH(vcf_record.get_format(true, true), "");
}

TEST(VCFRecordTest, get_format___genotyping_from_maximum_likelihood) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);

    std::string actual = vcf_record.get_format(true, false);

    std::string expected = "GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS";
    EXPECT_EQ(actual, expected);
}

TEST(VCFRecordTest, get_format___genotyping_from_coverage) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);

    std::string actual = vcf_record.get_format(false, true);

    std::string expected = "GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF";
    EXPECT_EQ(actual, expected);
}


class VCFRecordTest___to_string___Fixture : public ::testing::Test {
private:
    class VCFRecordMock : public VCFRecord {
    public:
        using VCFRecord::VCFRecord;
        MOCK_METHOD(std::string, alts_to_string, (), (const override));
        MOCK_METHOD(std::string, get_format, (bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage), (const override));
        MOCK_METHOD(std::string, sample_infos_to_string, (bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage), (const override));
    };

public:
    VCFRecordTest___to_string___Fixture() : vcf(create_VCF_with_default_parameters()), vcf_record(&vcf, "chrom", 10, "ref", "", "info", "graph_type_info") {}

    void SetUp() override {
    }

    void TearDown() override {
    }
    VCF vcf;
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
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf, "chrom1", 1, "ACGT", "A", ".", ".");

    size_t actual = vr.get_longest_allele_length();
    size_t expected = 4;

    EXPECT_EQ(actual, expected);
}

TEST(VCFRecordTest, get_longest_allele_length___longest_allele_is_alt) {
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf, "chrom1", 1, "ACGT", "A", ".", ".");
    vr.add_new_alt("ACGTTTTAC");
    vr.add_new_alt("C");

    size_t actual = vr.get_longest_allele_length();
    size_t expected = 9;

    EXPECT_EQ(actual, expected);
}

class VCFRecordTest___ref_allele_is_inside_given_interval______Fixture : public ::testing::Test {
public:
    VCFRecordTest___ref_allele_is_inside_given_interval______Fixture() :
        vcf(create_VCF_with_default_parameters()),
        vcf_record(&vcf, "chrom", 10, "REF", ""){}
    VCF vcf;
    VCFRecord vcf_record;

};
TEST_F(VCFRecordTest___ref_allele_is_inside_given_interval______Fixture, different_chroms) {
    EXPECT_FALSE(vcf_record.ref_allele_is_inside_given_interval("dif_chrom", 0, 0));
}
TEST_F(VCFRecordTest___ref_allele_is_inside_given_interval______Fixture, ends_1_position_short) {
    EXPECT_FALSE(vcf_record.ref_allele_is_inside_given_interval("chrom", 10, 12));
}
TEST_F(VCFRecordTest___ref_allele_is_inside_given_interval______Fixture, starts_1_position_after) {
    EXPECT_FALSE(vcf_record.ref_allele_is_inside_given_interval("chrom", 11, 13));
}
TEST_F(VCFRecordTest___ref_allele_is_inside_given_interval______Fixture, exact_interval) {
    EXPECT_TRUE(vcf_record.ref_allele_is_inside_given_interval("chrom", 10, 13));
}
TEST_F(VCFRecordTest___ref_allele_is_inside_given_interval______Fixture, larger_interval) {
    EXPECT_TRUE(vcf_record.ref_allele_is_inside_given_interval("chrom", 9, 14));
}



class VCFRecordTest___merge_record_into_this______Fixture : public ::testing::Test {
public:
    VCFRecordTest___merge_record_into_this______Fixture() :
            vcf(create_VCF_with_default_parameters()),
            vcf_record_only_ref_no_alts(&vcf, "1", 1, "A", ""),
            vcf_record_ref_A_alt_T(&vcf, "1", 1, "A", "T"),
            vcf_record_ref_A_alt_TTT(&vcf, "1", 1, "A", "TTT"),
            vcf_record_ref_A_alt_T_TT_TTT(&vcf, "1", 1, "A", "T")
    {
        vcf_record_ref_A_alt_T_TT_TTT.add_new_alt("TT");
        vcf_record_ref_A_alt_T_TT_TTT.add_new_alt("TTT");
    }
    VCF vcf;
    VCFRecord vcf_record_only_ref_no_alts;
    VCFRecord vcf_record_ref_A_alt_T;
    VCFRecord vcf_record_ref_A_alt_TTT;
    VCFRecord vcf_record_ref_A_alt_T_TT_TTT;
};

// TODO : mock this->sampleIndex_to_sampleInfo.merge_other_samples_infos_into_this()
TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_T_into_no_alts) {
    vcf_record_only_ref_no_alts.merge_record_into_this(vcf_record_ref_A_alt_T);

    EXPECT_EQ(1, vcf_record_only_ref_no_alts.get_alts().size());
    EXPECT_EQ("T", vcf_record_only_ref_no_alts.get_alts()[0]);
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_no_alts_into_T) {
    vcf_record_ref_A_alt_T.merge_record_into_this(vcf_record_only_ref_no_alts);

    EXPECT_EQ(1, vcf_record_ref_A_alt_T.get_alts().size());
    EXPECT_EQ("T", vcf_record_ref_A_alt_T.get_alts()[0]);
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_TTT_into_T) {
    vcf_record_ref_A_alt_T.merge_record_into_this(vcf_record_ref_A_alt_TTT);

    EXPECT_EQ(2, vcf_record_ref_A_alt_T.get_alts().size());
    EXPECT_EQ("T", vcf_record_ref_A_alt_T.get_alts()[0]);
    EXPECT_EQ("TTT", vcf_record_ref_A_alt_T.get_alts()[1]);
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_first_alt_is_common___expects_death) {
    EXPECT_DEATH(vcf_record_ref_A_alt_T_TT_TTT.merge_record_into_this(vcf_record_ref_A_alt_T), "");
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_last_alt_is_common___expects_death) {
    EXPECT_DEATH(vcf_record_ref_A_alt_T_TT_TTT.merge_record_into_this(vcf_record_ref_A_alt_TTT), "");
}

class VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture : public ::testing::Test {
public:
    VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture() :
        vcf(create_VCF_with_default_parameters()),
        vcf_record_only_ref_no_alts(&vcf, "1", 1, "A", ""),
        vcf_record_tri_allelic(&vcf, "1", 1, "A", "T"),
        vcf_record_null_reference(&vcf, "1", 1, ".", "T"),
        vcf_record_ref_A(&vcf, "1", 1, "A", "T"),
        vcf_record_ref_C(&vcf, "1", 1, "C", "CT"),
        vcf_record_ref_A_too_long_alt(&vcf, "1", 1, "A", "TTTTTT"),
        vcf_record_ref_A_same_chrom_different_pos(&vcf, "1", 2, "A", "T"),
        vcf_record_ref_A_different_chrom_same_pos(&vcf, "2", 1, "A", "T"),
        vcf_record_ref_A_different_alt(&vcf, "1", 1, "A", "C"),
        vcf_record_ref_A_same_alt(&vcf, "1", 1, "A", "T")
    {
        vcf_record_tri_allelic.add_new_alt("TT");
    }
    VCF vcf;
    VCFRecord vcf_record_only_ref_no_alts;
    VCFRecord vcf_record_tri_allelic;
    VCFRecord vcf_record_null_reference;
    VCFRecord vcf_record_ref_A;
    VCFRecord vcf_record_ref_C;
    VCFRecord vcf_record_ref_A_too_long_alt;
    VCFRecord vcf_record_ref_A_same_chrom_different_pos;
    VCFRecord vcf_record_ref_A_different_chrom_same_pos;
    VCFRecord vcf_record_ref_A_different_alt;
    VCFRecord vcf_record_ref_A_same_alt;
};

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, merge_only_ref_no_alts___expects_death) {
    EXPECT_DEATH(vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_only_ref_no_alts), "");
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, merge_triallelic___expects_death) {
    EXPECT_DEATH(vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_tri_allelic), "");
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, this_has_null_ref) {
    bool actual = vcf_record_null_reference.can_biallelic_record_be_merged_into_this(vcf_record_ref_A);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, other_has_null_ref) {
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_null_reference);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, different_refs) {
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_ref_C);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, this_has_alt_too_long) {
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_ref_A_too_long_alt, 5);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, other_has_alt_too_long) {
    bool actual = vcf_record_ref_A_too_long_alt.can_biallelic_record_be_merged_into_this(vcf_record_ref_A, 5);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, same_record) {
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_ref_A);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, same_chrom_different_pos) {
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_ref_A_same_chrom_different_pos);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, same_different_chrom_same_pos) {
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_ref_A_different_chrom_same_pos);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, can_be_merged) {
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_ref_A_different_alt);
    EXPECT_TRUE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, same_alt_different_record) {
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_ref_A_same_alt);
    EXPECT_FALSE(actual);
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
//    vector<string> new_formats = {};
//    vr.add_formats(new_formats);
//    vector<string> expected_formats = { "GT" };
//    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
//}
//
//TEST(VCFRecordTest, add_formats_some)
//{
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
//    vector<string> new_formats = { "hi", "there" };
//    vr.add_formats(new_formats);
//    vector<string> expected_formats = { "GT", "hi", "there" };
//    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
//}
//
//TEST(VCFRecordTest, add_formats_some_repeat)
//{
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
//    vector<string> new_formats = { "hi", "there" };
//    vr.add_formats(new_formats);
//    vr.add_formats(new_formats);
//    vector<string> expected_formats = { "GT", "hi", "there" };
//    EXPECT_ITERABLE_EQ(vector<string>, expected_formats, vr.format);
//}
//
//TEST(VCFRecordTest, add_formats_some_overlapping)
//{
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
//    uint16_t v = 20;
//    EXPECT_DEATH(vr.set_format(0, "hello", v), "");
//    float w = 20.0;
//    EXPECT_DEATH(vr.set_format(0, "hello", w), "");
//}
//
//TEST(VCFRecordTest, add_format_cap_too_big)
//{
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
//    uint32_t v = 60000000;
//    unordered_map<string, vector<uint16_t>> m;
//    vr.samples.push_back(m);
//    vr.set_format(0, "hello", v);
//    EXPECT_EQ(vr.get_format_u(0, "hello")[0], std::numeric_limits<uint16_t>::max() - 1);
//}
//
//TEST(VCFRecordTest, add_format_new_uint)
//{
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
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
//    VCFRecord vr(&vcf, "chrom1", 0, "A", "T");
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