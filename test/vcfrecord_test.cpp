#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "test_macro.cpp"
#include "vcfrecord.h"
#include <stdint.h>
#include <iostream>
#include <cmath>
#include "test_helpers_containers.h"
#include "test_helpers.h"
#include "vcf.h"

using namespace std;
using ::testing::Return;

TEST(VCFRecordTest, create_empty)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf);
    EXPECT_EQ(&vcf, vr.parent_vcf);
    EXPECT_EQ(".", vr.get_chrom());
    EXPECT_EQ((uint)0, vr.get_pos());
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ(".", vr.get_ref());
    EXPECT_EQ((uint)0, vr.get_alts().size());
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);
    EXPECT_EQ(1, vr.sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)1, vr.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST(VCFRecordTest, create_with_values)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
    EXPECT_EQ(&vcf, vr.parent_vcf);
    EXPECT_EQ("chrom1", vr.get_chrom());
    EXPECT_EQ((uint)3, vr.get_pos());
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.get_ref());
    EXPECT_EQ(1, vr.get_alts().size());
    EXPECT_EQ("T", vr.get_alts()[0]);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(VCF::VARIANT_CLASS_ID + "=SNP", vr.info);
    EXPECT_EQ(1, vr.sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vr.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST(VCFRecordTest, create_from_record)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord template_vr(&vcf, "chrom1", 3, "A", "T");
    VCFRecord vr(template_vr);
    EXPECT_EQ(&vcf, vr.parent_vcf);
    EXPECT_EQ("chrom1", vr.get_chrom());
    EXPECT_EQ((uint)3, vr.get_pos());
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.get_ref());
    EXPECT_EQ(1, vr.get_alts().size());
    EXPECT_EQ("T", vr.get_alts()[0]);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(VCF::VARIANT_CLASS_ID + "=SNP", vr.info);
    EXPECT_EQ(1, vr.sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vr.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST(VCFRecordTest, create_from_record_with_samples)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord template_vr(&vcf, "chrom1", 3, "A", "T");
    VCFRecord vr(template_vr);
    EXPECT_EQ(&vcf, vr.parent_vcf);
    EXPECT_EQ("chrom1", vr.get_chrom());
    EXPECT_EQ((uint)3, vr.get_pos());
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ("A", vr.get_ref());
    EXPECT_EQ(1, vr.get_alts().size());
    EXPECT_EQ("T", vr.get_alts()[0]);
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(VCF::VARIANT_CLASS_ID + "=SNP", vr.info);
    EXPECT_EQ(1, vr.sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)2, vr.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

class VCFRecordTest___default_VCF_Record___Fixture : public ::testing::Test {
public:
    VCFRecordTest___default_VCF_Record___Fixture()
        : vcf(create_VCF_with_default_parameters())
        , vcf_record(&vcf)
    {
    }
    VCF vcf;
    VCFRecord vcf_record;
};

TEST_F(
    VCFRecordTest___default_VCF_Record___Fixture, set_ref_and_clear_alts___ref_is_empty)
{
    vcf_record.set_ref_and_clear_alts("");

    EXPECT_EQ(".", vcf_record.get_ref());
    EXPECT_EQ(1, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(
    VCFRecordTest___default_VCF_Record___Fixture, set_ref_and_clear_alts___ref_is_dot)
{
    vcf_record.set_ref_and_clear_alts(".");

    EXPECT_EQ(".", vcf_record.get_ref());
    EXPECT_EQ(1, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(
    VCFRecordTest___default_VCF_Record___Fixture, set_ref_and_clear_alts___ref_is_valid)
{
    vcf_record.set_ref_and_clear_alts("ACGT");

    EXPECT_EQ("ACGT", vcf_record.get_ref());
    EXPECT_EQ(1, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture,
    set_ref_and_clear_alts___check_if_alts_are_cleared___and_sample_infos_are_resized)
{
    vcf_record.add_new_alt("G");
    vcf_record.add_new_alt("C");
    vcf_record.add_new_alt("T");
    EXPECT_EQ(4, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());

    vcf_record.set_ref_and_clear_alts("ACGT");

    EXPECT_EQ("ACGT", vcf_record.get_ref());
    EXPECT_EQ(0, vcf_record.get_alts().size());
    EXPECT_EQ(1, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture,
    add_new_alt___alt_is_empty___added_as_dot_allele)
{
    vcf_record.add_new_alt("");

    EXPECT_EQ(1, vcf_record.get_alts().size());
    EXPECT_EQ(".", vcf_record.get_alts()[0]);
    EXPECT_EQ(2, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, add_new_alt___alt_is_dot___added)
{
    vcf_record.add_new_alt(".");

    EXPECT_EQ(1, vcf_record.get_alts().size());
    EXPECT_EQ(".", vcf_record.get_alts()[0]);
    EXPECT_EQ(2, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture,
    add_new_alt___add_a_dot_alt___then_a_valid_alt___both_are_added)
{
    vcf_record.add_new_alt(".");
    vcf_record.add_new_alt("A");

    EXPECT_EQ(2, vcf_record.get_alts().size());
    EXPECT_EQ(".", vcf_record.get_alts()[0]);
    EXPECT_EQ("A", vcf_record.get_alts()[1]);
    EXPECT_EQ(3, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, add_new_alt___alt_is_valid)
{
    vcf_record.add_new_alt("AC");

    EXPECT_EQ(1, vcf_record.get_alts().size());
    EXPECT_EQ("AC", vcf_record.get_alts()[0]);
    EXPECT_EQ(2, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture, add_new_alt___add_two_valid_alts)
{
    vcf_record.add_new_alt("AC");
    vcf_record.add_new_alt("AG");

    EXPECT_EQ(2, vcf_record.get_alts().size());
    EXPECT_EQ("AC", vcf_record.get_alts()[0]);
    EXPECT_EQ("AG", vcf_record.get_alts()[1]);
    EXPECT_EQ(3, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture,
    add_new_alt___add_two_valid_alts_and_several_repeated_alts___expects_FatalRuntimeError)
{
    vcf_record.add_new_alt("AC");
    vcf_record.add_new_alt("AG");
    ASSERT_EXCEPTION(vcf_record.add_new_alt("AC"), FatalRuntimeError,
        "Error adding new ALT to a VCF record: ALT already exists");
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture,
    add_new_alts___add_two_valid_alts_and_a_empty_alt___all_are_added)
{
    std::vector<std::string> alts { "AC", "AG", "" };
    vcf_record.add_new_alts(alts.begin(), alts.end());

    EXPECT_EQ(3, vcf_record.get_alts().size());
    EXPECT_EQ("AC", vcf_record.get_alts()[0]);
    EXPECT_EQ("AG", vcf_record.get_alts()[1]);
    EXPECT_EQ(".", vcf_record.get_alts()[2]);
    EXPECT_EQ(4, vcf_record.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___default_VCF_Record___Fixture,
    add_new_alts___add_two_valid_alts_and_a_repeated_alt___expects_FatalRuntimeError)
{
    std::vector<std::string> alts { "AC", "AG", "",
        "." }; // NB: "" and "." are repeated because "" is translated to "."
    ASSERT_EXCEPTION(vcf_record.add_new_alts(alts.begin(), alts.end()),
        FatalRuntimeError, "Error adding new ALT to a VCF record: ALT already exists");
}

TEST(VCFRecordTest, clear_simple)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf, "chrom1", 3, "A", "T");
    vr.clear();
    EXPECT_EQ(&vcf, vr.parent_vcf);
    EXPECT_EQ(".", vr.get_chrom());
    EXPECT_EQ((uint)0, vr.get_pos());
    EXPECT_EQ(".", vr.id);
    EXPECT_EQ(".", vr.get_ref());
    EXPECT_EQ((uint)0, vr.get_alts().size());
    EXPECT_EQ(".", vr.qual);
    EXPECT_EQ(".", vr.filter);
    EXPECT_EQ(".", vr.info);
    EXPECT_EQ((uint)1, vr.sampleIndex_to_sampleInfo.size());
    EXPECT_EQ((uint)1, vr.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST(VCFRecordTest, equals)
{
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

TEST(VCFRecordTest, less_than)
{
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

TEST(VCFRecordTest, alts_to_string___no_alts)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);

    std::string actual = vcf_record.alts_to_string();

    EXPECT_EQ(".", actual);
}

TEST(VCFRecordTest, alts_to_string___three_alts)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);
    vcf_record.add_new_alt("A1");
    vcf_record.add_new_alt("A2");
    vcf_record.add_new_alt("A3");

    std::string actual = vcf_record.alts_to_string();

    EXPECT_EQ("A1,A2,A3", actual);
}

TEST(VCFRecordTest, get_format___no_flags_set___expects_FatalRuntimeError)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);
    ASSERT_EXCEPTION(vcf_record.get_format(false, false), FatalRuntimeError,
        "Error on getting format field from VCF record: incompatible genotyping "
        "options");
}

TEST(VCFRecordTest, get_format___both_flags_set___expects_FatalRuntimeError)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);
    ASSERT_EXCEPTION(vcf_record.get_format(true, true), FatalRuntimeError,
        "Error on getting format field from VCF record: incompatible genotyping "
        "options");
}

TEST(VCFRecordTest, get_format___genotyping_from_maximum_likelihood)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);

    std::string actual = vcf_record.get_format(true, false);

    std::string expected = "GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:"
                           "SUM_FWD_COVG:SUM_REV_COVG:GAPS";
    EXPECT_EQ(actual, expected);
}

TEST(VCFRecordTest, get_format___genotyping_from_coverage)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vcf_record(&vcf);

    std::string actual = vcf_record.get_format(false, true);

    std::string expected = "GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:"
                           "SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF";
    EXPECT_EQ(actual, expected);
}

class VCFRecordTest___to_string___Fixture : public ::testing::Test {
private:
    class VCFRecordMock : public VCFRecord {
    public:
        using VCFRecord::VCFRecord;
        MOCK_METHOD(std::string, alts_to_string, (), (const override));
        MOCK_METHOD(std::string, get_format,
            (bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage),
            (const override));
        MOCK_METHOD(std::string, sample_infos_to_string,
            (bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage),
            (const override));
    };

public:
    VCFRecordTest___to_string___Fixture()
        : vcf(create_VCF_with_default_parameters())
        , vcf_record(&vcf, "chrom", 10, "ref", "", "info", "graph_type_info")
    {
    }

    void SetUp() override { }

    void TearDown() override { }
    VCF vcf;
    VCFRecordMock vcf_record;
};

TEST_F(VCFRecordTest___to_string___Fixture, to_string)
{
    EXPECT_CALL(vcf_record, alts_to_string).Times(1).WillOnce(Return("alts"));
    EXPECT_CALL(vcf_record, get_format).Times(1).WillOnce(Return("format"));
    EXPECT_CALL(vcf_record, sample_infos_to_string)
        .Times(1)
        .WillOnce(Return("sample_infos"));

    std::string actual = vcf_record.to_string(true, false);
    std::string expected(
        "chrom\t11\t.\tref\talts\t.\t.\tinfo;graph_type_info\tformat\tsample_infos");

    EXPECT_EQ(actual, expected);
}

TEST(VCFRecordTest, get_longest_allele_length___longest_allele_is_ref)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf, "chrom1", 1, "ACGT", "A", ".", ".");

    size_t actual = vr.get_longest_allele_length();
    size_t expected = 4;

    EXPECT_EQ(actual, expected);
}

TEST(VCFRecordTest, get_longest_allele_length___longest_allele_is_alt)
{
    VCF vcf = create_VCF_with_default_parameters();
    VCFRecord vr(&vcf, "chrom1", 1, "ACGT", "A", ".", ".");
    vr.add_new_alt("ACGTTTTAC");
    vr.add_new_alt("C");

    size_t actual = vr.get_longest_allele_length();
    size_t expected = 9;

    EXPECT_EQ(actual, expected);
}

class VCFRecordTest___ref_allele_is_inside_given_interval______Fixture
    : public ::testing::Test {
public:
    VCFRecordTest___ref_allele_is_inside_given_interval______Fixture()
        : vcf(create_VCF_with_default_parameters())
        , vcf_record(&vcf, "chrom", 10, "REF", "")
    {
    }
    VCF vcf;
    VCFRecord vcf_record;
};
TEST_F(
    VCFRecordTest___ref_allele_is_inside_given_interval______Fixture, different_chroms)
{
    EXPECT_FALSE(vcf_record.ref_allele_is_inside_given_interval("dif_chrom", 0, 0));
}
TEST_F(VCFRecordTest___ref_allele_is_inside_given_interval______Fixture,
    ends_1_position_short)
{
    EXPECT_FALSE(vcf_record.ref_allele_is_inside_given_interval("chrom", 10, 12));
}
TEST_F(VCFRecordTest___ref_allele_is_inside_given_interval______Fixture,
    starts_1_position_after)
{
    EXPECT_FALSE(vcf_record.ref_allele_is_inside_given_interval("chrom", 11, 13));
}
TEST_F(VCFRecordTest___ref_allele_is_inside_given_interval______Fixture, exact_interval)
{
    EXPECT_TRUE(vcf_record.ref_allele_is_inside_given_interval("chrom", 10, 13));
}
TEST_F(
    VCFRecordTest___ref_allele_is_inside_given_interval______Fixture, larger_interval)
{
    EXPECT_TRUE(vcf_record.ref_allele_is_inside_given_interval("chrom", 9, 14));
}

class VCFRecordTest___merge_record_into_this______Fixture : public ::testing::Test {
public:
    VCFRecordTest___merge_record_into_this______Fixture()
        : vcf(create_VCF_with_default_parameters())
        , vcf_record_only_ref_no_alts(&vcf)
        , vcf_record_ref_A_alt_dot(&vcf, "1", 1, "A", ".")
        , vcf_record_ref_A_alt_T(&vcf, "1", 1, "A", "T")
        , vcf_record_ref_A_alt_T_dot(&vcf, "1", 1, "A", "T")
        , vcf_record_ref_A_alt_TTT(&vcf, "1", 1, "A", "TTT")
        , vcf_record_ref_A_alt_T_TT_TTT(&vcf, "1", 1, "A", "T")

    {
        vcf_record_ref_A_alt_T_dot.add_new_alt(".");

        vcf_record_ref_A_alt_T_TT_TTT.add_new_alt("TT");
        vcf_record_ref_A_alt_T_TT_TTT.add_new_alt("TTT");
    }
    VCF vcf;
    VCFRecord vcf_record_only_ref_no_alts;
    VCFRecord vcf_record_ref_A_alt_dot;
    VCFRecord vcf_record_ref_A_alt_T;
    VCFRecord vcf_record_ref_A_alt_T_dot;
    VCFRecord vcf_record_ref_A_alt_TTT;
    VCFRecord vcf_record_ref_A_alt_T_TT_TTT;
};

// TODO : mock this->sampleIndex_to_sampleInfo.merge_other_samples_infos_into_this()
TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_T_into_no_alts)
{
    vcf_record_only_ref_no_alts.merge_record_into_this(vcf_record_ref_A_alt_T);

    EXPECT_EQ(1, vcf_record_only_ref_no_alts.get_alts().size());
    EXPECT_EQ("T", vcf_record_only_ref_no_alts.get_alts()[0]);
    EXPECT_EQ(2,
        vcf_record_only_ref_no_alts.sampleIndex_to_sampleInfo[0]
            .get_number_of_alleles());
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_dot_into_no_alts)
{
    vcf_record_only_ref_no_alts.merge_record_into_this(vcf_record_ref_A_alt_dot);

    EXPECT_EQ(1, vcf_record_only_ref_no_alts.get_alts().size());
    EXPECT_EQ(".", vcf_record_only_ref_no_alts.get_alts()[0]);
    EXPECT_EQ(2,
        vcf_record_only_ref_no_alts.sampleIndex_to_sampleInfo[0]
            .get_number_of_alleles());
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_no_alts_into_T)
{
    vcf_record_ref_A_alt_T.merge_record_into_this(vcf_record_only_ref_no_alts);

    EXPECT_EQ(1, vcf_record_ref_A_alt_T.get_alts().size());
    EXPECT_EQ("T", vcf_record_ref_A_alt_T.get_alts()[0]);
    EXPECT_EQ(
        2, vcf_record_ref_A_alt_T.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_no_alts_into_dot)
{
    vcf_record_ref_A_alt_dot.merge_record_into_this(vcf_record_only_ref_no_alts);

    EXPECT_EQ(1, vcf_record_ref_A_alt_dot.get_alts().size());
    EXPECT_EQ(".", vcf_record_ref_A_alt_dot.get_alts()[0]);
    EXPECT_EQ(2,
        vcf_record_ref_A_alt_dot.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_dot_into_T)
{
    vcf_record_ref_A_alt_T.merge_record_into_this(vcf_record_ref_A_alt_dot);

    EXPECT_EQ(2, vcf_record_ref_A_alt_T.get_alts().size());
    EXPECT_EQ("T", vcf_record_ref_A_alt_T.get_alts()[0]);
    EXPECT_EQ(".", vcf_record_ref_A_alt_T.get_alts()[1]);
    EXPECT_EQ(
        3, vcf_record_ref_A_alt_T.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_TTT_into_T)
{
    vcf_record_ref_A_alt_T.merge_record_into_this(vcf_record_ref_A_alt_TTT);

    EXPECT_EQ(2, vcf_record_ref_A_alt_T.get_alts().size());
    EXPECT_EQ("T", vcf_record_ref_A_alt_T.get_alts()[0]);
    EXPECT_EQ("TTT", vcf_record_ref_A_alt_T.get_alts()[1]);
    EXPECT_EQ(
        3, vcf_record_ref_A_alt_T.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_T_into_TTT)
{
    vcf_record_ref_A_alt_TTT.merge_record_into_this(vcf_record_ref_A_alt_T);

    EXPECT_EQ(2, vcf_record_ref_A_alt_TTT.get_alts().size());
    EXPECT_EQ("TTT", vcf_record_ref_A_alt_TTT.get_alts()[0]);
    EXPECT_EQ("T", vcf_record_ref_A_alt_TTT.get_alts()[1]);
    EXPECT_EQ(3,
        vcf_record_ref_A_alt_TTT.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture, merge_T_dot_into_TTT)
{
    vcf_record_ref_A_alt_TTT.merge_record_into_this(vcf_record_ref_A_alt_T_dot);

    EXPECT_EQ(3, vcf_record_ref_A_alt_TTT.get_alts().size());
    EXPECT_EQ("TTT", vcf_record_ref_A_alt_TTT.get_alts()[0]);
    EXPECT_EQ("T", vcf_record_ref_A_alt_TTT.get_alts()[1]);
    EXPECT_EQ(".", vcf_record_ref_A_alt_TTT.get_alts()[2]);
    EXPECT_EQ(4,
        vcf_record_ref_A_alt_TTT.sampleIndex_to_sampleInfo[0].get_number_of_alleles());
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture,
    merge_first_alt_is_common___expects_FatalRuntimeError)
{
    ASSERT_EXCEPTION(
        vcf_record_ref_A_alt_T_TT_TTT.merge_record_into_this(vcf_record_ref_A_alt_T),
        FatalRuntimeError,
        "When merging two VCF records, they have common ALTs, this should not happen");
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture,
    merge_last_alt_is_common___expects_FatalRuntimeError)
{
    ASSERT_EXCEPTION(
        vcf_record_ref_A_alt_T_TT_TTT.merge_record_into_this(vcf_record_ref_A_alt_TTT),
        FatalRuntimeError,
        "When merging two VCF records, they have common ALTs, this should not happen");
}

TEST_F(VCFRecordTest___merge_record_into_this______Fixture,
    merge_both_have_dot_alleles___expects_FatalRuntimeError)
{
    ASSERT_EXCEPTION(
        vcf_record_ref_A_alt_dot.merge_record_into_this(vcf_record_ref_A_alt_T_dot),
        FatalRuntimeError,
        "When merging two VCF records, they have common ALTs, this should not happen");
}

class VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture
    : public ::testing::Test {
public:
    VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture()
        : vcf(create_VCF_with_default_parameters())
        , vcf_record_only_ref_no_alts(&vcf, "1", 1, "A", "")
        , vcf_record_tri_allelic(&vcf, "1", 1, "A", "T")
        , vcf_record_null_reference_alt_A(&vcf, "1", 1, ".", "T")
        , vcf_record_null_reference_alt_C(&vcf, "1", 1, ".", "C")
        , vcf_record_ref_A(&vcf, "1", 1, "A", "T")
        , vcf_record_ref_C(&vcf, "1", 1, "C", "CT")
        , vcf_record_ref_A_too_long_alt(&vcf, "1", 1, "A", "TTTTTT")
        , vcf_record_ref_A_same_chrom_different_pos(&vcf, "1", 2, "A", "T")
        , vcf_record_ref_A_different_chrom_same_pos(&vcf, "2", 1, "A", "T")
        , vcf_record_ref_A_different_alt(&vcf, "1", 1, "A", "C")
        , vcf_record_ref_A_same_alt(&vcf, "1", 1, "A", "T")
    {
        vcf_record_tri_allelic.add_new_alt("TT");
    }
    VCF vcf;
    VCFRecord vcf_record_only_ref_no_alts;
    VCFRecord vcf_record_tri_allelic;
    VCFRecord vcf_record_null_reference_alt_A;
    VCFRecord vcf_record_null_reference_alt_C;
    VCFRecord vcf_record_ref_A;
    VCFRecord vcf_record_ref_C;
    VCFRecord vcf_record_ref_A_too_long_alt;
    VCFRecord vcf_record_ref_A_same_chrom_different_pos;
    VCFRecord vcf_record_ref_A_different_chrom_same_pos;
    VCFRecord vcf_record_ref_A_different_alt;
    VCFRecord vcf_record_ref_A_same_alt;
};

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    merge_only_ref_no_alts___expects_FatalRuntimeError)
{
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(
        vcf_record_only_ref_no_alts);
    EXPECT_TRUE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    merge_triallelic___expects_FatalRuntimeError)
{
    ASSERT_EXCEPTION(vcf_record_ref_A.can_biallelic_record_be_merged_into_this(
                         vcf_record_tri_allelic),
        FatalRuntimeError,
        "When merging two biallelic records, one of them is not biallelic");
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    this_has_null_ref)
{
    bool actual
        = vcf_record_null_reference_alt_A.can_biallelic_record_be_merged_into_this(
            vcf_record_ref_A);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    other_has_null_ref)
{
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(
        vcf_record_null_reference_alt_A);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    both_have_null_ref)
{
    bool actual
        = vcf_record_null_reference_alt_A.can_biallelic_record_be_merged_into_this(
            vcf_record_null_reference_alt_C);
    EXPECT_TRUE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    different_refs)
{
    bool actual
        = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_ref_C);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    this_has_alt_too_long)
{
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(
        vcf_record_ref_A_too_long_alt, 5);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    other_has_alt_too_long)
{
    bool actual
        = vcf_record_ref_A_too_long_alt.can_biallelic_record_be_merged_into_this(
            vcf_record_ref_A, 5);
    EXPECT_FALSE(actual);
}

TEST_F(
    VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture, same_record)
{
    bool actual
        = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(vcf_record_ref_A);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    same_chrom_different_pos)
{
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(
        vcf_record_ref_A_same_chrom_different_pos);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    same_different_chrom_same_pos)
{
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(
        vcf_record_ref_A_different_chrom_same_pos);
    EXPECT_FALSE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    can_be_merged)
{
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(
        vcf_record_ref_A_different_alt);
    EXPECT_TRUE(actual);
}

TEST_F(VCFRecordTest___can_biallelic_record_be_merged_into_this______Fixture,
    same_alt_different_record)
{
    bool actual = vcf_record_ref_A.can_biallelic_record_be_merged_into_this(
        vcf_record_ref_A_same_alt);
    EXPECT_FALSE(actual);
}