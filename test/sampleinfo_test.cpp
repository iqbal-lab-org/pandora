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
using ::testing::DoubleNear;
using ::testing::DoubleEq;
using ::testing::Property;
using ::testing::_;


TEST(SampleInfoTest, constructor___zero_alleles___expects_death) {
    EXPECT_DEATH(SampleInfo(0, 0, &default_genotyping_options), "");
}

TEST(SampleInfoTest, constructor___one_allele) {
    SampleInfo sample_info(1, 1, &default_genotyping_options);

    EXPECT_EQ(sample_info.get_sample_index(), 1);
    EXPECT_EQ(sample_info.get_number_of_alleles(), 1);

    std::vector< std::vector<uint32_t> > default_coverages{{0}};
    EXPECT_ITERABLE_EQ(std::vector< std::vector<uint32_t> >, sample_info.get_allele_to_forward_coverages(),
                       default_coverages);
    EXPECT_ITERABLE_EQ(std::vector< std::vector<uint32_t> >, sample_info.get_allele_to_reverse_coverages(),
                       default_coverages);
    EXPECT_EQ(sample_info.get_genotyping_options(), &default_genotyping_options);
    EXPECT_EQ(sample_info.get_exp_depth_covg_for_this_sample(), default_genotyping_options.get_sample_index_to_exp_depth_covg()[1]);
}

TEST(SampleInfoTest, constructor___two_alleles) {
    SampleInfo sample_info(1, 2, &default_genotyping_options);

    EXPECT_EQ(sample_info.get_sample_index(), 1);
    EXPECT_EQ(sample_info.get_number_of_alleles(), 2);

    std::vector< std::vector<uint32_t> > default_coverages{{0}, {0}};
    EXPECT_ITERABLE_EQ(std::vector< std::vector<uint32_t> >, sample_info.get_allele_to_forward_coverages(),
                       default_coverages);
    EXPECT_ITERABLE_EQ(std::vector< std::vector<uint32_t> >, sample_info.get_allele_to_reverse_coverages(),
                       default_coverages);
    EXPECT_EQ(sample_info.get_genotyping_options(), &default_genotyping_options);
    EXPECT_EQ(sample_info.get_exp_depth_covg_for_this_sample(), default_genotyping_options.get_sample_index_to_exp_depth_covg()[1]);
}

TEST(SampleInfoTest, constructor___five_alleles) {
    SampleInfo sample_info(3, 5, &default_genotyping_options);

    EXPECT_EQ(sample_info.get_sample_index(), 3);
    EXPECT_EQ(sample_info.get_number_of_alleles(), 5);

    std::vector< std::vector<uint32_t> > default_coverages{{0}, {0}, {0}, {0}, {0}};
    EXPECT_ITERABLE_EQ(std::vector< std::vector<uint32_t> >, sample_info.get_allele_to_forward_coverages(),
                       default_coverages);
    EXPECT_ITERABLE_EQ(std::vector< std::vector<uint32_t> >, sample_info.get_allele_to_reverse_coverages(),
                       default_coverages);
    EXPECT_EQ(sample_info.get_genotyping_options(), &default_genotyping_options);
    EXPECT_EQ(sample_info.get_exp_depth_covg_for_this_sample(), default_genotyping_options.get_sample_index_to_exp_depth_covg()[3]);
}

class SampleInfoTest___Fixture : public ::testing::Test {
public:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        using SampleInfo::compute_likelihood;
    };

    SampleInfoTest___Fixture() :
            default_sample_info(0, 2, &default_genotyping_options),
            default_sample_info_two_alleles(0, 2, &default_genotyping_options),
            default_sample_info_three_alleles(0, 3, &default_genotyping_options) {}

    void SetUp() override {
        allele_to_coverage_one_allele = std::vector<std::vector<uint32_t> >({{1, 2}});

        allele_to_coverage_two_alleles = std::vector<std::vector<uint32_t>>({{1, 2}, {3, 4}});

        allele_to_coverage_three_alleles = std::vector<std::vector<uint32_t>>({{1, 2}, {3, 4}, {5, 6}});
    }

    void TearDown() override {
    }

    SampleInfoMock default_sample_info;
    SampleInfoMock default_sample_info_two_alleles;
    SampleInfoMock default_sample_info_three_alleles;
    std::vector<std::vector<uint32_t> > allele_to_coverage_empty;
    std::vector<std::vector<uint32_t> > allele_to_coverage_one_allele;
    std::vector<std::vector<uint32_t> > allele_to_coverage_two_alleles;
    std::vector<std::vector<uint32_t> > allele_to_coverage_three_alleles;
};

TEST_F(SampleInfoTest___Fixture, get_and_set_gt_from_max_likelihood_path___set_a_valid_gt) {
    default_sample_info.set_gt_from_max_likelihood_path(5);

    uint32_t actual = default_sample_info.get_gt_from_max_likelihood_path();
    uint32_t expected = 5;

    EXPECT_EQ(actual, expected);
}

TEST_F(SampleInfoTest___Fixture,
       get_and_set_gt_from_max_likelihood_path___set_a_valid_then_an_invalid_gt___throws_exception) {
    default_sample_info.set_gt_from_max_likelihood_path(5);
    default_sample_info.set_gt_from_max_likelihood_path(boost::none);

    EXPECT_FALSE(default_sample_info.is_gt_from_max_likelihood_path_valid());
    try {
        default_sample_info.get_gt_from_max_likelihood_path();
        FAIL() << "Should have thrown exception";
    } catch (std::runtime_error &error) {}
}

TEST_F(SampleInfoTest___Fixture, is_gt_from_max_likelihood_path_valid___default_sample_info___is_invalid) {
    EXPECT_FALSE(default_sample_info.is_gt_from_max_likelihood_path_valid());
}

TEST_F(SampleInfoTest___Fixture, is_gt_from_max_likelihood_path_valid___set_valid_gt___is_valid) {
    default_sample_info.set_gt_from_max_likelihood_path(5);
    EXPECT_TRUE(default_sample_info.is_gt_from_max_likelihood_path_valid());
}

TEST_F(SampleInfoTest___Fixture, is_gt_from_max_likelihood_path_valid___set_valid_then_invalid_gt___is_invalid) {
    default_sample_info.set_gt_from_max_likelihood_path(5);
    default_sample_info.set_gt_from_max_likelihood_path(boost::none);
    EXPECT_FALSE(default_sample_info.is_gt_from_max_likelihood_path_valid());
}


TEST_F(SampleInfoTest___Fixture, get_gt_from_max_likelihood_path___default_sample_info___throws_exception) {
    try {
        default_sample_info.get_gt_from_max_likelihood_path();
        FAIL() << "Should have thrown exception";
    } catch (std::runtime_error &error) {}
}

TEST_F(SampleInfoTest___Fixture, gt_from_max_likelihood_path_to_string___invalid_gt) {
    std::string actual = default_sample_info.gt_from_max_likelihood_path_to_string();
    std::string expected = ".";

    EXPECT_EQ(actual, expected);
}

TEST_F(SampleInfoTest___Fixture, gt_from_max_likelihood_path_to_string___valid_gt) {
    default_sample_info.set_gt_from_max_likelihood_path(5);

    std::string actual = default_sample_info.gt_from_max_likelihood_path_to_string();
    std::string expected = "5";

    EXPECT_EQ(actual, expected);
}

TEST_F(SampleInfoTest___Fixture, get_allele_to_forward_coverages___default_sample_info) {
    std::vector<std::vector<uint32_t> > actual = default_sample_info.get_allele_to_forward_coverages();
    std::vector<std::vector<uint32_t> > expected{ { 0 }, { 0 } };

    EXPECT_EQ(actual, expected);
}

TEST_F(SampleInfoTest___Fixture, get_allele_to_reverse_coverages___default_sample_info) {
    std::vector<std::vector<uint32_t> > actual = default_sample_info.get_allele_to_reverse_coverages();
    std::vector<std::vector<uint32_t> > expected{ { 0 }, { 0 } };

    EXPECT_EQ(actual, expected);
}

TEST_F(SampleInfoTest___Fixture, set_coverage_information___forward_coverage_has_no_alleles___expects_death) {
    EXPECT_DEATH(
            default_sample_info.set_coverage_information(allele_to_coverage_empty, allele_to_coverage_three_alleles),
            "");
}

TEST_F(SampleInfoTest___Fixture, set_coverage_information___forward_coverage_has_one_allele___expects_death) {
    EXPECT_DEATH(default_sample_info.set_coverage_information(allele_to_coverage_one_allele,
                                                              allele_to_coverage_three_alleles), "");
}

TEST_F(SampleInfoTest___Fixture, set_coverage_information___reverse_coverage_has_no_alleles___expects_death) {
    EXPECT_DEATH(
            default_sample_info.set_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_empty),
            "");
}

TEST_F(SampleInfoTest___Fixture, set_coverage_information___reverse_coverage_has_one_allele___expects_death) {
    EXPECT_DEATH(default_sample_info.set_coverage_information(allele_to_coverage_three_alleles,
                                                              allele_to_coverage_one_allele), "");
}

TEST_F(SampleInfoTest___Fixture,
       set_coverage_information___both_coverages_have_two_alleles___different_number_of_bases___expects_death) {
    EXPECT_DEATH(default_sample_info.set_coverage_information(allele_to_coverage_two_alleles, {{1, 2},
                                                                                               {3}}), "");
}

TEST_F(SampleInfoTest___Fixture, set_coverage_information___both_coverages_have_two_alleles) {
    default_sample_info.set_coverage_information(allele_to_coverage_two_alleles, allele_to_coverage_two_alleles);

    auto expected = allele_to_coverage_two_alleles;

    auto actual = default_sample_info.get_allele_to_forward_coverages();
    EXPECT_EQ(actual, expected);

    actual = default_sample_info.get_allele_to_reverse_coverages();
    EXPECT_EQ(actual, expected);
}

TEST_F(SampleInfoTest___Fixture, set_coverage_information___fwd_coverage_has_two_alleles___rev_coverage_has_three_alleles___sample_info_expects_three_alleles___expects_death) {
    EXPECT_DEATH(default_sample_info_three_alleles.set_coverage_information(allele_to_coverage_two_alleles, allele_to_coverage_three_alleles), "");
}

TEST_F(SampleInfoTest___Fixture, set_coverage_information___fwd_coverage_has_three_alleles___rev_coverage_has_two_alleles___sample_info_expects_three_alleles___expects_death) {
    EXPECT_DEATH(default_sample_info_three_alleles.set_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_two_alleles), "");
}

TEST_F(SampleInfoTest___Fixture, set_coverage_information___both_coverages_have_three_alleles) {
    default_sample_info_three_alleles.set_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_three_alleles);

    auto expected = allele_to_coverage_three_alleles;
    auto actual = default_sample_info_three_alleles.get_allele_to_forward_coverages();
    EXPECT_EQ(actual, expected);

    actual = default_sample_info_three_alleles.get_allele_to_reverse_coverages();
    EXPECT_EQ(actual, expected);
}


TEST_F(SampleInfoTest___Fixture,
       set_coverage_information___forward_covg_has_two_alleles___reverse_covg_has_three_alleles___expects_death) {
    EXPECT_DEATH(default_sample_info.set_coverage_information(allele_to_coverage_two_alleles,
                                                              allele_to_coverage_three_alleles), "");
}

TEST_F(SampleInfoTest___Fixture, set_number_of_alleles_and_resize_coverage_information___resize_to_zero_alleles___expects_death) {
    default_sample_info_three_alleles.set_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_three_alleles);
    EXPECT_DEATH(default_sample_info_three_alleles.set_number_of_alleles_and_resize_coverage_information(0), "");
}

TEST_F(SampleInfoTest___Fixture, set_number_of_alleles_and_resize_coverage_information___resize_to_one_allele___shrinks) {
    default_sample_info_three_alleles.set_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_three_alleles);

    default_sample_info_three_alleles.set_number_of_alleles_and_resize_coverage_information(1);

    EXPECT_EQ(1, default_sample_info_three_alleles.get_number_of_alleles());
    EXPECT_EQ(allele_to_coverage_one_allele, default_sample_info_three_alleles.get_allele_to_forward_coverages());
    EXPECT_EQ(allele_to_coverage_one_allele, default_sample_info_three_alleles.get_allele_to_reverse_coverages());
}

TEST_F(SampleInfoTest___Fixture, set_number_of_alleles_and_resize_coverage_information___resize_to_two_alleles___shrinks) {
    default_sample_info_three_alleles.set_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_three_alleles);

    default_sample_info_three_alleles.set_number_of_alleles_and_resize_coverage_information(2);

    EXPECT_EQ(2, default_sample_info_three_alleles.get_number_of_alleles());
    EXPECT_EQ(allele_to_coverage_two_alleles, default_sample_info_three_alleles.get_allele_to_forward_coverages());
    EXPECT_EQ(allele_to_coverage_two_alleles, default_sample_info_three_alleles.get_allele_to_reverse_coverages());
}

TEST_F(SampleInfoTest___Fixture, set_number_of_alleles_and_resize_coverage_information___resize_to_three_alleles___no_change) {
    default_sample_info_three_alleles.set_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_three_alleles);

    default_sample_info_three_alleles.set_number_of_alleles_and_resize_coverage_information(3);

    EXPECT_EQ(3, default_sample_info_three_alleles.get_number_of_alleles());
    EXPECT_EQ(allele_to_coverage_three_alleles, default_sample_info_three_alleles.get_allele_to_forward_coverages());
    EXPECT_EQ(allele_to_coverage_three_alleles, default_sample_info_three_alleles.get_allele_to_reverse_coverages());
}

TEST_F(SampleInfoTest___Fixture, set_number_of_alleles_and_resize_coverage_information___resize_to_four_alleles___expands) {
    default_sample_info_three_alleles.set_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_three_alleles);

    default_sample_info_three_alleles.set_number_of_alleles_and_resize_coverage_information(4);

    EXPECT_EQ(4, default_sample_info_three_alleles.get_number_of_alleles());
    auto allele_to_coverage_four_alleles = allele_to_coverage_three_alleles;
    allele_to_coverage_four_alleles.push_back({0});
    EXPECT_EQ(allele_to_coverage_four_alleles, default_sample_info_three_alleles.get_allele_to_forward_coverages());
    EXPECT_EQ(allele_to_coverage_four_alleles, default_sample_info_three_alleles.get_allele_to_reverse_coverages());
}



TEST_F(SampleInfoTest___Fixture, gt_coverages_compatible___default_sample_info___invalid_gt) {
    EXPECT_FALSE(default_sample_info.is_gt_from_coverages_valid());
    try {
        default_sample_info.get_gt_from_coverages();
        FAIL() << "Should have thrown exception";
    } catch (const std::runtime_error &e) {}
    EXPECT_EQ(".", default_sample_info.gt_from_coverages_compatible_to_string());
}


TEST_F(SampleInfoTest___Fixture, gt_coverages_compatible___default_sample_info___valid_gt) {
    default_sample_info.set_gt_from_coverages_compatible(5);
    EXPECT_TRUE(default_sample_info.is_gt_from_coverages_compatible_valid());
    EXPECT_EQ(5, default_sample_info.get_gt_from_coverages_compatible());
    EXPECT_EQ("5", default_sample_info.gt_from_coverages_compatible_to_string());
}

class SampleInfoTest___genotype_from_coverage___Fixture : public ::testing::Test {
private:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(bool, check_if_coverage_information_is_correct, (), (const override));
        MOCK_METHOD(boost::optional<GenotypeAndMaxLikelihood>, get_genotype_from_coverage, (), (const override));
    };

public:
    SampleInfoTest___genotype_from_coverage___Fixture() :
            default_sample_info(0, 2, &default_genotyping_options) {}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfoMock default_sample_info;
};

TEST_F(SampleInfoTest___genotype_from_coverage___Fixture, valid_genotype) {
    EXPECT_CALL(default_sample_info, check_if_coverage_information_is_correct)
            .Times(1)
            .WillOnce(Return(true));

    EXPECT_CALL(default_sample_info, get_genotype_from_coverage).
            Times(1).
            WillOnce(Return(std::make_pair<uint32_t, double>(2, -1.0)));

    default_sample_info.genotype_from_coverage();

    EXPECT_TRUE(default_sample_info.is_gt_from_coverages_valid());
    EXPECT_EQ(2, default_sample_info.get_gt_from_coverages());
    EXPECT_NEAR(-1.0, default_sample_info.get_likelihood_of_gt_from_coverages(), 0.000001);
    EXPECT_EQ(2, default_sample_info.get_gt_from_coverages_compatible());
}

TEST_F(SampleInfoTest___genotype_from_coverage___Fixture, invalid_genotype) {
    EXPECT_CALL(default_sample_info, check_if_coverage_information_is_correct)
            .Times(1)
            .WillOnce(Return(true));

    EXPECT_CALL(default_sample_info, get_genotype_from_coverage).
            Times(1).
            WillOnce(Return(boost::none));

    default_sample_info.genotype_from_coverage();

    EXPECT_FALSE(default_sample_info.is_gt_from_coverages_valid());
    try {
        default_sample_info.get_gt_from_coverages();
        FAIL() << "Should have thrown exception";
    } catch (const std::runtime_error &e) {}
    try {
        default_sample_info.get_likelihood_of_gt_from_coverages();
        FAIL() << "Should have thrown exception";
    } catch (const std::runtime_error &e) {}
}


class SampleInfoTest___merge_other_sample_info_into_this___Fixture : public ::testing::Test {
private:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(boost::optional<SampleInfo::GenotypeAndMaxLikelihood>, get_genotype_from_coverage, (), (const override));
    };


public:
    SampleInfoTest___merge_other_sample_info_into_this___Fixture() :
            sample_info_with_two_alleles(0, 2, &default_genotyping_options),
            sample_info_with_three_alleles(0, 3, &default_genotyping_options) {}

    void SetUp() override {
        allele_to_coverage_two_alleles.push_back({1, 2});
        allele_to_coverage_two_alleles.push_back({3, 4});
        sample_info_with_two_alleles.set_coverage_information(allele_to_coverage_two_alleles,
                                                              allele_to_coverage_two_alleles);

        allele_to_coverage_three_alleles.push_back({1, 2});
        allele_to_coverage_three_alleles.push_back({5, 6});
        allele_to_coverage_three_alleles.push_back({7, 8});
        sample_info_with_three_alleles.set_coverage_information(allele_to_coverage_three_alleles,
                                                                allele_to_coverage_three_alleles);
    }

    void TearDown() override {
    }

    SampleInfoMock sample_info_with_two_alleles;
    SampleInfoMock sample_info_with_three_alleles;
    std::vector<std::vector<uint32_t> > allele_to_coverage_two_alleles;
    std::vector<std::vector<uint32_t> > allele_to_coverage_three_alleles;
};

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture,
       merge_a_sample_with_three_alleles_into_one_with_two_alleles___same_gt_from_max_likelihood_path) {
    sample_info_with_two_alleles.set_gt_from_max_likelihood_path(0);
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(0);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    std::vector<std::vector<uint32_t> > allele_to_coverage_expected({{1, 2},
                                                                     {3, 4},
                                                                     {5, 6},
                                                                     {7, 8}});
    EXPECT_EQ(allele_to_coverage_expected, sample_info_with_two_alleles.get_allele_to_forward_coverages());
    EXPECT_EQ(allele_to_coverage_expected, sample_info_with_two_alleles.get_allele_to_reverse_coverages());
    EXPECT_EQ(0, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture,
       merge_a_sample_with_three_alleles_into_one_with_two_alleles___both_gt_from_max_likelihood_path_are_invalid) {
    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_FALSE(sample_info_with_two_alleles.is_gt_from_max_likelihood_path_valid());
}

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture,
       merge_a_sample_with_three_alleles_into_one_with_two_alleles___gt_merged_in_is_invalid) {
    sample_info_with_two_alleles.set_gt_from_max_likelihood_path(1);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(1, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}


TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture,
       merge_a_sample_with_three_alleles_into_one_with_two_alleles___original_gt_invalid___gt_merged_in_is_valid_zero) {
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(0);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(0, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture,
       merge_a_sample_with_three_alleles_into_one_with_two_alleles___original_gt_invalid___gt_merged_in_is_valid_one) {
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(1);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(2, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture,
       merge_a_sample_with_three_alleles_into_one_with_two_alleles___original_gt_invalid___gt_merged_in_is_valid_two) {
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(2);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(3, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture,
       merge_a_sample_with_three_alleles_into_one_with_two_alleles___both_gts_are_valid_and_first_is_not_zero___genotypes_from_coverage_to_solve_conflict) {
    EXPECT_CALL(sample_info_with_two_alleles, get_genotype_from_coverage()).
            Times(1).
            WillOnce(Return(SampleInfo::GenotypeAndMaxLikelihood(1, -1.0)));
    sample_info_with_two_alleles.set_gt_from_max_likelihood_path(1);
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(0);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(1, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}


TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture,
       merge_a_sample_with_three_alleles_into_one_with_two_alleles___both_gts_are_valid_and_second_is_not_zero___genotypes_from_coverage_to_solve_conflict) {
    EXPECT_CALL(sample_info_with_two_alleles, get_genotype_from_coverage()).
            Times(1).
            WillOnce(Return(SampleInfo::GenotypeAndMaxLikelihood(3, -1.0)));
    sample_info_with_two_alleles.set_gt_from_max_likelihood_path(0);
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(2);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(3, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}


class SampleInfoTest___get_gaps___Fixture : public ::testing::Test {
public:
    static GenotypingOptions genotyping_options_with_min_kmer_covg_10;

    SampleInfoTest___get_gaps___Fixture() :
            sample_info_with_min_kmer_covg_10(0, 2, &genotyping_options_with_min_kmer_covg_10) {}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfo sample_info_with_min_kmer_covg_10;
};

GenotypingOptions SampleInfoTest___get_gaps___Fixture::genotyping_options_with_min_kmer_covg_10(
        {10, 10, 10, 10, 10, 10, 10, 10, 10, 10}, 0.01, 0, 0, 0, 0, 0, 10, false);

TEST_F(SampleInfoTest___get_gaps___Fixture, get_gaps___coverages_just_below_threshold) {
    sample_info_with_min_kmer_covg_10.set_coverage_information({{5},
                                                                {5, 5, 5, 5}}, {{4},
                                                                                {4, 4, 4, 4}});

    double actual_allele_0 = sample_info_with_min_kmer_covg_10.get_gaps(0);
    double actual_allele_1 = sample_info_with_min_kmer_covg_10.get_gaps(1);

    EXPECT_TRUE(Maths::equals(1.0, actual_allele_0));
    EXPECT_TRUE(Maths::equals(1.0, actual_allele_1));
}

TEST_F(SampleInfoTest___get_gaps___Fixture, get_gaps___coverages_all_equal_threshold) {
    sample_info_with_min_kmer_covg_10.set_coverage_information({{5},
                                                                {5, 5, 5, 5}}, {{5},
                                                                                {5, 5, 5, 5}});

    double actual_allele_0 = sample_info_with_min_kmer_covg_10.get_gaps(0);
    double actual_allele_1 = sample_info_with_min_kmer_covg_10.get_gaps(1);

    EXPECT_TRUE(Maths::equals(0.0, actual_allele_0));
    EXPECT_TRUE(Maths::equals(0.0, actual_allele_1));
}

TEST_F(SampleInfoTest___get_gaps___Fixture, get_gaps___coverages_all_above_threshold) {
    sample_info_with_min_kmer_covg_10.set_coverage_information({{5},
                                                                {5, 5, 5, 5}}, {{6},
                                                                                {6, 6, 6, 6}});

    double actual_allele_0 = sample_info_with_min_kmer_covg_10.get_gaps(0);
    double actual_allele_1 = sample_info_with_min_kmer_covg_10.get_gaps(1);

    EXPECT_TRUE(Maths::equals(0.0, actual_allele_0));
    EXPECT_TRUE(Maths::equals(0.0, actual_allele_1));
}

TEST_F(SampleInfoTest___get_gaps___Fixture, get_gaps___coverages_below_equal_and_above_threshold) {
    sample_info_with_min_kmer_covg_10.set_coverage_information({{0, 0,  0},
                                                                {9, 10, 11, 9, 10, 9}}, {{9, 10, 11},
                                                                                         {0, 0,  0, 0, 0, 0}});

    double actual_allele_0 = sample_info_with_min_kmer_covg_10.get_gaps(0);
    double actual_allele_1 = sample_info_with_min_kmer_covg_10.get_gaps(1);

    EXPECT_TRUE(Maths::equals(1.0 / 3.0, actual_allele_0));
    EXPECT_TRUE(Maths::equals(0.5, actual_allele_1));
}


TEST(SampleInfoTest, get_min_coverage_threshold_for_this_sample___min_allele_covg_is_higher) {
    GenotypingOptions genotyping_options(
            {10, 5}, 0.01, 0, 100, 1.0, 0, 0, 00, false);
    SampleInfo sample_info(0, 2, &genotyping_options);

    EXPECT_EQ(100, sample_info.get_min_coverage_threshold_for_this_sample());
}

TEST(SampleInfoTest, get_min_coverage_threshold_for_this_sample___min_fraction_allele_covg_is_higher) {
    GenotypingOptions genotyping_options(
            {10, 100}, 0.01, 0, 40, 0.5, 0, 0, 00, false);
    SampleInfo sample_info(1, 2, &genotyping_options);

    EXPECT_EQ(50, sample_info.get_min_coverage_threshold_for_this_sample());
}


class SampleInfoTest___get_likelihoods_for_all_alleles___Fixture : public ::testing::Test {
public:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;

        MOCK_METHOD(size_t, get_number_of_alleles, (), (const override));
        MOCK_METHOD(uint32_t, get_min_coverage_threshold_for_this_sample, (), (const override));
        MOCK_METHOD(uint32_t, get_total_mean_coverage_over_all_alleles_given_a_minimum_threshold, (uint32_t
                minimum_threshold), (const override));
        MOCK_METHOD(uint32_t, get_total_mean_coverage_given_a_minimum_threshold, (uint32_t
                allele, uint32_t
                minimum_threshold), (const override));
        MOCK_METHOD(double, get_gaps, (uint32_t
                allele), (const override));
        MOCK_METHOD(double, compute_likelihood, (bool
                min_coverage_threshold_is_satisfied, double
                expected_depth_covg, double
                total_mean_coverage_of_allele_above_threshold,
                        double
                total_mean_coverage_of_all_other_alleles_above_threshold, double
                error_rate, double
                gaps), (const override));
    };

    SampleInfoTest___get_likelihoods_for_all_alleles___Fixture() :
            sample_info(0, 2, &default_genotyping_options) {}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfoMock sample_info;
};

TEST_F(SampleInfoTest___get_likelihoods_for_all_alleles___Fixture, get_likelihoods_for_all_alleles) {
    EXPECT_CALL(sample_info, get_number_of_alleles)
            .WillRepeatedly(Return(2));

    EXPECT_CALL(sample_info, get_min_coverage_threshold_for_this_sample)
            .Times(1)
            .WillOnce(Return(5));

    EXPECT_CALL(sample_info, get_total_mean_coverage_over_all_alleles_given_a_minimum_threshold(5))
            .Times(1)
            .WillOnce(Return(100));

    //allele 0 - min coverage threshold not satisfied
    EXPECT_CALL(sample_info, get_total_mean_coverage_given_a_minimum_threshold(0, 5))
            .Times(1)
            .WillOnce(Return(0));

    EXPECT_CALL(sample_info, get_gaps(0))
            .Times(1)
            .WillOnce(Return(0.3));

    EXPECT_CALL(sample_info,
                compute_likelihood(false, default_genotyping_options.get_sample_index_to_exp_depth_covg()[0],
                                   0, 100, DoubleEq(default_genotyping_options.get_error_rate()), DoubleEq(0.3)))
            .Times(1)
            .WillOnce(Return(-50.5));

    //allele 1 - min coverage threshold not satisfied
    EXPECT_CALL(sample_info, get_total_mean_coverage_given_a_minimum_threshold(1, 5))
            .Times(1)
            .WillOnce(Return(100));

    EXPECT_CALL(sample_info, get_gaps(1))
            .Times(1)
            .WillOnce(Return(0.5));

    EXPECT_CALL(sample_info,
                compute_likelihood(true, DoubleEq(default_genotyping_options.get_sample_index_to_exp_depth_covg()[1]),
                                   DoubleEq(100.0), DoubleEq(0.0),
                                   DoubleEq(default_genotyping_options.get_error_rate()), DoubleEq(0.5)))
            .Times(1)
            .WillOnce(Return(-2));

    std::vector<double> actual = sample_info.get_likelihoods_for_all_alleles();

    std::vector<double> expected = {-50.5, -2};
    EXPECT_ITERABLE_EQ(std::vector<double>, actual, expected);
}


class SampleInfoTest___get_confidence___Fixture : public ::testing::Test {
public:
    static GenotypingOptions genotyping_options_get_confidence;

    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(std::vector<double>, get_likelihoods_for_all_alleles, (), (const override));
        MOCK_METHOD(uint32_t, get_mean_coverage_both_alleles, (uint32_t
                allele), (const override));
    };

    SampleInfoTest___get_confidence___Fixture() :
            sample_info(0, 2, &genotyping_options_get_confidence), default_sample_info(0, 2, &default_genotyping_options) {}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfoMock sample_info;
    SampleInfoMock default_sample_info;
};

GenotypingOptions SampleInfoTest___get_confidence___Fixture::genotyping_options_get_confidence
        ({10, 10}, 0.01, 0, 0, 0.0, 50, 100, 0, 0);

TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___not_enough_total_covg) {
    std::vector<double> likelihoods({-10, -3, -5});
    EXPECT_CALL(sample_info, get_likelihoods_for_all_alleles)
            .Times(1)
            .WillOnce(Return(likelihoods));

    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(1))
            .Times(1)
            .WillOnce(Return(30));

    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(2))
            .Times(1)
            .WillOnce(Return(10));

    auto actual = sample_info.get_confidence();

    EXPECT_EQ(boost::none, actual);
}

TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___not_enough_difference_in_covg) {
    std::vector<double> likelihoods({-10, -3, -5});
    EXPECT_CALL(sample_info, get_likelihoods_for_all_alleles)
            .Times(1)
            .WillOnce(Return(likelihoods));

    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(1))
            .Times(1)
            .WillOnce(Return(100));

    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(2))
            .Times(1)
            .WillOnce(Return(199));

    auto actual = sample_info.get_confidence();

    EXPECT_EQ(boost::none, actual);
}


TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___returns_valid_confidence) {
    std::vector<double> likelihoods({-10, -3, -5});
    EXPECT_CALL(sample_info, get_likelihoods_for_all_alleles)
            .Times(1)
            .WillOnce(Return(likelihoods));

    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(1))
            .Times(1)
            .WillOnce(Return(100));

    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(2))
            .Times(1)
            .WillOnce(Return(200));

    auto actual = sample_info.get_confidence();

    EXPECT_EQ(1, std::get<0>(*actual));
    EXPECT_NEAR(2.0, std::get<1>(*actual), 0.000001);
    EXPECT_NEAR(-3.0, std::get<2>(*actual), 0.000001);
}


class SampleInfoTest___get_confidence_to_string___Fixture : public ::testing::Test {
public:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(boost::optional<IndexAndConfidenceAndMaxLikelihood>, get_confidence, (), (const override));
    };

    SampleInfoTest___get_confidence_to_string___Fixture() :
            sample_info(0, 2, &default_genotyping_options) {}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfoMock sample_info;
};

TEST_F(SampleInfoTest___get_confidence_to_string___Fixture, invalid_confidence) {
    EXPECT_CALL(sample_info, get_confidence)
            .Times(1)
            .WillOnce(Return(boost::none));

    EXPECT_EQ(".", sample_info.get_confidence_to_string());
}


TEST_F(SampleInfoTest___get_confidence_to_string___Fixture, valid_confidence) {
    EXPECT_CALL(sample_info, get_confidence)
            .Times(1)
            .WillOnce(Return(std::make_tuple((size_t) 1, 50.5, -200.2)));

    EXPECT_EQ("50.5", sample_info.get_confidence_to_string());
}


class SampleInfoTest___get_genotype_from_coverage___Fixture : public ::testing::Test {
public:
    static GenotypingOptions genotyping_options_high_confidence_threshold;

    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(boost::optional<SampleInfo::IndexAndConfidenceAndMaxLikelihood>, get_confidence, (), (const override));
    };

    SampleInfoTest___get_genotype_from_coverage___Fixture() :
            sample_info(0, 2, &genotyping_options_high_confidence_threshold) {}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfoMock sample_info;
};

GenotypingOptions SampleInfoTest___get_genotype_from_coverage___Fixture::genotyping_options_high_confidence_threshold
        ({10, 10}, 0.01, 100, 0, 0.0, 0, 0, 0, 0);

TEST_F(SampleInfoTest___get_genotype_from_coverage___Fixture, invalid_confidence) {
    EXPECT_CALL(sample_info, get_confidence)
            .Times(1)
            .WillOnce(Return(boost::none));

    auto actual = sample_info.get_genotype_from_coverage();

    EXPECT_EQ(boost::none, actual);
}


TEST_F(SampleInfoTest___get_genotype_from_coverage___Fixture, valid_confidence_but_below_threshold) {
    EXPECT_CALL(sample_info, get_confidence)
            .Times(1)
            .WillOnce(Return(std::make_tuple((size_t) 1, 99.0, -50.5)));

    auto actual = sample_info.get_genotype_from_coverage();

    EXPECT_EQ(boost::none, actual);
}

TEST_F(SampleInfoTest___get_genotype_from_coverage___Fixture, valid_confidence_and_above_threshold) {
    EXPECT_CALL(sample_info, get_confidence)
            .Times(1)
            .WillOnce(Return(std::make_tuple((size_t) 1, 105.0, -50.5)));

    auto actual = sample_info.get_genotype_from_coverage();

    EXPECT_EQ(1, actual->first);
    EXPECT_NEAR(-50.5, actual->second, 0.000001);
}


TEST_F(SampleInfoTest___Fixture, to_string___no_flags_set___expects_death) {
    EXPECT_DEATH(default_sample_info.to_string(false, false), "");
}

TEST_F(SampleInfoTest___Fixture, to_string___both_flags_set___expects_death) {
    EXPECT_DEATH(default_sample_info.to_string(true, true), "");
}

TEST_F(SampleInfoTest___Fixture, to_string___genotyping_from_maximum_likelihood) {
    default_sample_info_three_alleles.set_gt_from_max_likelihood_path(1);
    default_sample_info_three_alleles.set_coverage_information({{10},
                                                  {20, 30},
                                                  {40, 50, 70}},
                                                 {{70},
                                                  {80,  90},
                                                  {100, 120, 130}});

    std::string actual = default_sample_info_three_alleles.to_string(true, false);

    std::string expected = "1:10,25,53:70,85,116:10,25,50:70,85,120:10,50,160:70,170,350:0,0,0";
    EXPECT_EQ(actual, expected);
}

TEST_F(SampleInfoTest___Fixture, to_string___genotyping_from_compatible_coverage) {
    default_sample_info_three_alleles.set_gt_from_coverages_compatible(2);
    default_sample_info_three_alleles.set_coverage_information({{10},
                                                  {20, 30},
                                                  {40, 50, 70}},
                                                 {{70},
                                                  {80,  90},
                                                  {100, 120, 130}});

    std::string actual = default_sample_info_three_alleles.to_string(false, true);

    std::string expected = "2:10,25,53:70,85,116:10,25,50:70,85,120:10,50,160:70,170,350:0,0,0:-1559.97,-1558.47,-1577.88:1.50545";
    EXPECT_EQ(actual, expected);
}


TEST_F(SampleInfoTest___Fixture, compute_likelihood___min_coverage_threshold_is_satisfied) {
    bool min_coverage_threshold_is_satisfied = true;
    double expected_depth_covg = 10;
    double total_mean_coverage_of_allele_above_threshold = 4;
    double total_mean_coverage_of_all_other_alleles_above_threshold = 100;
    double error_rate = 0.05;
    double gaps = 0.02;
    double actual = default_sample_info.compute_likelihood(min_coverage_threshold_is_satisfied, expected_depth_covg,
                                                           total_mean_coverage_of_allele_above_threshold,
                                                           total_mean_coverage_of_all_other_alleles_above_threshold,
                                                           error_rate, gaps);

    double expected = -expected_depth_covg
                      + total_mean_coverage_of_allele_above_threshold * log(expected_depth_covg)
                      - Maths::logfactorial(total_mean_coverage_of_allele_above_threshold)
                      + total_mean_coverage_of_all_other_alleles_above_threshold * log(error_rate)
                      - expected_depth_covg * gaps
                      + log(1 - exp(-expected_depth_covg)) * (1 - gaps);
    EXPECT_DOUBLE_EQ(actual, expected);
}


TEST_F(SampleInfoTest___Fixture, compute_likelihood___min_coverage_threshold_is_not_satisfied) {
    bool min_coverage_threshold_is_satisfied = false;
    double expected_depth_covg = 10;
    double total_mean_coverage_of_allele_above_threshold = 4;
    double total_mean_coverage_of_all_other_alleles_above_threshold = 100;
    double error_rate = 0.05;
    double gaps = 0.02;
    double actual = default_sample_info.compute_likelihood(min_coverage_threshold_is_satisfied, expected_depth_covg,
                                                           total_mean_coverage_of_allele_above_threshold,
                                                           total_mean_coverage_of_all_other_alleles_above_threshold,
                                                           error_rate, gaps);

    double expected = -expected_depth_covg
                      + total_mean_coverage_of_all_other_alleles_above_threshold * log(error_rate)
                      - expected_depth_covg * gaps
                      + log(1 - exp(-expected_depth_covg)) * (1 - gaps);
    EXPECT_DOUBLE_EQ(actual, expected);
}



















template<class SAMPLE_TYPE>
class SampleIndexToSampleInfoTemplateAllVisible : public SampleIndexToSampleInfoTemplate<SAMPLE_TYPE> {
public:
    using SampleIndexToSampleInfoTemplate<SAMPLE_TYPE>::SampleIndexToSampleInfoTemplate;
    using SampleIndexToSampleInfoTemplate<SAMPLE_TYPE>::emplace_back_several_empty_sample_infos;
    using SampleIndexToSampleInfoTemplate<SAMPLE_TYPE>::set_number_of_alleles_and_resize_coverage_information_for_all_samples;
    using SampleIndexToSampleInfoTemplate<SAMPLE_TYPE>::merge_other_samples_infos_into_this;
};

class SampleIndexToSampleInfoTemplate___Fixture : public ::testing::Test {
public:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;

        SampleInfoMock(const SampleInfoMock &other) : SampleInfoMock(other.get_sample_index(), 2,
                                                                     other.genotyping_options) {}

        MOCK_METHOD(void, merge_other_sample_info_into_this, (const SampleInfo &other), (override));
        MOCK_METHOD(std::string, to_string, (bool
                genotyping_from_maximum_likelihood, bool
                genotyping_from_compatible_coverage), (const override));
    };


    void SetUp() override {
        sample_index_to_sample_info.emplace_back_several_empty_sample_infos(2, 2, &default_genotyping_options);
    }

    void TearDown() override {
    }

    SampleIndexToSampleInfoTemplateAllVisible<SampleInfoMock> sample_index_to_sample_info;
};

TEST_F(SampleIndexToSampleInfoTemplate___Fixture, emplace_back_several_empty_sample_infos___same_number_of_alleles) {
    sample_index_to_sample_info.emplace_back_several_empty_sample_infos(3, 2, &default_genotyping_options);

    std::vector<std::vector<uint32_t>> empty_coverage_vector{{0},{0}};

    EXPECT_EQ(5, sample_index_to_sample_info.size());
    for (size_t index = 0; index < sample_index_to_sample_info.size(); ++index) {
        EXPECT_EQ(index, sample_index_to_sample_info[index].get_sample_index());
        EXPECT_EQ(2, sample_index_to_sample_info[index].get_number_of_alleles());
        EXPECT_EQ(empty_coverage_vector, sample_index_to_sample_info[index].get_allele_to_forward_coverages());
        EXPECT_EQ(empty_coverage_vector, sample_index_to_sample_info[index].get_allele_to_reverse_coverages());
        EXPECT_EQ(&default_genotyping_options, sample_index_to_sample_info[index].get_genotyping_options());
    }
}

TEST_F(SampleIndexToSampleInfoTemplate___Fixture, emplace_back_several_empty_sample_infos___different_number_of_alleles___all_samples_get_fixed) {
    sample_index_to_sample_info.emplace_back_several_empty_sample_infos(2, 5, &default_genotyping_options);

    std::vector<std::vector<uint32_t>> empty_coverage_vector{{0},{0},{0},{0},{0}};

    EXPECT_EQ(4, sample_index_to_sample_info.size());
    for (size_t index = 0; index < sample_index_to_sample_info.size(); ++index) {
        EXPECT_EQ(index, sample_index_to_sample_info[index].get_sample_index());
        EXPECT_EQ(5, sample_index_to_sample_info[index].get_number_of_alleles());
        EXPECT_EQ(empty_coverage_vector, sample_index_to_sample_info[index].get_allele_to_forward_coverages());
        EXPECT_EQ(empty_coverage_vector, sample_index_to_sample_info[index].get_allele_to_reverse_coverages());
        EXPECT_EQ(&default_genotyping_options, sample_index_to_sample_info[index].get_genotyping_options());
    }
}


TEST_F(SampleIndexToSampleInfoTemplate___Fixture,
       merge_other_samples_infos_into_this___different_nb_of_samples___expects_death) {
    SampleIndexToSampleInfoTemplateAllVisible<SampleInfoMock> another_sample_index_to_sample_info;
    another_sample_index_to_sample_info.emplace_back_several_empty_sample_infos(5, 2, &default_genotyping_options);

    EXPECT_DEATH(sample_index_to_sample_info.merge_other_samples_infos_into_this(another_sample_index_to_sample_info),
                 "");
}


TEST_F(SampleIndexToSampleInfoTemplate___Fixture, merge_other_samples_infos_into_this___successful_merge) {
    SampleIndexToSampleInfoTemplateAllVisible<SampleInfoMock> another_sample_index_to_sample_info;
    another_sample_index_to_sample_info.emplace_back_several_empty_sample_infos(2, 2, &default_genotyping_options);

    EXPECT_CALL(sample_index_to_sample_info[0],
                merge_other_sample_info_into_this(Property(&SampleInfoMock::get_sample_index, 0)))
            .Times(1);
    EXPECT_CALL(sample_index_to_sample_info[1],
                merge_other_sample_info_into_this(Property(&SampleInfoMock::get_sample_index, 1)))
            .Times(1);

    sample_index_to_sample_info.merge_other_samples_infos_into_this(another_sample_index_to_sample_info);
}


TEST_F(SampleIndexToSampleInfoTemplate___Fixture, to_string___no_samples) {
    SampleIndexToSampleInfoTemplateAllVisible<SampleInfoMock> no_samples;
    std::string actual = no_samples.to_string(true, false);
    EXPECT_EQ("", actual);
}

TEST_F(SampleIndexToSampleInfoTemplate___Fixture, to_string___two_samples) {
    EXPECT_CALL(sample_index_to_sample_info[0], to_string)
            .Times(1)
            .WillOnce(Return("to_string_0"));
    EXPECT_CALL(sample_index_to_sample_info[1], to_string)
            .Times(1)
            .WillOnce(Return("to_string_1"));

    std::string actual = sample_index_to_sample_info.to_string(true, false);
    EXPECT_EQ("to_string_0\tto_string_1", actual);
}


class SampleInfoTest___solve_incompatible_gt_conflict_with___Fixture : public ::testing::Test {
public:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(bool, is_gt_from_coverages_compatible_valid, (), (const override));
        MOCK_METHOD(uint32_t, get_gt_from_coverages_compatible, (), (const override));
        MOCK_METHOD(double, get_likelihood_of_gt_from_coverages_compatible, (), (const override));
        MOCK_METHOD(void, set_gt_from_coverages_compatible, (const boost::optional<uint32_t> &gt), (override));
    };

    SampleInfoTest___solve_incompatible_gt_conflict_with___Fixture() :
            sample_info_invalid_gt(0, 2, &default_genotyping_options),
            sample_info_gt_0_with_likelihood_minus_10(0, 2, &default_genotyping_options),
            sample_info_gt_0_with_likelihood_minus_5(0, 2, &default_genotyping_options),
            sample_info_gt_3_with_likelihood_minus_10(0, 2, &default_genotyping_options),
            sample_info_gt_3_with_likelihood_minus_5(0, 2, &default_genotyping_options) {}


    void SetUp() override {
        ON_CALL(sample_info_invalid_gt, is_gt_from_coverages_compatible_valid)
        .WillByDefault(Return(false));
        ON_CALL(sample_info_gt_0_with_likelihood_minus_10, is_gt_from_coverages_compatible_valid)
                .WillByDefault(Return(true));
        ON_CALL(sample_info_gt_0_with_likelihood_minus_5, is_gt_from_coverages_compatible_valid)
                .WillByDefault(Return(true));
        ON_CALL(sample_info_gt_3_with_likelihood_minus_10, is_gt_from_coverages_compatible_valid)
                .WillByDefault(Return(true));
        ON_CALL(sample_info_gt_3_with_likelihood_minus_5, is_gt_from_coverages_compatible_valid)
                .WillByDefault(Return(true));

        ON_CALL(sample_info_gt_0_with_likelihood_minus_10, get_gt_from_coverages_compatible)
                .WillByDefault(Return(0));
        ON_CALL(sample_info_gt_0_with_likelihood_minus_5, get_gt_from_coverages_compatible)
                .WillByDefault(Return(0));
        ON_CALL(sample_info_gt_3_with_likelihood_minus_10, get_gt_from_coverages_compatible)
                .WillByDefault(Return(3));
        ON_CALL(sample_info_gt_3_with_likelihood_minus_5, get_gt_from_coverages_compatible)
                .WillByDefault(Return(3));

        ON_CALL(sample_info_gt_0_with_likelihood_minus_10, get_likelihood_of_gt_from_coverages_compatible)
                .WillByDefault(Return(-10.0));
        ON_CALL(sample_info_gt_0_with_likelihood_minus_5, get_likelihood_of_gt_from_coverages_compatible)
                .WillByDefault(Return(-5.0));
        ON_CALL(sample_info_gt_3_with_likelihood_minus_10, get_likelihood_of_gt_from_coverages_compatible)
                .WillByDefault(Return(-10.0));
        ON_CALL(sample_info_gt_3_with_likelihood_minus_5, get_likelihood_of_gt_from_coverages_compatible)
                .WillByDefault(Return(-5.0));
    }

    void TearDown() override {
    }

    SampleInfoMock sample_info_invalid_gt;
    SampleInfoMock sample_info_gt_0_with_likelihood_minus_10;
    SampleInfoMock sample_info_gt_0_with_likelihood_minus_5;
    SampleInfoMock sample_info_gt_3_with_likelihood_minus_10;
    SampleInfoMock sample_info_gt_3_with_likelihood_minus_5;
};

TEST_F(SampleInfoTest___solve_incompatible_gt_conflict_with___Fixture, first_sample_info_has_invalid_gt) {
    EXPECT_CALL(sample_info_invalid_gt, set_gt_from_coverages_compatible)
    .Times(0);
    EXPECT_CALL(sample_info_gt_0_with_likelihood_minus_10, set_gt_from_coverages_compatible)
        .Times(0);

    sample_info_invalid_gt.solve_incompatible_gt_conflict_with(
            sample_info_gt_0_with_likelihood_minus_10);
}

TEST_F(SampleInfoTest___solve_incompatible_gt_conflict_with___Fixture, second_sample_info_has_invalid_gt) {
    EXPECT_CALL(sample_info_gt_0_with_likelihood_minus_10, set_gt_from_coverages_compatible)
            .Times(0);
    EXPECT_CALL(sample_info_invalid_gt, set_gt_from_coverages_compatible)
            .Times(0);

    sample_info_gt_0_with_likelihood_minus_10.solve_incompatible_gt_conflict_with(
            sample_info_invalid_gt);
}

TEST_F(SampleInfoTest___solve_incompatible_gt_conflict_with___Fixture, both_gts_are_to_ref) {
    EXPECT_CALL(sample_info_gt_0_with_likelihood_minus_10, set_gt_from_coverages_compatible)
            .Times(0);
    EXPECT_CALL(sample_info_gt_0_with_likelihood_minus_5, set_gt_from_coverages_compatible)
            .Times(0);

    sample_info_gt_0_with_likelihood_minus_10.solve_incompatible_gt_conflict_with(
            sample_info_gt_0_with_likelihood_minus_5);
}

TEST_F(SampleInfoTest___solve_incompatible_gt_conflict_with___Fixture, first_gt_is_0_second_is_3_gt_0_higher_likelihood) {
    EXPECT_CALL(sample_info_gt_0_with_likelihood_minus_5, set_gt_from_coverages_compatible)
            .Times(0);
    EXPECT_CALL(sample_info_gt_3_with_likelihood_minus_10, set_gt_from_coverages_compatible(boost::optional<uint32_t>(0)))
            .Times(1);

    sample_info_gt_0_with_likelihood_minus_5.solve_incompatible_gt_conflict_with(
            sample_info_gt_3_with_likelihood_minus_10);
}

TEST_F(SampleInfoTest___solve_incompatible_gt_conflict_with___Fixture, first_gt_is_0_second_is_3_gt_3_higher_likelihood) {
    EXPECT_CALL(sample_info_gt_0_with_likelihood_minus_10, set_gt_from_coverages_compatible(boost::optional<uint32_t>(boost::none)))
            .Times(1);
    EXPECT_CALL(sample_info_gt_3_with_likelihood_minus_5, set_gt_from_coverages_compatible)
            .Times(0);

    sample_info_gt_0_with_likelihood_minus_10.solve_incompatible_gt_conflict_with(
            sample_info_gt_3_with_likelihood_minus_5);
}

TEST_F(SampleInfoTest___solve_incompatible_gt_conflict_with___Fixture, first_gt_is_3_second_is_0_gt_3_higher_likelihood) {
    EXPECT_CALL(sample_info_gt_3_with_likelihood_minus_5, set_gt_from_coverages_compatible)
            .Times(0);
    EXPECT_CALL(sample_info_gt_0_with_likelihood_minus_10, set_gt_from_coverages_compatible(boost::optional<uint32_t>(boost::none)))
            .Times(1);

    sample_info_gt_3_with_likelihood_minus_5.solve_incompatible_gt_conflict_with(
            sample_info_gt_0_with_likelihood_minus_10);
}

TEST_F(SampleInfoTest___solve_incompatible_gt_conflict_with___Fixture, first_gt_is_3_second_is_0_gt_0_higher_likelihood) {
    EXPECT_CALL(sample_info_gt_3_with_likelihood_minus_10, set_gt_from_coverages_compatible(boost::optional<uint32_t>(0)))
            .Times(1);
    EXPECT_CALL(sample_info_gt_0_with_likelihood_minus_5, set_gt_from_coverages_compatible)
            .Times(0);

    sample_info_gt_3_with_likelihood_minus_10.solve_incompatible_gt_conflict_with(
            sample_info_gt_0_with_likelihood_minus_5);
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
// LIKELIHOOD TESTS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// REASON COMMENTED OUT: THERE IS NO WAY TO TRY TO COMPUTE THE LIKELIHOOD WITHOUT HAVING A SAMPLE ANYMORE
//TEST(VCFRecordLikelihoodTest, does_not_crash_with_no_samples) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    EXPECT_NO_FATAL_FAILURE(vr.likelihood({}, 0.01, 0));
//}
//

// REASON COMMENTED OUT: INFO IS NOW FIXED: LIKELIHOOD IS ALWAYS IN INFO
//TEST(VCFRecordLikelihoodTest, does_not_run_if_info_missing) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    SampleInfo<uint16_t> sample_info;
//    sample_info["nothing"] = {0};
//    vr.sampleIndex_to_sampleInfo.push_back(sample_info);
//    assert(vr.sampleIndex_to_sampleInfo.size() > 0);
//    std::vector<float> f = {0.0, 0.0};
//    vr.set_format(0,"GAPS", f);
//    vr.likelihood({1}, 0.01, 0);
//    bool found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
//    EXPECT_FALSE(found_likelihood);
//
//    vr.sampleIndex_to_sampleInfo[0]["GT"] = {1};
//    vr.likelihood({1}, 0.01, 0);
//    found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
//    EXPECT_FALSE(found_likelihood);
//
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 1};
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {1};
//    vr.likelihood({1}, 0.01, 0);
//    found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
//    EXPECT_FALSE(found_likelihood);
//
//    vr.sampleIndex_to_sampleInfo[0].erase("MEAN_FWD_COVG");
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 1};
//    vr.likelihood({1}, 0.01, 0);
//    found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
//    EXPECT_FALSE(found_likelihood);
//
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1};
//    vr.likelihood({1}, 0.01, 0);
//    found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
//    EXPECT_FALSE(found_likelihood);
//
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 1};
//    vr.sampleIndex_to_sampleInfo[0].erase("MEAN_REV_COVG");
//    vr.likelihood({1}, 0.01, 0);
//    found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
//    EXPECT_FALSE(found_likelihood);
//}
//

// REASON COMMENTED OUT: IT IS NOT POSSIBLE TO ADD CUSTOM FIELDS TO INFO ANYMORE
//TEST(VCFRecordLikelihoodTest, adds_likelihood_with_info) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vr.sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1);
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 2};
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 2};
//    std::vector<float> f = {0.0, 0.0};
//    vr.set_format(0,"GAPS", f);
//    vr.likelihood({1}, 0.01, 0);
//    bool found_likelihood = !vr.get_format_f(0, "LIKELIHOOD").empty();
//    EXPECT_TRUE(found_likelihood);
//}
//
TEST_F(SampleInfoTest___Fixture, get_likelihoods_for_all_alleles___gets_correct_likelihood_simple_case) {
    default_sample_info.set_coverage_information({{1},
                                                  {2}}, {{1},
                                                         {2}});

    double actual = default_sample_info.get_likelihoods_for_all_alleles()[0];
    double expected = -1.0 - log(2.0) + 4.0 * log(0.01) + log(1 - exp(-1.0));
    EXPECT_NEAR(actual, expected, 0.00001);

    actual = default_sample_info.get_likelihoods_for_all_alleles()[1];
    expected = -1 - log(4) - log(3) - log(2) + 2 * log(0.01) + log(1 - exp(-(float(1))));
    EXPECT_NEAR(actual, expected, 0.00001);
}

TEST(SampleInfoTest, get_likelihoods_for_all_alleles___gets_correct_likelihood_with_min_covg_threshold) {
    GenotypingOptions genotyping_options({1}, 0.01, 0, 3, 0, 0, 0, 0, false);
    SampleInfo sample_info(0, 2, &genotyping_options);
    sample_info.set_coverage_information({{1},
                                          {2}}, {{1},
                                                 {2}});


    double actual = sample_info.get_likelihoods_for_all_alleles()[0];
    double expected = 4 * log(0.01) - 1 + log(1 - exp(-(float(1))));
    EXPECT_NEAR(actual, expected, 0.00001);

    actual = sample_info.get_likelihoods_for_all_alleles()[1];
    expected = -1 - log(4) - log(3) - log(2) + log(1 - exp(-(float(1))));
    EXPECT_NEAR(actual, expected, 0.00001);
}


TEST_F(SampleInfoTest___Fixture, get_likelihoods_for_all_alleles___handles_ref_covg_0) {
    default_sample_info.set_coverage_information({{0},
                                                  {2}}, {{0},
                                                         {2}});

    double actual = default_sample_info.get_likelihoods_for_all_alleles()[0];
    double expected = -1 + 4 * log(0.01) + log(1 - exp(-(float(1))));
    EXPECT_NEAR(actual, expected, 0.00001);

    actual = default_sample_info.get_likelihoods_for_all_alleles()[1];
    expected = -1 - log(4) - log(3) - log(2) + log(1 - exp(-(float(1))));
    EXPECT_NEAR(actual, expected, 0.00001);
}

TEST_F(SampleInfoTest___Fixture, get_likelihoods_for_all_alleles___handles_alt_covg_0) {
    default_sample_info.set_coverage_information({{1},
                                                  {0}}, {{1},
                                                         {0}});

    double actual = default_sample_info.get_likelihoods_for_all_alleles()[0];
    double expected = -1 - log(2) + log(1 - exp(-(float(1))));
    EXPECT_NEAR(actual, expected, 0.00001);

    actual = default_sample_info.get_likelihoods_for_all_alleles()[1];
    expected = -1 + 2 * log(0.01) + log(1 - exp(-(float(1))));
    EXPECT_NEAR(actual, expected, 0.00001);
}


class SampleInfoTest___gets_correct_likelihood_gaps___Fixture : public ::testing::Test {
public:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(double, get_gaps, (uint32_t
                allele), (const override));
    };

    SampleInfoTest___gets_correct_likelihood_gaps___Fixture() :
            sample_info(0, 2, &default_genotyping_options) {}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfoMock sample_info;
};

TEST_F(SampleInfoTest___gets_correct_likelihood_gaps___Fixture,
       get_likelihoods_for_all_alleles___gets_correct_likelihood_gaps) {
    EXPECT_CALL(sample_info, get_gaps(0))
            .WillRepeatedly(Return(0.5));

    EXPECT_CALL(sample_info, get_gaps(1))
            .WillRepeatedly(Return(0.8));

    sample_info.set_coverage_information({{1},
                                          {2}}, {{1},
                                                 {2}});

    double actual = sample_info.get_likelihoods_for_all_alleles()[0];
    double expected = -1 - log(2) + 4 * log(0.01) + 0.5 * log(1 - exp(-(float(1)))) - 0.5;
    EXPECT_NEAR(actual, expected, 0.00001);

    actual = sample_info.get_likelihoods_for_all_alleles()[1];
    expected = -1 - log(4) - log(3) - log(2) + 2 * log(0.01) + 0.2 * log(1 - exp(-(float(1)))) - 0.8;
    EXPECT_NEAR(actual, expected, 0.00001);
}

// REASON COMMENTED OUT: ALREADY TESTED, SEE TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___not_enough_total_covg)
//TEST(VCFRecordLikelihoodTest, death_not_enough_covgs) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vr.sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(2);
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 2};
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 2};
//    vr.sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"] = {1, 2};
//    vr.sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"] = {1, 2};
//    std::vector<float> f = {0.5, 0.8};
//    vr.set_format(0,"GAPS", f);
//    vr.set_format(1,"GAPS", f);
//    EXPECT_DEATH(vr.likelihood({1}, 0.01, 0), "");
//}
//
// REASON COMMENTED OUT: THERE IS NO REASON FOR THIS TEST AS DIFFERENT DEPTHS IN DIFFERENT SAMPLES DO NOT AFFECT LIKELIHOOD COMPUTATION
//TEST(VCFRecordLikelihoodTest, samples_with_different_depths) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vr.sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(2);
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {1, 2};
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 2};
//    vr.sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"] = {1, 2};
//    vr.sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"] = {1, 2};
//    std::vector<float> f = {0.5, 0.8};
//    vr.set_format(0,"GAPS", f);
//    vr.set_format(1,"GAPS", f);
//    vr.likelihood({1,2}, 0.01, 0);
//
//    float exp_likelihood = -1 - log(2) + 4 * log(0.01) + 0.5*log(1-exp(-(float(1)))) - 0.5;
//    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_sampleInfo[0]["LIKELIHOOD"][0]);
//    exp_likelihood = -1 - log(4) - log(3) - log(2) + 2 * log(0.01) + 0.2*log(1-exp(-(float(1)))) - 0.8;
//    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_sampleInfo[0]["LIKELIHOOD"][1]);
//    exp_likelihood = 2*log(2) -2 - log(2) + 4 * log(0.01) + 0.5*log(1-exp(-(float(2)))) - 2*0.5;
//    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_sampleInfo[1]["LIKELIHOOD"][0]);
//    exp_likelihood = 4*log(2) -2 - log(4) - log(3) - log(2) + 2 * log(0.01) + 0.2*log(1-exp(-(float(2)))) - 2*0.8;
//    EXPECT_FLOAT_EQ(exp_likelihood, vr.sampleIndex_to_sampleInfo[1]["LIKELIHOOD"][1]);
//}
//
//

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CONFIDENCE TESTS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// REASON COMMENTED OUT: THERE IS NO WAY TO TRY TO COMPUTE THE CONFIDENCE WITHOUT HAVING INFO ANYMORE
//TEST(VCFRecordConfidenceTest, does_not_run_if_info_missing) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vr.sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1);
//    vr.confidence();
//    bool found_confidence = !vr.get_format_f(0,"GT_CONF").empty();
//    EXPECT_FALSE(found_confidence);
//    std::vector<float> f = {-1.0};
//    vr.set_format(0,"LIKELIHOOD", f);
//    EXPECT_DEATH(vr.confidence(), "");
//}
//
// REASON COMMENTED OUT: INFO IS NOW FIXED: CONFIDENCE IS ALWAYS IN INFO
//TEST(VCFRecordConfidenceTest, adds_confidence_with_info) {
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vr.sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1);
//    vr.sampleIndex_to_sampleInfo[0]["LIKELIHOOD"] = {-1.0, -2.5};
//    vr.sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(1);
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 0};
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {0, 0};
//    vr.confidence();
//    bool found_confidence = !vr.get_format_f(0,"GT_CONF").empty();
//    EXPECT_TRUE(found_confidence);
//}



TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___gets_correct_confidence_simple_case) {
    EXPECT_CALL(default_sample_info, get_likelihoods_for_all_alleles)
    .WillRepeatedly(Return(std::vector<double>({-1.0, 0.0})));
    EXPECT_CALL(default_sample_info, get_mean_coverage_both_alleles)
    .WillRepeatedly(Return(0));

    size_t index;
    double confidence, max_likelihood;
    std::tie(index, confidence, max_likelihood) = *(default_sample_info.get_confidence());
    EXPECT_EQ(1, index);
    EXPECT_NEAR(1.0, confidence, 0.000001);
    EXPECT_NEAR(0.0, max_likelihood, 0.000001);
}


TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___gets_correct_confidence_two_alts) {
    EXPECT_CALL(default_sample_info, get_likelihoods_for_all_alleles)
            .WillRepeatedly(Return(std::vector<double>({-14.0, -6.0, -3.0})));
    EXPECT_CALL(default_sample_info, get_mean_coverage_both_alleles)
            .WillRepeatedly(Return(0));

    size_t index;
    double confidence, max_likelihood;
    std::tie(index, confidence, max_likelihood) = *(default_sample_info.get_confidence());
    EXPECT_EQ(2, index);
    EXPECT_NEAR(3.0, confidence, 0.000001);
    EXPECT_NEAR(-3.0, max_likelihood, 0.000001);
}


TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___gets_correct_confidence_min_total___confidence_is_invalid) {
    GenotypingOptions genotyping_options({10, 10}, 0.01, 0, 0, 0.0, 3, 0, 0, 0);
    SampleInfoTest___get_confidence___Fixture::SampleInfoMock sample_info(0, 2, &genotyping_options);
    EXPECT_CALL(sample_info, get_likelihoods_for_all_alleles)
            .WillRepeatedly(Return(std::vector<double>({-14.0, -6.0, -3.0})));
    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(_))
            .WillRepeatedly(Return(0));
    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(2))
            .WillRepeatedly(Return(2));

    auto invalid_confidence = sample_info.get_confidence();
    EXPECT_EQ(boost::none, invalid_confidence);
}

TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___gets_correct_confidence_min_total___confidence_is_valid) {
    GenotypingOptions genotyping_options({10, 10}, 0.01, 0, 0, 0.0, 2, 0, 0, 0);
    SampleInfoTest___get_confidence___Fixture::SampleInfoMock sample_info(0, 2, &genotyping_options);
    EXPECT_CALL(sample_info, get_likelihoods_for_all_alleles)
            .WillRepeatedly(Return(std::vector<double>({-14.0, -6.0, -3.0})));
    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(_))
            .WillRepeatedly(Return(0));
    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(2))
            .WillRepeatedly(Return(2));

    size_t index;
    double confidence, max_likelihood;
    std::tie(index, confidence, max_likelihood) = *(sample_info.get_confidence());
    EXPECT_EQ(2, index);
    EXPECT_NEAR(3.0, confidence, 0.000001);
    EXPECT_NEAR(-3.0, max_likelihood, 0.000001);
}


TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___gets_correct_confidence_min_diff) {
    GenotypingOptions genotyping_options({10, 10}, 0.01, 0, 0, 0.0, 0, 3, 0, 0);
    SampleInfoTest___get_confidence___Fixture::SampleInfoMock sample_info(0, 2, &genotyping_options);
    EXPECT_CALL(sample_info, get_likelihoods_for_all_alleles)
            .WillRepeatedly(Return(std::vector<double>({-14.0, -6.0, -3.0})));
    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(0))
            .WillRepeatedly(Return(0));
    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(1))
            .WillRepeatedly(Return(2));
    EXPECT_CALL(sample_info, get_mean_coverage_both_alleles(2))
            .WillRepeatedly(Return(5));

    size_t index;
    double confidence, max_likelihood;
    std::tie(index, confidence, max_likelihood) = *(sample_info.get_confidence());
    EXPECT_EQ(2, index);
    EXPECT_NEAR(3.0, confidence, 0.000001);
    EXPECT_NEAR(-3.0, max_likelihood, 0.000001);
}


TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___handles_ref_covg_0) {
    EXPECT_CALL(default_sample_info, get_likelihoods_for_all_alleles)
            .WillRepeatedly(Return(std::vector<double>({std::numeric_limits<double>::lowest(), -1.5})));
    EXPECT_CALL(default_sample_info, get_mean_coverage_both_alleles)
            .WillRepeatedly(Return(0));

    size_t index;
    double confidence, max_likelihood;
    std::tie(index, confidence, max_likelihood) = *(default_sample_info.get_confidence());
    EXPECT_EQ(1, index);
    EXPECT_NEAR(-std::numeric_limits<double>::lowest() - 1.5, confidence, 0.000001);
    EXPECT_NEAR(-1.5, max_likelihood, 0.000001);
}


TEST_F(SampleInfoTest___get_confidence___Fixture, get_confidence___handles_alt_covg_0) {
    EXPECT_CALL(default_sample_info, get_likelihoods_for_all_alleles)
            .WillRepeatedly(Return(std::vector<double>({-1.5, std::numeric_limits<double>::lowest()})));
    EXPECT_CALL(default_sample_info, get_mean_coverage_both_alleles)
            .WillRepeatedly(Return(0));

    size_t index;
    double confidence, max_likelihood;
    std::tie(index, confidence, max_likelihood) = *(default_sample_info.get_confidence());
    EXPECT_EQ(0, index);
    EXPECT_NEAR(-std::numeric_limits<double>::lowest() - 1.5, confidence, 0.000001);
    EXPECT_NEAR(-1.5, max_likelihood, 0.000001);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// REGENOTYPE TESTS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// REASON COMMENTED OUT: ALL TESTS ARE EITHER REDUNDANT OR NOT APPLICABLE ANYMORE
// sample 0 missing confidence: done in TEST_F(SampleInfoTest___get_genotype_from_coverage___Fixture, invalid_confidence)
// sample 1 confidence below threshold: done in TEST_F(SampleInfoTest___get_genotype_from_coverage___Fixture, valid_confidence_but_below_threshold)
// sample 2 confidence above threshold, but has correct GT 0 already: done in TEST_F(SampleInfoTest___get_genotype_from_coverage___Fixture, valid_confidence_and_above_threshold)
// sample 3 confidence above threshold, but has correct GT 1 already: done in TEST_F(SampleInfoTest___get_genotype_from_coverage___Fixture, valid_confidence_and_above_threshold)
// sample 4 confidence above threshold, has incorrect GT 0: there is no way to have an incorrect GT from coverage: it is derived from the coverage
// sample 5 confidence above threshold, has incorrect GT 1: there is no way to have an incorrect GT from coverage: it is derived from the coverage
//TEST(VCFRecordRegenotypeTest, correctly_genotypes) {
//    // sample 0 missing confidence
//    // sample 1 confidence below threshold
//    // sample 2 confidence above threshold, but has correct GT 0 already
//    // sample 3 confidence above threshold, but has correct GT 1 already
//    // sample 4 confidence above threshold, has incorrect GT 0
//    // sample 5 confidence above threshold, has incorrect GT 1
//
//    VCFRecord vr("chrom1", 3, "A", "T");
//    vr.sampleIndex_to_sampleInfo.emplace_back_several_empty_sample_infos(6);
//
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"] = {0, 2};
//    vr.sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"] = {1, 3};
//    vr.sampleIndex_to_sampleInfo[0]["LIKELIHOOD"] = {4, 5};
//    vr.sampleIndex_to_sampleInfo[0]["GT"] = {1};
//
//    vr.sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"] = {0, 2};
//    vr.sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"] = {1, 3};
//    vr.sampleIndex_to_sampleInfo[1]["LIKELIHOOD"] = {4, 5};
//    vr.sampleIndex_to_sampleInfo[1]["GT"] = {1};
//    vr.sampleIndex_to_sampleInfo[1]["GT_CONF"] = {1};
//
//    vr.sampleIndex_to_sampleInfo[2]["MEAN_FWD_COVG"] = {0, 2};
//    vr.sampleIndex_to_sampleInfo[2]["MEAN_REV_COVG"] = {1, 3};
//    vr.sampleIndex_to_sampleInfo[2]["LIKELIHOOD"] = {6, 4};
//    vr.sampleIndex_to_sampleInfo[2]["GT"] = {0};
//    vr.sampleIndex_to_sampleInfo[2]["GT_CONF"] = {2};
//
//    vr.sampleIndex_to_sampleInfo[3]["MEAN_FWD_COVG"] = {0, 2};
//    vr.sampleIndex_to_sampleInfo[3]["MEAN_REV_COVG"] = {1, 3};
//    vr.sampleIndex_to_sampleInfo[3]["LIKELIHOOD"] = {4, 6};
//    vr.sampleIndex_to_sampleInfo[3]["GT"] = {1};
//    vr.sampleIndex_to_sampleInfo[3]["GT_CONF"] = {2};
//
//    vr.sampleIndex_to_sampleInfo[4]["MEAN_FWD_COVG"] = {0, 2};
//    vr.sampleIndex_to_sampleInfo[4]["MEAN_REV_COVG"] = {1, 3};
//    vr.sampleIndex_to_sampleInfo[4]["LIKELIHOOD"] = {6, 4};
//    vr.sampleIndex_to_sampleInfo[4]["GT"] = {1};
//    vr.sampleIndex_to_sampleInfo[4]["GT_CONF"] = {2};
//
//    vr.sampleIndex_to_sampleInfo[5]["MEAN_FWD_COVG"] = {0, 2};
//    vr.sampleIndex_to_sampleInfo[5]["MEAN_REV_COVG"] = {1, 3};
//    vr.sampleIndex_to_sampleInfo[5]["LIKELIHOOD"] = {4, 6};
//    vr.sampleIndex_to_sampleInfo[5]["GT"] = {0};
//    vr.sampleIndex_to_sampleInfo[5]["GT_CONF"] = {2};
//
//    vr.genotype(1);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"][0], 0);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"][0], 1);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[0]["MEAN_FWD_COVG"][1], 2);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[0]["MEAN_REV_COVG"][1], 3);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[0]["LIKELIHOOD"][0], 4);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[0]["LIKELIHOOD"][1], 5);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[0]["GT"].size(), (uint) 0);
//
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"][0], 0);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"][0], 1);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[1]["MEAN_FWD_COVG"][1], 2);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[1]["MEAN_REV_COVG"][1], 3);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[1]["LIKELIHOOD"][0], 4);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[1]["LIKELIHOOD"][1], 5);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[1]["GT"].size(), (uint) 0);
//
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[2]["MEAN_FWD_COVG"][0], 0);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[2]["MEAN_REV_COVG"][0], 1);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[2]["MEAN_FWD_COVG"][1], 2);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[2]["MEAN_REV_COVG"][1], 3);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[2]["LIKELIHOOD"][0], 6);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[2]["LIKELIHOOD"][1], 4);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[2]["GT"][0], 0);
//
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[3]["MEAN_FWD_COVG"][0], 0);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[3]["MEAN_REV_COVG"][0], 1);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[3]["MEAN_FWD_COVG"][1], 2);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[3]["MEAN_REV_COVG"][1], 3);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[3]["LIKELIHOOD"][0], 4);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[3]["LIKELIHOOD"][1], 6);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[3]["GT"][0], 1);
//
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[4]["MEAN_FWD_COVG"][0], 0);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[4]["MEAN_REV_COVG"][0], 1);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[4]["MEAN_FWD_COVG"][1], 2);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[4]["MEAN_REV_COVG"][1], 3);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[4]["LIKELIHOOD"][0], 6);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[4]["LIKELIHOOD"][1], 4);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[4]["GT"][0], 0);
//
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[5]["MEAN_FWD_COVG"][0], 0);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[5]["MEAN_REV_COVG"][0], 1);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[5]["MEAN_FWD_COVG"][1], 2);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[5]["MEAN_REV_COVG"][1], 3);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[5]["LIKELIHOOD"][0], 4);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[5]["LIKELIHOOD"][1], 6);
//    EXPECT_EQ(vr.sampleIndex_to_sampleInfo[5]["GT"][0], 1);
//}
