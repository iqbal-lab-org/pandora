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

class SampleInfoTest___Fixture : public ::testing::Test {
public:
    SampleInfoTest___Fixture() :
            default_sample_info(0, &default_genotyping_options){}

    void SetUp() override {
        allele_to_coverage_one_allele.push_back({1, 2});

        allele_to_coverage_two_alleles.push_back({1, 2});
        allele_to_coverage_two_alleles.push_back({3, 4});

        allele_to_coverage_three_alleles.push_back({1, 2});
        allele_to_coverage_three_alleles.push_back({3, 4});
        allele_to_coverage_three_alleles.push_back({5, 6});
    }

    void TearDown() override {
    }

    SampleInfo default_sample_info;
    std::vector< std::vector<uint32_t> > allele_to_coverage_empty;
    std::vector< std::vector<uint32_t> > allele_to_coverage_one_allele;
    std::vector< std::vector<uint32_t> > allele_to_coverage_two_alleles;
    std::vector< std::vector<uint32_t> > allele_to_coverage_three_alleles;
};

TEST_F(SampleInfoTest___Fixture, get_and_set_gt_from_max_likelihood_path___set_a_valid_gt) {
    default_sample_info.set_gt_from_max_likelihood_path(5);

    uint32_t actual = default_sample_info.get_gt_from_max_likelihood_path();
    uint32_t expected = 5;

    EXPECT_EQ(actual, expected);
}

TEST_F(SampleInfoTest___Fixture, get_and_set_gt_from_max_likelihood_path___set_a_valid_then_an_invalid_gt___throws_exception) {
    default_sample_info.set_gt_from_max_likelihood_path(5);
    default_sample_info.set_gt_from_max_likelihood_path(boost::none);

    EXPECT_FALSE(default_sample_info.is_gt_from_max_likelihood_path_valid());
    try {
        default_sample_info.get_gt_from_max_likelihood_path();
        FAIL() << "Should have thrown exception";
    }catch (std::runtime_error &error) {}
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
    }catch (std::runtime_error &error) {}
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
    std::vector< std::vector<uint32_t> > actual = default_sample_info.get_allele_to_forward_coverages();
    std::vector< std::vector<uint32_t> > expected;

    EXPECT_EQ(actual, expected);
}

TEST_F(SampleInfoTest___Fixture, get_allele_to_reverse_coverages___default_sample_info) {
    std::vector< std::vector<uint32_t> > actual = default_sample_info.get_allele_to_reverse_coverages();
    std::vector< std::vector<uint32_t> > expected;

    EXPECT_EQ(actual, expected);
}

TEST_F(SampleInfoTest___Fixture, add_coverage_information___forward_coverage_has_no_alleles___expects_death) {
    EXPECT_DEATH(default_sample_info.add_coverage_information(allele_to_coverage_empty, allele_to_coverage_three_alleles), "");
}

TEST_F(SampleInfoTest___Fixture, add_coverage_information___forward_coverage_has_one_allele___expects_death) {
    EXPECT_DEATH(default_sample_info.add_coverage_information(allele_to_coverage_one_allele, allele_to_coverage_three_alleles), "");
}

TEST_F(SampleInfoTest___Fixture, add_coverage_information___reverse_coverage_has_no_alleles___expects_death) {
    EXPECT_DEATH(default_sample_info.add_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_empty), "");
}

TEST_F(SampleInfoTest___Fixture, add_coverage_information___reverse_coverage_has_one_allele___expects_death) {
    EXPECT_DEATH(default_sample_info.add_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_one_allele), "");
}

TEST_F(SampleInfoTest___Fixture, add_coverage_information___both_coverages_have_two_alleles) {
    default_sample_info.add_coverage_information(allele_to_coverage_two_alleles, allele_to_coverage_two_alleles);

    auto expected = allele_to_coverage_two_alleles;

    auto actual = default_sample_info.get_allele_to_forward_coverages();
    EXPECT_EQ(actual, expected);

    actual = default_sample_info.get_allele_to_reverse_coverages();
    EXPECT_EQ(actual, expected);
}


TEST_F(SampleInfoTest___Fixture, add_coverage_information___forward_covg_has_two_alleles___reverse_covg_has_three_alleles___expects_death) {
    EXPECT_DEATH(default_sample_info.add_coverage_information(allele_to_coverage_two_alleles, allele_to_coverage_three_alleles), "");
}

TEST_F(SampleInfoTest___Fixture, gt_coverages_compatible___default_sample_info___invalid_gt) {
    EXPECT_FALSE(default_sample_info.is_gt_from_coverages_valid());
    try {
        default_sample_info.get_gt_from_coverages();
        FAIL() << "Should have thrown exception";
    }catch (const std::runtime_error &e){ }
    EXPECT_EQ(".", default_sample_info.gt_from_coverages_compatible_to_string());
}


TEST_F(SampleInfoTest___Fixture, gt_coverages_compatible___default_sample_info___valid_gt) {
    default_sample_info.set_gt_coverages_compatible(5);
    EXPECT_TRUE(default_sample_info.is_gt_from_coverages_compatible_valid());
    EXPECT_EQ(5, default_sample_info.get_gt_coverages_compatible());
    EXPECT_EQ("5", default_sample_info.gt_from_coverages_compatible_to_string());
}

class SampleInfoTest___genotype_from_coverage___Fixture : public ::testing::Test {
private:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(void, check_if_coverage_information_is_correct, (), (const override));
        MOCK_METHOD(boost::optional<GenotypeAndMaxLikelihood>, get_genotype_from_coverage, (), (const override));
    };

public:
    SampleInfoTest___genotype_from_coverage___Fixture() :
            default_sample_info(0, &default_genotyping_options){}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfoMock default_sample_info;
};

TEST_F(SampleInfoTest___genotype_from_coverage___Fixture, valid_genotype) {
    EXPECT_CALL(default_sample_info, check_if_coverage_information_is_correct)
    .Times(1);

    EXPECT_CALL(default_sample_info, get_genotype_from_coverage).
    Times(1).
    WillOnce(Return(std::make_pair<uint32_t, double>(0, -1.0)));

    default_sample_info.genotype_from_coverage();

    EXPECT_TRUE(default_sample_info.is_gt_from_coverages_valid());
    EXPECT_EQ(0, default_sample_info.get_gt_from_coverages());
    EXPECT_EQ(-1.0, default_sample_info.get_likelihood_of_gt_from_coverages());
}

TEST_F(SampleInfoTest___genotype_from_coverage___Fixture, invalid_genotype) {
    EXPECT_FALSE(default_sample_info.is_gt_from_coverages_valid());
    try {
        default_sample_info.get_gt_from_coverages();
        FAIL() << "Should have thrown exception";
    }catch(const std::runtime_error &e){}
    try {
        default_sample_info.get_likelihood_of_gt_from_coverages();
        FAIL() << "Should have thrown exception";
    }catch(const std::runtime_error &e){}
}



class SampleInfoTest___merge_other_sample_info_into_this___Fixture : public ::testing::Test {
private:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(void, genotype_from_coverage, (), ());
        MOCK_METHOD(uint32_t, get_gt_from_coverages, (), (const override));
    };


public:
    SampleInfoTest___merge_other_sample_info_into_this___Fixture() :
            sample_info_with_two_alleles(0, &default_genotyping_options),
            sample_info_with_three_alleles(0, &default_genotyping_options)
            {}

    void SetUp() override {
        allele_to_coverage_two_alleles.push_back({1, 2});
        allele_to_coverage_two_alleles.push_back({3, 4});
        sample_info_with_two_alleles.add_coverage_information(allele_to_coverage_two_alleles, allele_to_coverage_two_alleles);

        allele_to_coverage_three_alleles.push_back({1, 2});
        allele_to_coverage_three_alleles.push_back({5, 6});
        allele_to_coverage_three_alleles.push_back({7, 8});
        sample_info_with_three_alleles.add_coverage_information(allele_to_coverage_three_alleles, allele_to_coverage_three_alleles);
    }

    void TearDown() override {
    }

    SampleInfoMock sample_info_with_two_alleles;
    SampleInfoMock sample_info_with_three_alleles;
    std::vector< std::vector<uint32_t> > allele_to_coverage_two_alleles;
    std::vector< std::vector<uint32_t> > allele_to_coverage_three_alleles;
};

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture, merge_a_sample_with_three_alleles_into_one_with_two_alleles___same_gt_from_max_likelihood_path) {
    sample_info_with_two_alleles.set_gt_from_max_likelihood_path(0);
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(0);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    std::vector< std::vector<uint32_t> > allele_to_coverage_expected({{1, 2}, {3,4}, {5,6}, {7,8}});
    EXPECT_EQ(allele_to_coverage_expected, sample_info_with_two_alleles.get_allele_to_forward_coverages());
    EXPECT_EQ(allele_to_coverage_expected, sample_info_with_two_alleles.get_allele_to_reverse_coverages());
    EXPECT_EQ(0, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture, merge_a_sample_with_three_alleles_into_one_with_two_alleles___both_gt_from_max_likelihood_path_are_invalid) {
    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_FALSE(sample_info_with_two_alleles.is_gt_from_max_likelihood_path_valid());
}

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture, merge_a_sample_with_three_alleles_into_one_with_two_alleles___gt_merged_in_is_invalid) {
    sample_info_with_two_alleles.set_gt_from_max_likelihood_path(1);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(1, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}


TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture, merge_a_sample_with_three_alleles_into_one_with_two_alleles___original_gt_invalid___gt_merged_in_is_valid_zero) {
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(0);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(0, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture, merge_a_sample_with_three_alleles_into_one_with_two_alleles___original_gt_invalid___gt_merged_in_is_valid_one) {
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(1);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(2, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture, merge_a_sample_with_three_alleles_into_one_with_two_alleles___original_gt_invalid___gt_merged_in_is_valid_two) {
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(2);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(3, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}

TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture, merge_a_sample_with_three_alleles_into_one_with_two_alleles___original_gt_invalid___both_gts_are_valid_and_first_is_not_zero___genotypes_from_coverage_to_solve_conflict) {
    EXPECT_CALL(sample_info_with_two_alleles, get_gt_from_coverages()).
            Times(1).
            WillOnce(Return(1));
    sample_info_with_two_alleles.set_gt_from_max_likelihood_path(1);
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(0);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(1, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}


TEST_F(SampleInfoTest___merge_other_sample_info_into_this___Fixture, merge_a_sample_with_three_alleles_into_one_with_two_alleles___original_gt_invalid___both_gts_are_valid_and_second_is_not_zero___genotypes_from_coverage_to_solve_conflict) {
    EXPECT_CALL(sample_info_with_two_alleles, get_gt_from_coverages()).
    Times(1).
    WillOnce(Return(3));
    sample_info_with_two_alleles.set_gt_from_max_likelihood_path(0);
    sample_info_with_three_alleles.set_gt_from_max_likelihood_path(2);

    sample_info_with_two_alleles.merge_other_sample_info_into_this(sample_info_with_three_alleles);

    EXPECT_EQ(3, sample_info_with_two_alleles.get_gt_from_max_likelihood_path());
}


class SampleInfoTest___get_gaps___Fixture : public ::testing::Test {
public:
    static GenotypingOptions genotyping_options_with_min_kmer_covg_10;
    SampleInfoTest___get_gaps___Fixture() :
            sample_info_with_min_kmer_covg_10(0, &genotyping_options_with_min_kmer_covg_10){}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfo sample_info_with_min_kmer_covg_10;
};

GenotypingOptions SampleInfoTest___get_gaps___Fixture::genotyping_options_with_min_kmer_covg_10(
        {10,10,10,10,10,10,10,10,10,10}, 0.01, 0, 0, 0, 0, 0, 10, false);

TEST_F(SampleInfoTest___get_gaps___Fixture, get_gaps___coverages_just_below_threshold) {
    sample_info_with_min_kmer_covg_10.add_coverage_information({{5}, {5,5,5,5}}, {{4}, {4,4,4,4}});

    double actual_allele_0 = sample_info_with_min_kmer_covg_10.get_gaps(0);
    double actual_allele_1 = sample_info_with_min_kmer_covg_10.get_gaps(1);

    EXPECT_TRUE(Maths::equals(1.0, actual_allele_0));
    EXPECT_TRUE(Maths::equals(1.0, actual_allele_1));
}

TEST_F(SampleInfoTest___get_gaps___Fixture, get_gaps___coverages_all_equal_threshold) {
    sample_info_with_min_kmer_covg_10.add_coverage_information({{5}, {5,5,5,5}}, {{5}, {5,5,5,5}});

    double actual_allele_0 = sample_info_with_min_kmer_covg_10.get_gaps(0);
    double actual_allele_1 = sample_info_with_min_kmer_covg_10.get_gaps(1);

    EXPECT_TRUE(Maths::equals(0.0, actual_allele_0));
    EXPECT_TRUE(Maths::equals(0.0, actual_allele_1));
}

TEST_F(SampleInfoTest___get_gaps___Fixture, get_gaps___coverages_all_above_threshold) {
    sample_info_with_min_kmer_covg_10.add_coverage_information({{5}, {5,5,5,5}}, {{6}, {6,6,6,6}});

    double actual_allele_0 = sample_info_with_min_kmer_covg_10.get_gaps(0);
    double actual_allele_1 = sample_info_with_min_kmer_covg_10.get_gaps(1);

    EXPECT_TRUE(Maths::equals(0.0, actual_allele_0));
    EXPECT_TRUE(Maths::equals(0.0, actual_allele_1));
}

TEST_F(SampleInfoTest___get_gaps___Fixture, get_gaps___coverages_below_equal_and_above_threshold) {
    sample_info_with_min_kmer_covg_10.add_coverage_information({{0,0,0}, {9,10,11,9,10,9}}, {{9,10,11}, {0,0,0,0,0,0}});

    double actual_allele_0 = sample_info_with_min_kmer_covg_10.get_gaps(0);
    double actual_allele_1 = sample_info_with_min_kmer_covg_10.get_gaps(1);

    EXPECT_TRUE(Maths::equals(1.0/3.0, actual_allele_0));
    EXPECT_TRUE(Maths::equals(0.5, actual_allele_1));
}


TEST(SampleInfoTest, get_min_coverage_threshold_for_this_sample___min_allele_covg_is_higher) {
    GenotypingOptions genotyping_options(
            {10,5}, 0.01, 0, 100, 1.0, 0, 0, 00, false);
    SampleInfo sample_info(0, &genotyping_options);

    EXPECT_EQ(100, sample_info.get_min_coverage_threshold_for_this_sample());
}

TEST(SampleInfoTest, get_min_coverage_threshold_for_this_sample___min_fraction_allele_covg_is_higher) {
    GenotypingOptions genotyping_options(
            {10,100}, 0.01, 0, 40, 0.5, 0, 0, 00, false);
    SampleInfo sample_info(1, &genotyping_options);

    EXPECT_EQ(50, sample_info.get_min_coverage_threshold_for_this_sample());
}


class SampleInfoTest___get_likelihoods_for_all_alleles___Fixture : public ::testing::Test {
public:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;

        MOCK_METHOD(size_t, get_number_of_alleles, (), (const));
        MOCK_METHOD(uint32_t, get_min_coverage_threshold_for_this_sample, (), (const));
        MOCK_METHOD(uint32_t, get_total_mean_coverage_over_all_alleles_given_a_mininum_threshold, (uint32_t minimum_threshold), (const));
        MOCK_METHOD(uint32_t, get_total_mean_coverage_given_a_mininum_threshold, (uint32_t allele, uint32_t minimum_threshold), (const));
        MOCK_METHOD(double, get_gaps, (uint32_t allele), (const));
        MOCK_METHOD(double, compute_likelihood, (bool min_coverage_threshold_is_satisfied, uint32_t expected_depth_covg, uint32_t total_mean_coverage_of_allele_above_threshold,
                uint32_t total_mean_coverage_of_all_other_alleles_above_threshold, double error_rate, double gaps), (const));
    };

    SampleInfoTest___get_likelihoods_for_all_alleles___Fixture() :
            sample_info(0, &default_genotyping_options){}

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

    EXPECT_CALL(sample_info, get_total_mean_coverage_over_all_alleles_given_a_mininum_threshold(5))
    .Times(1)
    .WillOnce(Return(100));

    //allele 0 - min coverage threshold not satisfied
    EXPECT_CALL(sample_info, get_total_mean_coverage_given_a_mininum_threshold(0, 5))
            .Times(1)
            .WillOnce(Return(0));

    EXPECT_CALL(sample_info, get_gaps(0))
            .Times(1)
            .WillOnce(Return(0.3));

    EXPECT_CALL(sample_info, compute_likelihood(false, default_genotyping_options.get_sample_index_to_exp_depth_covg()[0],
            0, 100, DoubleEq(default_genotyping_options.get_error_rate()), DoubleEq(0.3)))
            .Times(1)
            .WillOnce(Return(-50.5));

    //allele 1 - min coverage threshold not satisfied
    EXPECT_CALL(sample_info, get_total_mean_coverage_given_a_mininum_threshold(1, 5))
            .Times(1)
            .WillOnce(Return(100));

    EXPECT_CALL(sample_info, get_gaps(1))
            .Times(1)
            .WillOnce(Return(0.5));

    EXPECT_CALL(sample_info, compute_likelihood(true, default_genotyping_options.get_sample_index_to_exp_depth_covg()[1],
                                                100, 0, DoubleEq(default_genotyping_options.get_error_rate()), DoubleEq(0.5)))
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
        MOCK_METHOD(std::vector<double>, get_likelihoods_for_all_alleles, (), (const));
        MOCK_METHOD(uint32_t, get_mean_coverage_both_alleles, (uint32_t allele), (const));
    };

    SampleInfoTest___get_confidence___Fixture() :
            sample_info(0, &genotyping_options_get_confidence){}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfoMock sample_info;
};
GenotypingOptions SampleInfoTest___get_confidence___Fixture::genotyping_options_get_confidence
    ({10,10}, 0.01, 0, 0, 0.0, 50, 100, 0, 0);

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

    SampleInfo::IndexAndConfidenceAndMaxLikelihood expected = std::make_tuple(1, 2.0, -3.0);
    EXPECT_EQ(expected, actual);
}


class SampleInfoTest___get_confidence_to_string___Fixture : public ::testing::Test {
public:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(boost::optional<IndexAndConfidenceAndMaxLikelihood>, get_confidence, (), (const));
    };

    SampleInfoTest___get_confidence_to_string___Fixture() :
            sample_info(0, &default_genotyping_options){}

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
            .WillOnce(Return(std::make_tuple((size_t)1, 50.5, -200.2)));

    EXPECT_EQ("50.5", sample_info.get_confidence_to_string());
}


class SampleInfoTest___get_genotype_from_coverage___Fixture : public ::testing::Test {
public:
    static GenotypingOptions genotyping_options_high_confidence_threshold;
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(boost::optional<SampleInfo::IndexAndConfidenceAndMaxLikelihood>, get_confidence, (), (const));
    };

    SampleInfoTest___get_genotype_from_coverage___Fixture() :
            sample_info(0, &genotyping_options_high_confidence_threshold){}

    void SetUp() override {
    }

    void TearDown() override {
    }

    SampleInfoMock sample_info;
};
GenotypingOptions SampleInfoTest___get_genotype_from_coverage___Fixture::genotyping_options_high_confidence_threshold
        ({10,10}, 0.01, 100, 0, 0.0, 0, 0, 0, 0);

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
            .WillOnce(Return(std::make_tuple((size_t)1, 99.0, -50.5)));

    auto actual = sample_info.get_genotype_from_coverage();

    EXPECT_EQ(boost::none, actual);
}

TEST_F(SampleInfoTest___get_genotype_from_coverage___Fixture, valid_confidence_and_above_threshold) {
    EXPECT_CALL(sample_info, get_confidence)
            .Times(1)
            .WillOnce(Return(std::make_tuple((size_t)1, 105.0, -50.5)));

    auto actual = sample_info.get_genotype_from_coverage();

    std::pair<uint32_t, double> expected(1, -50.5);
    EXPECT_EQ(expected, actual);
}