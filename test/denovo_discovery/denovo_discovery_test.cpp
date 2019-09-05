#include "gtest/gtest.h"
#include "denovo_discovery/denovo_discovery.h"


TEST(ExpectedKmerCoverage, IncreasingErrorRatesDecreasingExpectedCovg) {
    const uint32_t read_covg { 50 };
    const uint32_t ref_length { 500 };
    const uint32_t k { 11 };
    const std::vector<double> error_rates = { 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
                                              0.7, 0.75, 0.8, 0.85, 0.9, 0.95 };

    const std::vector<double> expected = { 49.0, 27.871204521546524, 15.376719208410002, 8.199818940791095,
                                           4.209067950080002, 2.06952166557312, 0.9688901040699994, 0.4287883755264308,
                                           0.17777055743999992, 0.06826304619110846, 0.02392578125,
                                           0.0075081636759814375, 0.002055208960000001, 0.0004730908711279294,
                                           8.680203000000014e-05, 1.1682510375976562e-05, 1.0035199999999976e-06,
                                           4.238380371093757e-08, 4.899999999999988e-10, 2.3925781250000235e-13 };

    for (uint32_t i = 0; i < error_rates.size(); i++) {
        const auto &e { error_rates[i] };
        const DenovoDiscovery denovo { k, e };
        const auto result { denovo.calculate_kmer_coverage(read_covg, ref_length) };
        EXPECT_DOUBLE_EQ(result, expected[i]);
    }
}


TEST(ExpectedKmerCoverage, ZeroReadCovgReturnsZero) {
    const uint32_t read_covg { 0 };
    const uint32_t ref_length { 500 };
    const uint32_t k { 11 };
    const double error_rate { 0.1 };
    const DenovoDiscovery denovo { k, error_rate };

    const auto result { denovo.calculate_kmer_coverage(read_covg, ref_length) };
    const double expected { 0 };

    EXPECT_DOUBLE_EQ(result, expected);
}


TEST(ExpectedKmerCoverage, ZeroRefLengthThrowsException) {
    const uint32_t read_covg { 50 };
    const uint32_t ref_length { 0 };
    const uint32_t k { 11 };
    const double error_rate { 0.1 };
    const DenovoDiscovery denovo { k, error_rate };

    const std::string expected_msg { "ref_length should be greater than 0." };

    try {
        const auto result { denovo.calculate_kmer_coverage(read_covg, ref_length) };
        FAIL() << "Expected std::invalid_argument. Instead, got result: " << std::to_string(result);
    } catch (std::invalid_argument const &err) {
        EXPECT_EQ(err.what(), expected_msg);
    } catch (...) {
        FAIL() << "Expected std::invalid_argument";
    }
}


TEST(ExpectedKmerCoverage, ZeroKThrowsException) {
    const uint32_t read_covg { 50 };
    const uint32_t ref_length { 500 };
    const uint32_t k { 0 };
    const double error_rate { 0.1 };
    const DenovoDiscovery denovo { k, error_rate };

    const std::string expected_msg { "K should be greater than 0." };

    try {
        const auto result { denovo.calculate_kmer_coverage(read_covg, ref_length) };
        FAIL() << "Expected std::invalid_argument. Instead, got result: " << std::to_string(result);
    } catch (std::invalid_argument const &err) {
        EXPECT_EQ(err.what(), expected_msg);
    } catch (...) {
        FAIL() << "Expected std::invalid_argument";
    }
}


TEST(ExpectedKmerCoverage, negativeErrorRateThrowsException) {
    const uint32_t read_covg { 50 };
    const uint32_t ref_length { 500 };
    const uint32_t k { 5 };
    const double error_rate { -0.1 };
    const DenovoDiscovery denovo { k, error_rate };

    const std::string expected_msg { "error_rate should not be a negative value." };

    try {
        const auto result { denovo.calculate_kmer_coverage(read_covg, ref_length) };
        FAIL() << "Expected std::invalid_argument. Instead, got result: " << std::to_string(result);
    } catch (std::invalid_argument const &err) {
        EXPECT_EQ(err.what(), expected_msg);
    } catch (...) {
        FAIL() << "Expected std::invalid_argument";
    }
}


TEST(FindPathsThroughCandidateRegionTest, emptyPileupReturnsEmpty) {
    const int k { 9 };
    const double error_rate { 0.11 };
    DenovoDiscovery denovo { k, error_rate };
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "ATGCGCTGAGAGTCGGACT";

    denovo.find_paths_through_candidate_region(candidate_region);

    const DenovoPaths expected;
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);

}


TEST(FindPathsThroughCandidateRegionTest, kmerSizeBiggerThanCandidateReturnsEmpty) {
    const int k { 99 };
    const double error_rate { 0.11 };
    DenovoDiscovery denovo { k, error_rate };
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "ATGCGCTGAGAGTCGGACT";
    candidate_region.pileup = { "FOO", "BAR" };

    denovo.find_paths_through_candidate_region(candidate_region);

    const DenovoPaths expected;
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);

}


TEST(FindPathsThroughCandidateRegionTest, passInDataThatCausesGatbErrorNoFailAndReturnsEmpty) {
    const int k { 9 };
    const double error_rate { 0.11 };
    DenovoDiscovery denovo { k, error_rate };
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "ATGCGCTGAGAGTCGGACT";
    candidate_region.pileup = { "FOO", "BAR" };

    denovo.find_paths_through_candidate_region(candidate_region);

    const DenovoPaths expected;
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);

}


TEST(FindPathsThroughCandidateRegionTest, startKmersDontExistInGraphReturnEmpty) {
    const int k { 9 };
    const double error_rate { 0.11 };
    DenovoDiscovery denovo { k, error_rate };
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "GGGGGGGGGGAGTCGGACT";
    candidate_region.pileup = { "ATGCGCTGAGAGTCGGACT", "ATGCGCTGAGAGTCGGACT" };

    denovo.find_paths_through_candidate_region(candidate_region);

    const DenovoPaths expected;
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);
}


TEST(FindPathsThroughCandidateRegionTest, endKmersDontExistInGraphReturnEmpty) {
    const int k { 9 };
    const double error_rate { 0.11 };
    DenovoDiscovery denovo { k, error_rate };
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "ATGCGCTGAGCCCCCCCCC";
    candidate_region.pileup = { "ATGCGCTGAGAGTCGGACT", "ATGCGCTGAGAGTCGGACT" };

    denovo.find_paths_through_candidate_region(candidate_region);

    const DenovoPaths expected;
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);
}


TEST(FindPathsThroughCandidateRegionTest, endKmerExistsInStartKmersFindPathAndCycles) {
    const int k { 9 };
    const double error_rate { 0.11 };
    const uint8_t max_insertion_size = 50;
    DenovoDiscovery denovo { k, error_rate, max_insertion_size };
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "ATGCGCTGAGATGCGCTGA";
    candidate_region.pileup = { "ATGCGCTGACATGCGCTGA", "ATGCGCTGACATGCGCTGA" };

    denovo.find_paths_through_candidate_region(candidate_region);

    const DenovoPaths expected { "ATGCGCTGACATGCGCTGA", "ATGCGCTGACATGCGCTGACATGCGCTGA",
                                 "ATGCGCTGACATGCGCTGACATGCGCTGACATGCGCTGA",
                                 "ATGCGCTGACATGCGCTGACATGCGCTGACATGCGCTGACATGCGCTGA",
                                 "ATGCGCTGACATGCGCTGACATGCGCTGACATGCGCTGACATGCGCTGACATGCGCTGA",
                                 "ATGCGCTGACATGCGCTGACATGCGCTGACATGCGCTGACATGCGCTGACATGCGCTGACATGCGCTGA" };
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);
}


TEST(FindPathsThroughCandidateRegionTest, doGraphCleaningtwoIdenticalReadsPlusNoiseReturnOnePath) {
    const int k { 9 };
    const double error_rate { 0.11 };
    DenovoDiscovery denovo { k, error_rate };
    denovo.clean_assembly_graph = true;
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "ATGCGCTGAGAGTCGGACT";
    candidate_region.pileup = { "ATGCGCTGAGAGTCGGACT", "ATGCGCTGAGAGTCGGACT", "AAATAAA", "GCGGCGCGGCC" };

    denovo.find_paths_through_candidate_region(candidate_region);

    const DenovoPaths expected { "ATGCGCTGAGAGTCGGACT" };
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);
}


TEST(FindPathsThroughCandidateRegionTest, twoIdenticalReadsReturnOnePath) {
    const int k { 9 };
    const double error_rate { 0.11 };
    DenovoDiscovery denovo { k, error_rate };
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "ATGCGCTGAGAGTCGGACT";
    candidate_region.pileup = { "ATGCGCTGAGAGTCGGACT", "ATGCGCTGAGAGTCGGACT" };

    denovo.find_paths_through_candidate_region(candidate_region);

    const DenovoPaths expected { "ATGCGCTGAGAGTCGGACT" };
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);
}


TEST(FindPathsThroughCandidateRegionTest, twoPossiblePathsWithLowCovgOnBothReturnsNoPath) {
    const int k { 9 };
    const double error_rate { 0.11 };
    DenovoDiscovery denovo { k, error_rate };
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "ATGCGCTGAGAGTCGGACT";
    candidate_region.pileup = { "ATGCGCTGAGAGTCGGACT", "ATGCGCTGATAGTCGGACT" };

    denovo.find_paths_through_candidate_region(candidate_region);

    const DenovoPaths expected;
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);

}


TEST(FindPathsThroughCandidateRegionTest, twoPossiblePathsOneWithLowCovgOnOneReturnsOnePath) {
    const int k { 9 };
    const double error_rate { 0.11 };
    DenovoDiscovery denovo { k, error_rate };
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "ATGCGCTGAGAGTCGGACT";
    candidate_region.pileup = { "ATGCGCTGAGAGTCGGACT", "ATGCGCTGAGAGTCGGACT", "ATGCGCTGATAGTCGGACT" };

    denovo.find_paths_through_candidate_region(candidate_region);

    const DenovoPaths expected { "ATGCGCTGAGAGTCGGACT" };
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);

}


TEST(FindPathsThroughCandidateRegionTest, twoPossiblePathsWithGoodCovgReturnsTwoPaths) {
    const int k { 9 };
    const double error_rate { 0.11 };
    DenovoDiscovery denovo { k, error_rate };
    CandidateRegion candidate_region { Interval(0, 1), "test" };
    candidate_region.max_likelihood_sequence = "ATGCGCTGAGAGTCGGACT";
    candidate_region.pileup = { "ATGCGCTGAGAGTCGGACT", "ATGCGCTGAGAGTCGGACT", "ATGCGCTGATAGTCGGACT",
                                "ATGCGCTGATAGTCGGACT" };

    denovo.find_paths_through_candidate_region(candidate_region);
    std::sort(candidate_region.denovo_paths.begin(), candidate_region.denovo_paths.end());

    const DenovoPaths expected { "ATGCGCTGAGAGTCGGACT", "ATGCGCTGATAGTCGGACT" };
    const auto &actual { candidate_region.denovo_paths };

    EXPECT_EQ(actual, expected);

}