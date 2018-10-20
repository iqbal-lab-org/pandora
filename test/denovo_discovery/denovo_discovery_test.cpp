#include "gtest/gtest.h"
#include "denovo_discovery/denovo_discovery.h"


TEST(ExpectedKmerCoverage, IncreasingErrorRates_DecreasingExpectedCovg) {
    const uint32_t read_covg{50};
    const uint32_t ref_length{500};
    const uint32_t k{11};
    const std::vector<double> error_rates = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
                                             0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95};

    const std::vector<double> expected = {49.0, 27.871204521546524, 15.376719208410002, 8.199818940791095,
                                          4.209067950080002, 2.06952166557312, 0.9688901040699994,
                                          0.4287883755264308, 0.17777055743999992, 0.06826304619110846,
                                          0.02392578125, 0.0075081636759814375, 0.002055208960000001,
                                          0.0004730908711279294, 8.680203000000014e-05, 1.1682510375976562e-05,
                                          1.0035199999999976e-06, 4.238380371093757e-08, 4.899999999999988e-10,
                                          2.3925781250000235e-13};

    assert(error_rates.size() == expected.size());

    for (uint32_t i = 0; i < error_rates.size(); i++) {
        const auto &e{error_rates[i]};
        const auto result{denovo_discovery::calculate_kmer_coverage(read_covg, ref_length, k, e)};
        EXPECT_DOUBLE_EQ(result, expected[i]);
    }
}


TEST(ExpectedKmerCoverage, ZeroReadCovg_ReturnZero) {
    const uint32_t read_covg{0};
    const uint32_t ref_length{500};
    const uint32_t k{11};
    const double error_rate{0.1};

    const auto result{denovo_discovery::calculate_kmer_coverage(read_covg, ref_length, k, error_rate)};
    const double expected{0};

    EXPECT_DOUBLE_EQ(result, expected);
}


TEST(ExpectedKmerCoverage, ZeroRefLength_ThrowException) {
    const uint32_t read_covg{50};
    const uint32_t ref_length{0};
    const uint32_t k{11};
    const double error_rate{0.1};

    const std::string expected_msg{"ref_length should be greater than 0."};

    try {
        const auto result{denovo_discovery::calculate_kmer_coverage(read_covg, ref_length, k, error_rate)};
        FAIL() << "Expected std::invalid_argument. Instead, got result: " << std::to_string(result);
    }
    catch (std::invalid_argument const &err) {
        EXPECT_EQ(err.what(), expected_msg);
    }
    catch (...) {
        FAIL() << "Expected std::invalid_argument";
    }
}


TEST(ExpectedKmerCoverage, ZeroK_ThrowException) {
    const uint32_t read_covg{50};
    const uint32_t ref_length{500};
    const uint32_t k{0};
    const double error_rate{0.1};

    const std::string expected_msg{"K should be greater than 0."};

    try {
        const auto result{denovo_discovery::calculate_kmer_coverage(read_covg, ref_length, k, error_rate)};
        FAIL() << "Expected std::invalid_argument. Instead, got result: " << std::to_string(result);
    }
    catch (std::invalid_argument const &err) {
        EXPECT_EQ(err.what(), expected_msg);
    }
    catch (...) {
        FAIL() << "Expected std::invalid_argument";
    }
}