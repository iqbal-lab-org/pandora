#include <stdint.h>
#include <iostream>
#include "gtest/gtest.h"
#include "pangenome/pangraph.h"
#include "estimate_parameters.h"
#include "test_helpers_containers.h"
#include "test_helpers.h"

using namespace std;

const std::string TEST_CASE_DIR = "../../test/test_cases/";

TEST(EstimateParameters_FitMeanCovg, Empty)
{
    std::vector<uint> v = {};
    EXPECT_NEAR(0, fit_mean_covg(v, 0), 0.000001);
}

TEST(EstimateParameters_FitMeanCovg, SimpleMean)
{
    std::vector<uint> v = { 30, 24, 12, 3, 6, 2, 14, 15, 16, 18 };
    EXPECT_NEAR(4.071428571, fit_mean_covg(v, 0), 0.000001);
}

TEST(EstimateParameters_FitMeanCovg, SimpleMean_ZeroThreshOne)
{
    std::vector<uint> v = { 30, 24, 12, 3, 6, 2, 14, 15, 16, 18 };
    EXPECT_NEAR(5.181818182, fit_mean_covg(v, 1), 0.000001);
}

TEST(EstimateParameters_FitMeanCovg, SimpleMean_ZeroThreshTwo)
{
    std::vector<uint> v = { 30, 24, 12, 3, 6, 2, 14, 15, 16, 18 };
    EXPECT_NEAR(6.348837209, fit_mean_covg(v, 2), 0.000001);
}

TEST(EstimateParameters_FitVarianceCovg, Empty)
{
    std::vector<uint> v = {};
    double mean = 0;
    EXPECT_NEAR(0, fit_variance_covg(v, mean, 0), 0.000001);
}

TEST(EstimateParameters_FitVarianceCovg, SimpleVar)
{
    std::vector<uint> v = { 30, 24, 12, 3, 6, 2, 14, 15, 16, 18 };
    double mean = fit_mean_covg(v, 0);
    EXPECT_NEAR(11.75204082, fit_variance_covg(v, mean, 0), 0.000001);
}

TEST(EstimateParameters_FitVarianceCovg, SimpleVar_ZeroThreshOne)
{
    std::vector<uint> v = { 30, 24, 12, 3, 6, 2, 14, 15, 16, 18 };
    double mean = fit_mean_covg(v, 1);
    EXPECT_NEAR(9.203305785, fit_variance_covg(v, mean, 1), 0.000001);
}

TEST(EstimateParameters_FitVarianceCovg, SimpleVar_ZeroThreshTwo)
{
    std::vector<uint> v = { 30, 24, 12, 3, 6, 2, 14, 15, 16, 18 };
    double mean = fit_mean_covg(v, 2);
    EXPECT_NEAR(5.529475392, fit_variance_covg(v, mean, 2), 0.000001);
}

TEST(EstimateParameters_FitNegativeBinomial, MeanZero_FatalRuntimeError)
{
    double mean = 0, variance = 1;
    float p, r;
    ASSERT_EXCEPTION(fit_negative_binomial(mean, variance, p, r), FatalRuntimeError,
        "Negative binomial parameters are invalid");
}

TEST(EstimateParameters_FitNegativeBinomial, VarianceZero_FatalRuntimeError)
{
    double mean = 1, variance = 0;
    float p, r;
    ASSERT_EXCEPTION(fit_negative_binomial(mean, variance, p, r), FatalRuntimeError,
        "Negative binomial parameters are invalid");
}

TEST(EstimateParameters_FitNegativeBinomial, MeanVarianceEqual_FatalRuntimeError)
{
    double mean = 1, variance = 1;
    float p, r;
    ASSERT_EXCEPTION(fit_negative_binomial(mean, variance, p, r), FatalRuntimeError,
        "Negative binomial parameters are invalid");
}

TEST(EstimateParameters_FitNegativeBinomial, MeanGreaterThanVariance_FatalRuntimeError)
{
    double mean = 2, variance = 1;
    float p, r;
    ASSERT_EXCEPTION(fit_negative_binomial(mean, variance, p, r), FatalRuntimeError,
        "Negative binomial parameters are invalid");
}

TEST(EstimateParameters_FitNegativeBinomial, SimpleFit)
{
    double mean = 1, variance = 2;
    float p, r;
    fit_negative_binomial(mean, variance, p, r);
    EXPECT_FLOAT_EQ(p, 0.5);
    EXPECT_FLOAT_EQ(r, 1);
}

TEST(EstimateParameters_FitNegativeBinomial, Fit)
{
    double mean = 1, variance = 4;
    float p, r;
    fit_negative_binomial(mean, variance, p, r);
    EXPECT_FLOAT_EQ(p, 0.25);
    EXPECT_FLOAT_EQ(r, (float)1 / 3);
}

TEST(EstimateParameters_FindMeanCovg, Examples)
{
    // NB this finds the position in vector at which max of the second peak occurs
    std::vector<uint> v1 = { 30, 24, 12, 3, 6, 2, 14, 15, 16, 18, 40, 26, 35, 14 };
    EXPECT_EQ(uint(10), find_mean_covg(v1));
    // not thrown by a single weird one not near peak
    std::vector<uint> v2 = { 30, 24, 12, 3, 70, 2, 14, 15, 16, 18, 40, 26, 35, 14 };
    EXPECT_EQ(uint(10), find_mean_covg(v2));
    // doesn't matter if second peak much lower
    std::vector<uint> v3 = { 30, 24, 12, 3, 6, 2, 14, 15, 16, 18, 14, 8, 9, 1 };
    EXPECT_EQ(uint(9), find_mean_covg(v3));
    // do need an increase three times
    std::vector<uint> v4 = { 30, 24, 12, 3, 6, 2, 11, 10, 9, 8, 4, 3, 2, 1 };
    EXPECT_EQ(uint(0), find_mean_covg(v4));
    // a too clean one with run of zeroes doing it
    std::vector<uint> v5 = { 30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 5 };
    EXPECT_EQ(uint(12), find_mean_covg(v5));
    // insufficient zeroes with a peak
    std::vector<uint> v6
        = { 30, 0, 0, 0, 0, 0, 0, 0, 0, 2, 14, 15, 16, 18, 14, 8, 9, 1 };
    EXPECT_EQ(uint(13), find_mean_covg(v6));
}

TEST(EstimateParameters_FindProbThresh, Empty_Zero)
{
    // NB this finds the position in vector at which min occurs between 2 peaks
    std::vector<uint> v = {};
    EXPECT_EQ(0, find_prob_thresh(v));
}

TEST(EstimateParameters_FindProbThresh, SimpleCases)
{
    // NB this finds the position in vector at which min occurs between 2 peaks
    std::vector<uint> v1
        = { 30, 24, 18, 16, 12, 3, 6, 2, 1, 15, 16, 18, 12, 26, 35, 40 };
    EXPECT_EQ(8 - 200, find_prob_thresh(v1));
    // not thrown by low values outside of valley
    std::vector<uint> v2 = { 1, 30, 24, 12, 3, 6, 2, 0, 15, 16, 18, 12, 26, 35, 40, 0 };
    EXPECT_EQ(7 - 200, find_prob_thresh(v2));
}

TEST(EstimateParameters_FindProbThresh, SecondPeakClose)
{
    std::vector<uint> v = { 30, 24, 12, 3, 6, 2, 1, 15, 16, 18, 12, 26, 35, 40 };
    EXPECT_EQ(6 - 200, find_prob_thresh(v));
}

TEST(EstimateParameters_FindProbThresh, OnePeak)
{
    std::vector<uint> v = { 30, 24, 12, 3, 3, 2, 0, 0 };
    EXPECT_EQ(5 - 200, find_prob_thresh(v));
}

void setup_index_estimate_parameters(
    std::vector<std::shared_ptr<LocalPRG>>& prgs, std::shared_ptr<Index>& index)
{
    auto s = std::make_shared<LocalPRG>(LocalPRG(
        0, "nested varsite", "AAA 5 GGGG 7 CCC 8 TTT 7  6 GG 5 TTTTTACAGACGT"));
    prgs.emplace_back(s);
    s = std::make_shared<LocalPRG>(LocalPRG(
        1, "much more complex", "TCATTC 5 ACTC 7 TAGTCA 8 TTGTGA 7  6 AACTAG 5 AGCTG"));
    prgs.emplace_back(s);
    s = std::make_shared<LocalPRG>(LocalPRG(2,
        "one with lots of null at start and end, and a long stretch in between",
        " 5  7  9  11 AGTTCTGAAACATTGCGCGTGAGATCTCTG 12 T 11  10 A 9  8 C 7  6 G 5 "));
    prgs.emplace_back(s);
    prgs[0]->minimizer_sketch(index.get(), 1, 6);
    prgs[1]->minimizer_sketch(index.get(), 1, 6);
    prgs[2]->minimizer_sketch(index.get(), 1, 6);
}

TEST(EstimateParameters_EstimateParameters, NoPangraphNodes)
{
    auto outdir = "test";
    uint32_t k = 6, covg = 10, sample_id = 0;
    float e_rate = 0.01;
    bool bin = false;

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto expected_depth_covg
        = estimate_parameters(pangraph, outdir, k, e_rate, covg, bin, sample_id);
    EXPECT_EQ(expected_depth_covg, (uint)10);
}

TEST(EstimateParameters_EstimateParameters, PangraphWithNodes_SimpleBinomial)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    auto index = std::make_shared<Index>();
    setup_index_estimate_parameters(prgs, index);

    auto outdir = "estimate_parameters_test";
    uint32_t w = 1, k = 6, covg = 10, sample_id = 0, min_cluster_size = 1,
             genome_size = 22;
    float e_rate = 0.01;
    bool bin = true, illumina = true;

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    const auto filepath = TEST_CASE_DIR + "estimate_parameters_reads.fa";
    const std::string sid{"sample"};
    const SampleData sample = std::make_pair(sid, filepath);

    fs::create_directories(outdir);
    pangraph_from_read_file(sample, pangraph, index, prgs, w, k, 1, e_rate,
        outdir, min_cluster_size, genome_size, illumina);
    pangraph->add_hits_to_kmergraphs();
    fs::remove_all(outdir);

    auto expected_depth_covg
        = estimate_parameters(pangraph, outdir, k, e_rate, covg, bin, sample_id);
    EXPECT_NEAR(expected_depth_covg, (uint)4, 1);
}

TEST(
    EstimateParameters_EstimateParameters, PangraphWithNodes_SimpleBinomial_LowCoverage)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    auto index = std::make_shared<Index>();
    setup_index_estimate_parameters(prgs, index);

    auto outdir = "estimate_parameters_test";
    uint32_t w = 1, k = 6, covg = 10, sample_id = 0, min_cluster_size = 1,
             genome_size = 22;
    float e_rate = 0.01;
    bool bin = true, illumina = true;

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    const auto filepath = TEST_CASE_DIR + "estimate_parameters_reads3.fa";
    SampleData sample_data{"estimate_parameters_reads3", filepath};

    fs::create_directories(outdir);
    pangraph_from_read_file(sample_data, pangraph, index, prgs, w, k, 1, e_rate,
        outdir, min_cluster_size, genome_size, illumina);
    pangraph->add_hits_to_kmergraphs();
    fs::remove_all(outdir);

    auto expected_depth_covg
        = estimate_parameters(pangraph, outdir, k, e_rate, covg, bin, sample_id);
    EXPECT_NEAR(expected_depth_covg, (uint)2, 1);
}

TEST(EstimateParameters_EstimateParameters, PangraphWithNodes_SimpleNegativeBinomial)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    auto index = std::make_shared<Index>();
    setup_index_estimate_parameters(prgs, index);

    auto outdir = "estimate_parameters_test";
    uint32_t w = 1, k = 6, covg = 10, sample_id = 0, min_cluster_size = 1,
             genome_size = 22;
    float e_rate = 0.01;
    bool bin = false, illumina = true;

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    const auto filepath = TEST_CASE_DIR + "estimate_parameters_reads.fa";
    SampleData sample_data{"estimate_parameters_reads", filepath};

    fs::create_directories(outdir);
    pangraph_from_read_file(sample_data, pangraph, index, prgs, w, k, 1, e_rate,
        outdir, min_cluster_size, genome_size, illumina);
    pangraph->add_hits_to_kmergraphs();
    fs::remove_all(outdir);

    auto expected_depth_covg
        = estimate_parameters(pangraph, outdir, k, e_rate, covg, bin, sample_id);
    EXPECT_NEAR(expected_depth_covg, (uint)4, 1);
}

TEST(EstimateParameters_EstimateParameters,
    PangraphWithNodes_SimpleNegativeBinomial_HighCoverage)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    auto index = std::make_shared<Index>();
    setup_index_estimate_parameters(prgs, index);

    auto outdir = "estimate_parameters_test";
    uint32_t w = 1, k = 6, covg = 32, sample_id = 0, min_cluster_size = 1,
             genome_size = 22;
    float e_rate = 0.01;
    bool bin = false, illumina = true;

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    const auto filepath = TEST_CASE_DIR + "estimate_parameters_reads4.fa";
    SampleData sample_data{"estimate_parameters_reads4", filepath};

    fs::create_directories(outdir);
    pangraph_from_read_file(sample_data, pangraph, index, prgs, w, k, 1, e_rate,
        outdir, min_cluster_size, genome_size, illumina);
    pangraph->add_hits_to_kmergraphs();
    fs::remove_all(outdir);

    auto expected_depth_covg
        = estimate_parameters(pangraph, outdir, k, e_rate, covg, bin, sample_id);
    EXPECT_NEAR(
        expected_depth_covg, (uint)40, 1); // NB method overestimates covg, true is 32
}

TEST(EstimateParameters_EstimateParameters,
    PangraphWithNodes_SimpleBinomial_HighCoverage)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    auto index = std::make_shared<Index>();
    setup_index_estimate_parameters(prgs, index);

    auto outdir = "estimate_parameters_test";
    uint32_t w = 1, k = 6, covg = 32, sample_id = 0, min_cluster_size = 1,
             genome_size = 22;
    float e_rate = 0.01;
    bool bin = true, illumina = true;

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    const auto filepath = TEST_CASE_DIR + "estimate_parameters_reads4.fa";
    SampleData sample_data{"estimate_parameters_reads4_2", filepath};

    fs::create_directories(outdir);
    pangraph_from_read_file(sample_data, pangraph, index, prgs, w, k, 1, e_rate,
        outdir, min_cluster_size, genome_size, illumina);
    pangraph->add_hits_to_kmergraphs();
    fs::remove_all(outdir);

    auto expected_depth_covg
        = estimate_parameters(pangraph, outdir, k, e_rate, covg, bin, sample_id);
    EXPECT_NEAR(
        expected_depth_covg, (uint)32, 1); // NB method overestimates covg, true is 32
}

TEST(EstimateParameters_EstimateParameters, PangraphWithNodes_NoiseReads)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    auto index = std::make_shared<Index>();
    setup_index_estimate_parameters(prgs, index);

    auto outdir = "estimate_parameters_test";
    uint32_t w = 1, k = 6, covg = 10, sample_id = 0, min_cluster_size = 1,
             genome_size = 22;
    float e_rate = 0.01;
    bool bin = false, illumina = true;

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    const auto filepath = TEST_CASE_DIR + "estimate_parameters_reads2.fa";
    SampleData sample_data{"estimate_parameters_reads_2", filepath};

    fs::create_directories(outdir);
    pangraph_from_read_file(sample_data, pangraph, index, prgs, w, k, 1, e_rate,
        outdir, min_cluster_size, genome_size, illumina);
    pangraph->add_hits_to_kmergraphs();
    fs::remove_all(outdir);

    auto expected_depth_covg
        = estimate_parameters(pangraph, outdir, k, e_rate, covg, bin, sample_id);
    EXPECT_NEAR(expected_depth_covg, (uint)5, 1);
}
