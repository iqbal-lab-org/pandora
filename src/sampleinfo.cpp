#include "sampleinfo.h"

double SampleInfo::get_gaps (uint32_t allele, uint32_t threshold_to_consider_a_gap) const {
    double gaps = 0.0;
    size_t size = allele_to_forward_coverages[allele].size();

    for (size_t i = 0; i < size; ++i) {
        if (allele_to_forward_coverages[allele][i] + allele_to_reverse_coverages[allele][i] < threshold_to_consider_a_gap)
            gaps++;
    }

    return gaps / size;
}


std::vector<double> SampleInfo::get_likelihoods_for_all_alleles (uint32_t expected_depth_covg, double error_rate, uint32_t min_allele_covg,
                                     double min_fraction_allele_covg, uint32_t min_kmer_covg) const {
    std::vector<double> likelihoods;

    uint32_t min_coverage_threshold = std::max(min_allele_covg, uint(min_fraction_allele_covg*expected_depth_covg));
    double total_mean_coverage_over_all_alleles_above_threshold = get_total_mean_coverage_over_all_alleles_given_a_mininum_threshold(min_coverage_threshold);

    for (size_t allele = 0; allele < get_number_of_alleles(); ++allele) {
        uint32_t total_mean_coverage_of_allele_above_threshold = get_total_mean_coverage_given_a_mininum_threshold(allele, min_coverage_threshold);
        uint32_t total_mean_coverage_of_all_other_alleles_above_threshold = total_mean_coverage_over_all_alleles_above_threshold - total_mean_coverage_of_allele_above_threshold;
        bool min_coverage_threshold_is_satisfied = total_mean_coverage_of_allele_above_threshold > 0;
        double gaps = get_gaps(allele, min_kmer_covg);

        double likelihood = compute_likelihood(min_coverage_threshold_is_satisfied, expected_depth_covg, total_mean_coverage_of_allele_above_threshold,
                                               total_mean_coverage_of_all_other_alleles_above_threshold, error_rate, gaps);

        likelihoods.push_back(likelihood);
    }

    return likelihoods;
}


uint32_t SampleInfo::get_total_mean_coverage_given_a_mininum_threshold(uint32_t allele, uint32_t minimum_threshold) const {
    uint32_t total_mean_coverage = this->get_total_mean_coverage(allele);
    if (total_mean_coverage >= minimum_threshold)
        return total_mean_coverage;
    else
        return 0;
}


uint32_t SampleInfo::get_total_mean_coverage_over_all_alleles_given_a_mininum_threshold(uint32_t minimum_threshold) const {
    uint32_t total_mean_coverage_over_all_alleles = 0;
    for (uint32_t allele = 0; allele < get_number_of_alleles(); ++allele)
        total_mean_coverage_over_all_alleles += get_total_mean_coverage_given_a_mininum_threshold(allele, minimum_threshold);
    return total_mean_coverage_over_all_alleles;
}


double SampleInfo::compute_likelihood(bool min_coverage_threshold_is_satisfied, uint32_t expected_depth_covg, uint32_t total_mean_coverage_of_allele_above_threshold,
                                 uint32_t total_mean_coverage_of_all_other_alleles_above_threshold, double error_rate,
                                 double gaps) {
    // For more details about this, look at page 66 of the thesis
    if (min_coverage_threshold_is_satisfied) {
        return       - expected_depth_covg
                     + total_mean_coverage_of_allele_above_threshold * log(expected_depth_covg)
                     - Maths::logfactorial(total_mean_coverage_of_allele_above_threshold)
                     + total_mean_coverage_of_all_other_alleles_above_threshold * log(error_rate)
                     - expected_depth_covg * gaps
                     + log(1 - exp(-(double)expected_depth_covg)) * (1 - gaps);
    } else {
        return       - expected_depth_covg
                     + total_mean_coverage_of_all_other_alleles_above_threshold * log(error_rate)
                     - expected_depth_covg * gaps
                     + log(1 - exp(-(double)expected_depth_covg)) * (1 - gaps);
    }
}


boost::optional<SampleInfo::IndexAndConfidence> SampleInfo::get_confidence (uint32_t expected_depth_covg, double error_rate, uint32_t min_allele_covg,
                       double min_fraction_allele_covg, uint32_t min_kmer_covg,
                       const uint32_t &min_total_covg, const uint32_t &min_diff_covg) const {
    std::vector<double> likelihoods_for_all_alleles = get_likelihoods_for_all_alleles(expected_depth_covg, error_rate, min_allele_covg,
                                                                                      min_fraction_allele_covg, min_kmer_covg);

    size_t index_of_max_likelihood = Maths::arg_max(likelihoods_for_all_alleles.begin(), likelihoods_for_all_alleles.end());
    size_t index_of_second_max_likelihood; //guaranteed to exist, since we have at least two alleles
    {
        std::vector<double> likelihoods_for_all_alleles_without_max_element(likelihoods_for_all_alleles.begin(), likelihoods_for_all_alleles.end());
        likelihoods_for_all_alleles_without_max_element.erase(likelihoods_for_all_alleles.begin() + index_of_max_likelihood);
        index_of_second_max_likelihood = Maths::arg_max(likelihoods_for_all_alleles_without_max_element.begin(), likelihoods_for_all_alleles_without_max_element.end());
    }

    double max_likelihood = likelihoods_for_all_alleles[index_of_max_likelihood];
    double second_max_likelihood = likelihoods_for_all_alleles[index_of_second_max_likelihood];
    double confidence = std::abs(max_likelihood - second_max_likelihood);

    uint32_t total_mean_coverage_of_max_likelihood_allele = get_total_mean_coverage(index_of_max_likelihood);
    uint32_t total_mean_coverage_of_second_max_likelihood_allele = get_total_mean_coverage(index_of_second_max_likelihood);
    bool enough_total_covg = (total_mean_coverage_of_max_likelihood_allele + total_mean_coverage_of_second_max_likelihood_allele >= min_total_covg);
    bool enough_difference_in_covg = (std::abs((int64_t)(total_mean_coverage_of_max_likelihood_allele) -
                                               (int64_t)(total_mean_coverage_of_second_max_likelihood_allele))
                                      >= min_diff_covg);

    if (enough_total_covg and enough_difference_in_covg)
        return std::make_pair(index_of_max_likelihood, confidence);
    else
        return boost::none;
}


void SamplesInfos::merge_other_samples_infos_into_this(const SamplesInfos &other) {
    bool samples_infos_to_be_merged_have_the_same_size = this->size() == other.size();
    assert(samples_infos_to_be_merged_have_the_same_size);

    for (size_t sample_index = 0; sample_index < this->size(); ++sample_index) {
        (*this)[sample_index].merge_other_sample_info_into_this(other[sample_index]);
    }
}