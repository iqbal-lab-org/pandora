#include "sampleinfo.h"

void SampleInfo::add_coverage_information (const std::vector< std::vector<uint32_t> > &allele_to_forward_coverages,
                                           const std::vector< std::vector<uint32_t> > &allele_to_reverse_coverages) {
    this->allele_to_forward_coverages = allele_to_forward_coverages;
    this->allele_to_reverse_coverages = allele_to_reverse_coverages;
    check_if_coverage_information_is_correct();
}

void SampleInfo::genotype_from_coverage () {
    check_if_coverage_information_is_correct();

    auto genotype_and_max_likelihood_optional = get_genotype_from_coverage();
    if (genotype_and_max_likelihood_optional) {
        std::tie(GT_from_coverages, likelihood_of_GT_from_coverages) = *genotype_and_max_likelihood_optional;
    }else {
        GT_from_coverages = boost::none;
        likelihood_of_GT_from_coverages = boost::none;
    }
}


bool SampleInfo::check_if_coverage_information_is_correct() const {
    // at least two alleles because a VCF record has a ref and at least one alt
    bool there_are_at_least_two_alleles = allele_to_forward_coverages.size() >= 2 and allele_to_reverse_coverages.size() >= 2;
    bool same_length_vectors = allele_to_forward_coverages.size() == allele_to_reverse_coverages.empty();
    assert(there_are_at_least_two_alleles and same_length_vectors);
}

double SampleInfo::get_gaps (uint32_t allele) const {
    double gaps = 0.0;
    size_t size = allele_to_forward_coverages[allele].size();

    for (size_t i = 0; i < size; ++i) {
        if (allele_to_forward_coverages[allele][i] + allele_to_reverse_coverages[allele][i] < genotyping_options->get_min_kmer_covg())
            gaps++;
    }

    return gaps / size;
}


std::vector<double> SampleInfo::get_likelihoods_for_all_alleles () const {
    std::vector<double> likelihoods;

    uint32_t min_coverage_threshold = std::max(genotyping_options->get_min_allele_covg(), uint(genotyping_options->get_min_fraction_allele_covg() * exp_depth_covg_for_this_sample));
    double total_mean_coverage_over_all_alleles_above_threshold = get_total_mean_coverage_over_all_alleles_given_a_mininum_threshold(min_coverage_threshold);

    for (size_t allele = 0; allele < get_number_of_alleles(); ++allele) {
        uint32_t total_mean_coverage_of_allele_above_threshold = get_total_mean_coverage_given_a_mininum_threshold(allele, min_coverage_threshold);
        uint32_t total_mean_coverage_of_all_other_alleles_above_threshold = total_mean_coverage_over_all_alleles_above_threshold - total_mean_coverage_of_allele_above_threshold;
        bool min_coverage_threshold_is_satisfied = total_mean_coverage_of_allele_above_threshold > 0;
        double gaps = get_gaps(allele);

        double likelihood = compute_likelihood(min_coverage_threshold_is_satisfied, exp_depth_covg_for_this_sample, total_mean_coverage_of_allele_above_threshold,
                                               total_mean_coverage_of_all_other_alleles_above_threshold, genotyping_options->get_error_rate(), gaps);

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


boost::optional<SampleInfo::IndexAndConfidenceAndMaxLikelihood> SampleInfo::get_confidence () const {
    std::vector<double> likelihoods_for_all_alleles = get_likelihoods_for_all_alleles();

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
    bool enough_total_covg = (total_mean_coverage_of_max_likelihood_allele + total_mean_coverage_of_second_max_likelihood_allele >= genotyping_options->get_min_site_total_covg());
    bool enough_difference_in_covg = (std::abs((int64_t)(total_mean_coverage_of_max_likelihood_allele) -
                                               (int64_t)(total_mean_coverage_of_second_max_likelihood_allele))
                                      >= genotyping_options->get_min_site_diff_covg());

    if (enough_total_covg and enough_difference_in_covg)
        return std::make_tuple(index_of_max_likelihood, confidence, max_likelihood);
    else
        return boost::none;
}


std::string SampleInfo::get_confidence_to_string () const {
    std::stringstream ss;
    auto index_and_confidence_and_max_likelihood_optional = get_confidence();

    if (index_and_confidence_and_max_likelihood_optional) {
        ss << std::get<1>(*index_and_confidence_and_max_likelihood_optional);
    } else {
        ss << ".";
    }

    return ss.str();
}

boost::optional<SampleInfo::GenotypeAndMaxLikelihood> SampleInfo::get_genotype_from_coverage () {
    auto index_and_confidence_and_max_likelihood_optional = get_confidence();

    if (index_and_confidence_and_max_likelihood_optional) {
        size_t index_of_max_likelihood;
        double confidence;
        double max_likelihood;
        std::tie(index_of_max_likelihood, confidence, max_likelihood) = *index_and_confidence_and_max_likelihood_optional;

        bool satisfy_confidence_threshold = confidence  > genotyping_options->get_confidence_threshold();
        if (satisfy_confidence_threshold) {
            return std::make_pair((uint32_t)index_of_max_likelihood, max_likelihood);
        }
    }
    return boost::none;
}


std::string SampleInfo::to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const {
    bool only_one_flag_is_set = ((int)(genotyping_from_maximum_likelihood) + (int)(genotyping_from_coverage)) == 1;
    assert(only_one_flag_is_set);

    std::stringstream out;
    for (uint32_t allele = 0; allele < get_number_of_alleles(); ++allele) {
        if (genotyping_from_maximum_likelihood) {
            out << gt_from_max_likelihood_path_to_string();
        }
        if (genotyping_from_coverage) {
            out << gt_from_coverages_compatible_to_string();
        }

        out << ":" << get_mean_forward_coverage(allele)
            << ":" << get_mean_reverse_coverage(allele)
            << ":" << get_median_forward_coverage(allele)
            << ":" << get_median_reverse_coverage(allele)
            << ":" << get_sum_forward_coverage(allele)
            << ":" << get_sum_reverse_coverage(allele)
            << ":" << get_gaps(allele);

        if (genotyping_from_coverage) {
            out << ":" << get_likelihoods_for_all_alleles()[allele]
                << ":" << get_confidence_to_string();
        }
    }

    return out.str();
}


void SamplesInfos::merge_other_samples_infos_into_this(const SamplesInfos &other, uint32_t allele_offset) {
    bool samples_infos_to_be_merged_have_the_same_size = size() == other.size();
    assert(samples_infos_to_be_merged_have_the_same_size);

    for (size_t sample_index = 0; sample_index < size(); ++sample_index) {
        (*this)[sample_index].merge_other_sample_info_into_this(other[sample_index], allele_offset);
    }
}


std::string SamplesInfos::to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const {
    bool only_one_flag_is_set = ((int)(genotyping_from_maximum_likelihood) + (int)(genotyping_from_coverage)) == 1;
    assert(only_one_flag_is_set);

    std::stringstream out;

    for (uint32_t sample_info_index = 0; sample_info_index < size(); ++sample_info_index) {
        const SampleInfo &sample_info = (*this)[sample_info_index];
        out << sample_info.to_string(genotyping_from_maximum_likelihood, genotyping_from_coverage);

        bool is_the_last_sample_info = sample_info_index == size()-1;
        if (not is_the_last_sample_info)
            out << "\t";
    }

    return out.str();
}


void SampleInfo::merge_other_sample_info_into_this (const SampleInfo &other, uint32_t allele_offset) {
    allele_to_forward_coverages.insert(allele_to_forward_coverages.end(),
                                       other.allele_to_forward_coverages.begin(), other.allele_to_forward_coverages.end());
    allele_to_reverse_coverages.insert(allele_to_reverse_coverages.end(),
                                       other.allele_to_reverse_coverages.begin(), other.allele_to_reverse_coverages.end());

    merge_other_sample_gt_from_max_likelihood_path_into_this(other, allele_offset);
}

void SampleInfo::merge_other_sample_gt_from_max_likelihood_path_into_this (const SampleInfo &other, uint32_t allele_offset) {
    if (not other.is_gt_from_max_likelihood_path_valid())
        return;

    if (not this->is_gt_from_max_likelihood_path_valid()) {
        if (other.get_gt_from_max_likelihood_path() == 0) {
            this->set_gt_from_max_likelihood_path(0);
        } else {
            uint32_t new_gt_from_max_likelihood_path = other.get_gt_from_max_likelihood_path() + allele_offset;
            this->set_gt_from_max_likelihood_path(new_gt_from_max_likelihood_path);
        }
    } else if (this->get_gt_from_max_likelihood_path() != 0 or other.get_gt_from_max_likelihood_path() != 0) {
        //conflict, genotype using the coverages to solve
        this->genotype_from_coverage();
        this->set_gt_from_max_likelihood_path(this->get_gt_coverages());
    }
}