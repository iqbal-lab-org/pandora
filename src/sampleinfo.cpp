#include "sampleinfo.h"

void SampleInfo::add_coverage_information (const std::vector< std::vector<uint32_t> > &allele_to_forward_coverages,
                                           const std::vector< std::vector<uint32_t> > &allele_to_reverse_coverages) {
    this->allele_to_forward_coverages = allele_to_forward_coverages;
    this->allele_to_reverse_coverages = allele_to_reverse_coverages;
    assert(check_if_coverage_information_is_correct());
}

void SampleInfo::genotype_from_coverage () {
    assert(check_if_coverage_information_is_correct());

    auto genotype_and_max_likelihood_optional = get_genotype_from_coverage();
    if (genotype_and_max_likelihood_optional) {
        std::tie(GT_from_coverages, likelihood_of_GT_from_coverages) = *genotype_and_max_likelihood_optional;
    }else {
        GT_from_coverages = boost::none;
        likelihood_of_GT_from_coverages = boost::none;
    }

    // TODO: I don't really like this side-effect - refactor this
    set_gt_from_coverages_compatible(GT_from_coverages);
}


bool SampleInfo::check_if_coverage_information_is_correct() const {
    // at least two alleles because a VCF record has a ref and at least one alt
    bool there_are_at_least_two_alleles = allele_to_forward_coverages.size() >= 2 and allele_to_reverse_coverages.size() >= 2;
    bool forward_and_reverse_coverages_have_the_same_number_of_alleles = allele_to_forward_coverages.size() == allele_to_reverse_coverages.size();

    bool all_alleles_in_forward_and_reverse_have_the_same_number_of_bases = true;
    for (size_t allele_index = 0; allele_index < allele_to_forward_coverages.size(); ++allele_index) {
        bool alleles_in_forward_and_reverse_have_the_same_number_of_bases =
                allele_to_forward_coverages[allele_index].size() == allele_to_reverse_coverages[allele_index].size();
        if (not alleles_in_forward_and_reverse_have_the_same_number_of_bases) {
            all_alleles_in_forward_and_reverse_have_the_same_number_of_bases = false;
            break;
        }
    }

    return there_are_at_least_two_alleles and forward_and_reverse_coverages_have_the_same_number_of_alleles and all_alleles_in_forward_and_reverse_have_the_same_number_of_bases;
}

double SampleInfo::get_gaps (uint32_t allele) const {
    double gaps = 0.0;
    size_t number_of_bases_in_allele = allele_to_forward_coverages[allele].size();

    for (size_t base_index = 0; base_index < number_of_bases_in_allele; ++base_index) {
        if (allele_to_forward_coverages[allele][base_index] + allele_to_reverse_coverages[allele][base_index] < genotyping_options->get_min_kmer_covg())
            gaps++;
    }

    return gaps / number_of_bases_in_allele;
}

std::vector<double> SampleInfo::get_likelihoods_for_all_alleles () const {
    std::vector<double> likelihoods;

    uint32_t min_coverage_threshold = get_min_coverage_threshold_for_this_sample();
    uint32_t total_mean_coverage_over_all_alleles_above_threshold = get_total_mean_coverage_over_all_alleles_given_a_minimum_threshold(
            min_coverage_threshold);

    for (size_t allele = 0; allele < get_number_of_alleles(); ++allele) {
        uint32_t total_mean_coverage_of_allele_above_threshold = get_total_mean_coverage_given_a_minimum_threshold(
                allele, min_coverage_threshold);
        bool min_coverage_threshold_is_satisfied = total_mean_coverage_of_allele_above_threshold > 0;

        uint32_t total_mean_coverage_of_all_other_alleles_above_threshold = total_mean_coverage_over_all_alleles_above_threshold - total_mean_coverage_of_allele_above_threshold;

        double gaps = get_gaps(allele);

        double likelihood = compute_likelihood(min_coverage_threshold_is_satisfied, exp_depth_covg_for_this_sample, total_mean_coverage_of_allele_above_threshold,
                                               total_mean_coverage_of_all_other_alleles_above_threshold, genotyping_options->get_error_rate(), gaps);

        likelihoods.push_back(likelihood);
    }

    return likelihoods;
}


uint32_t SampleInfo::get_total_mean_coverage_given_a_minimum_threshold(uint32_t allele, uint32_t minimum_threshold) const {
    uint32_t total_mean_coverage = this->get_mean_coverage_both_alleles(allele);
    if (total_mean_coverage >= minimum_threshold)
        return total_mean_coverage;
    else
        return 0;
}


uint32_t SampleInfo::get_total_mean_coverage_over_all_alleles_given_a_minimum_threshold(uint32_t minimum_threshold) const {
    uint32_t total_mean_coverage_over_all_alleles = 0;
    for (uint32_t allele = 0; allele < get_number_of_alleles(); ++allele)
        total_mean_coverage_over_all_alleles += get_total_mean_coverage_given_a_minimum_threshold(allele,
                                                                                                  minimum_threshold);
    return total_mean_coverage_over_all_alleles;
}


double SampleInfo::compute_likelihood(bool min_coverage_threshold_is_satisfied, double expected_depth_covg, double total_mean_coverage_of_allele_above_threshold,
                                 double total_mean_coverage_of_all_other_alleles_above_threshold, double error_rate,
                                 double gaps) const {
    // For more details about this, look at page 66 of the thesis
    if (min_coverage_threshold_is_satisfied) {
        return       - expected_depth_covg
                     + total_mean_coverage_of_allele_above_threshold * log(expected_depth_covg)
                     - Maths::logfactorial(total_mean_coverage_of_allele_above_threshold)
                     + total_mean_coverage_of_all_other_alleles_above_threshold * log(error_rate)
                     - expected_depth_covg * gaps
                     + log(1 - exp(-expected_depth_covg)) * (1 - gaps);
    } else {
        return       - expected_depth_covg
                     + total_mean_coverage_of_all_other_alleles_above_threshold * log(error_rate)
                     - expected_depth_covg * gaps
                     + log(1 - exp(-expected_depth_covg)) * (1 - gaps);
    }
}


boost::optional<SampleInfo::IndexAndConfidenceAndMaxLikelihood> SampleInfo::get_confidence () const {
    std::vector<double> likelihoods_for_all_alleles = get_likelihoods_for_all_alleles();

    size_t index_of_max_likelihood = Maths::arg_max(likelihoods_for_all_alleles.begin(), likelihoods_for_all_alleles.end());
    size_t index_of_second_max_likelihood; //guaranteed to exist, since we have at least two alleles
    {
        std::vector<double> likelihoods_for_all_alleles_without_max_element(likelihoods_for_all_alleles.begin(), likelihoods_for_all_alleles.end());
        likelihoods_for_all_alleles_without_max_element.erase(likelihoods_for_all_alleles_without_max_element.begin()+index_of_max_likelihood);
        index_of_second_max_likelihood = Maths::arg_max(likelihoods_for_all_alleles_without_max_element.begin(), likelihoods_for_all_alleles_without_max_element.end());

        bool index_of_second_max_likelihood_is_past_index_of_max_likelihood = index_of_second_max_likelihood >= index_of_max_likelihood;
        if (index_of_second_max_likelihood_is_past_index_of_max_likelihood) {
            ++index_of_second_max_likelihood;
        }
    }

    double max_likelihood = likelihoods_for_all_alleles[index_of_max_likelihood];
    double second_max_likelihood = likelihoods_for_all_alleles[index_of_second_max_likelihood];
    double confidence = std::abs(max_likelihood - second_max_likelihood);

    uint32_t total_mean_coverage_of_max_likelihood_allele = get_mean_coverage_both_alleles(index_of_max_likelihood);
    uint32_t total_mean_coverage_of_second_max_likelihood_allele = get_mean_coverage_both_alleles(
            index_of_second_max_likelihood);
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

boost::optional<SampleInfo::GenotypeAndMaxLikelihood> SampleInfo::get_genotype_from_coverage () const {
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


std::string SampleInfo::to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_compatible_coverage) const {
    bool only_one_flag_is_set = ((int)(genotyping_from_maximum_likelihood) + (int)(genotyping_from_compatible_coverage)) == 1;
    assert(only_one_flag_is_set);

    std::vector<double> likelihoods_for_all_alleles = get_likelihoods_for_all_alleles();

    std::stringstream out;
    if (genotyping_from_maximum_likelihood) {
        out << gt_from_max_likelihood_path_to_string();
    }
    if (genotyping_from_compatible_coverage) {
        out << gt_from_coverages_compatible_to_string();
    }

    auto uint32_methods_to_call = {&SampleInfo::get_mean_forward_coverage, &SampleInfo::get_mean_reverse_coverage,
                                   &SampleInfo::get_median_forward_coverage, &SampleInfo::get_median_reverse_coverage,
                                   &SampleInfo::get_sum_forward_coverage, &SampleInfo::get_sum_reverse_coverage};
    for (const auto &method : uint32_methods_to_call) {
        print_the_output_of_method_for_each_allele(out, this, method);
    }

    auto double_methods_to_call = {&SampleInfo::get_gaps};
    for (const auto &method : double_methods_to_call) {
        print_the_output_of_method_for_each_allele(out, this, method);
    }

    if (genotyping_from_compatible_coverage) {
        out << ":";
        for (uint32_t allele = 0; allele < get_number_of_alleles(); ++allele) {
            if (allele > 0) {
                out << ",";
            }
            out << likelihoods_for_all_alleles[allele];
        }

        out << ":" << get_confidence_to_string();
    }

    return out.str();
}


void SampleInfo::merge_other_sample_info_into_this (const SampleInfo &other) {
    uint32_t allele_offset = this->get_number_of_alleles();

    allele_to_forward_coverages.insert(allele_to_forward_coverages.end(),
                                       other.allele_to_forward_coverages.begin()+1, other.allele_to_forward_coverages.end());
    allele_to_reverse_coverages.insert(allele_to_reverse_coverages.end(),
                                       other.allele_to_reverse_coverages.begin()+1, other.allele_to_reverse_coverages.end());

    merge_other_sample_gt_from_max_likelihood_path_into_this(other, allele_offset);

    //TODO: We do not merge GT_from_coverages_compatible as it is not needed
}

void SampleInfo::merge_other_sample_gt_from_max_likelihood_path_into_this (const SampleInfo &other, uint32_t allele_offset) {
    if (not other.is_gt_from_max_likelihood_path_valid())
        return;

    if (not this->is_gt_from_max_likelihood_path_valid()) {
        if (other.get_gt_from_max_likelihood_path() == 0) {
            this->set_gt_from_max_likelihood_path(0);
        } else {
            uint32_t new_gt_from_max_likelihood_path = other.get_gt_from_max_likelihood_path() + allele_offset - 1;
            this->set_gt_from_max_likelihood_path(new_gt_from_max_likelihood_path);
        }
    } else if (this->get_gt_from_max_likelihood_path() != 0 or other.get_gt_from_max_likelihood_path() != 0) {
        //conflict, genotype using the coverages to solve
        boost::optional<SampleInfo::GenotypeAndMaxLikelihood> genotype_and_max_likelihood_from_coverage_optional = this->get_genotype_from_coverage();
        boost::optional<uint32_t> genotype_from_coverage = boost::none;
        if (genotype_and_max_likelihood_from_coverage_optional) {
            genotype_from_coverage = genotype_and_max_likelihood_from_coverage_optional->first;
        }
        this->set_gt_from_max_likelihood_path(genotype_from_coverage);
    }
}


void SampleInfo::solve_incompatible_gt_conflict_with (SampleInfo &other) {
    bool any_of_gts_are_invalid_thus_no_conflict = not this->is_gt_from_coverages_compatible_valid() or not other.is_gt_from_coverages_compatible_valid();
    if (any_of_gts_are_invalid_thus_no_conflict)
        return;

    bool both_gts_are_to_ref_thus_no_conflict =
            this->get_gt_from_coverages_compatible() == 0 and other.get_gt_from_coverages_compatible() == 0;
    if (both_gts_are_to_ref_thus_no_conflict)
        return;

    if (this->get_likelihood_of_gt_from_coverages_compatible() > other.get_likelihood_of_gt_from_coverages_compatible()) {
            if (this->get_gt_from_coverages_compatible() == 0) {
                other.set_gt_from_coverages_compatible(0);
            }
            else {
                other.set_gt_from_coverages_compatible(boost::none);
            }
    } else {
        if (other.get_gt_from_coverages_compatible() == 0) {
            this->set_gt_from_coverages_compatible(0);
        }
        else {
            this->set_gt_from_coverages_compatible(boost::none);
        }
    }
}
