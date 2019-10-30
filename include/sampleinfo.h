#ifndef PANDORA_SAMPLEINFO_H
#define PANDORA_SAMPLEINFO_H

#include <unordered_map>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <boost/optional.hpp>
#include "Maths.h"
#include "OptionsAggregator.h"


// TODO: use memoization to speed up everything here
// TODO: this class is doing too much. There is the concept of an allele info which can be factored out to another class
// TODO: also there is SampleInfo hierarchy hidden here, where two subclasses could be derived, one dealing with genotyping from max likelihood path
// TODO: and the other with genotyping from coverage
class SampleInfo {
public:
    SampleInfo(uint32_t sample_index, GenotypingOptions const *  genotyping_options) :
    genotyping_options(genotyping_options),
    exp_depth_covg_for_this_sample(genotyping_options->get_sample_index_to_exp_depth_covg()[sample_index])
    {}



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // methods related to genotyping from the max likelihood path
    inline void set_gt_from_max_likelihood_path (const boost::optional<uint32_t> &gt) {
        this->GT_from_maximum_likelihood_path = gt;
    }

    inline bool is_gt_from_max_likelihood_path_valid () const {
        return (bool)(this->GT_from_maximum_likelihood_path);
    }

    inline uint32_t get_gt_from_max_likelihood_path () const {
        if (not is_gt_from_max_likelihood_path_valid()) {
            throw std::runtime_error("GT_from_maximum_likelihood_path is not valid at SampleInfo::get_gt_from_max_likelihood_path_valid()");
        }
        return *(this->GT_from_maximum_likelihood_path);
    }

    inline std::string gt_from_max_likelihood_path_to_string() const {
        std::stringstream ss;
        if (is_gt_from_max_likelihood_path_valid())
            ss << get_gt_from_max_likelihood_path();
        else
            ss << ".";
        return ss.str();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // methods related to genotyping from the coverage
    void add_coverage_information (const std::vector< std::vector<uint32_t> > &allele_to_forward_coverages,
                                   const std::vector< std::vector<uint32_t> > &allele_to_reverse_coverages);

    void genotype_from_coverage ();

    inline bool is_gt_from_coverages_valid () const {
        return (bool)(this->GT_from_coverages);
    }

    inline uint32_t get_gt_coverages () const {
        if (not is_gt_from_coverages_valid()) {
            throw std::runtime_error("GT_from_coverages is not valid at SampleInfo::get_gt_coverages()");
        }
        return *(this->GT_from_coverages);
    }

    inline double get_likelihood_of_GT_from_coverages () const {
        if (not is_gt_from_coverages_valid()) {
            throw std::runtime_error("GT_from_coverages is not valid at SampleInfo::get_gt_coverages()");
        }
        return *(this->likelihood_of_GT_from_coverages);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // methods related to genotyping from the coverage - compatible
    inline bool is_gt_from_coverages_compatible_valid () const {
        return (bool)(this->GT_from_coverages_compatible);
    }

    inline void set_gt_coverages_compatible (const boost::optional<uint32_t> &gt) {
        this->GT_from_coverages_compatible = gt;
    }

    inline uint32_t get_gt_coverages_compatible () const {
        if (not is_gt_from_coverages_compatible_valid()) {
            throw std::runtime_error("GT_from_coverages_compatible is not valid at SampleInfo::get_gt_coverages_compatible()");
        }
        return *(this->GT_from_coverages_compatible);
    }

    inline std::string gt_from_coverages_compatible_to_string() const {
        std::stringstream ss;
        if (is_gt_from_coverages_compatible_valid())
            ss << get_gt_coverages_compatible();
        else
            ss << ".";
        return ss.str();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // merge-related methods
    inline void merge_other_sample_info_into_this (const SampleInfo &other, uint32_t allele_offset);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::string to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const;

private:
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // genotyping options
    GenotypingOptions const * genotyping_options;
    uint32_t exp_depth_covg_for_this_sample;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // A VCF record can have several genotypings in pandora
    // The first comes by comparing the maximum likelihood multisample path vs a maximum likelihood singlesample path
    // The second is the genotyping based on coverage
    // This second genotyping can be not compatible, so then we make it compatible, having thus the third genotyping
    // See https://github.com/iqbal-lab/pandora1_paper/issues/105#issuecomment-547121546 for details

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // genotyping from maximum likelihood path
    boost::optional<uint32_t> GT_from_maximum_likelihood_path;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // genotyping from coverage
    std::vector< std::vector<uint32_t> > allele_to_forward_coverages;
    std::vector< std::vector<uint32_t> > allele_to_reverse_coverages;
    // this is the inferred GT from the coverages
    boost::optional<uint32_t> GT_from_coverages;
    // this is the likelihood of the inferred GT from the coverages, used for making GT compatible later
    boost::optional<double> likelihood_of_GT_from_coverages;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // genotyping from coverage - compatible
    boost::optional<uint32_t> GT_from_coverages_compatible;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // trivial getters over format fields
    inline size_t get_number_of_alleles() const {
        return allele_to_forward_coverages.size();
    }

    inline uint32_t get_mean_forward_coverage(uint32_t allele) const {
        return Maths::mean(allele_to_forward_coverages[allele].begin(), allele_to_forward_coverages[allele].end());
    }

    inline uint32_t get_median_forward_coverage(uint32_t allele) const {
        return Maths::median(allele_to_forward_coverages[allele].begin(), allele_to_forward_coverages[allele].end());
    }

    inline uint32_t get_sum_forward_coverage(uint32_t allele) const {
        return Maths::sum(allele_to_forward_coverages[allele].begin(), allele_to_forward_coverages[allele].end());
    }

    inline uint32_t get_mean_reverse_coverage(uint32_t allele) const {
        return Maths::mean(allele_to_reverse_coverages[allele].begin(), allele_to_reverse_coverages[allele].end());
    }

    inline uint32_t get_median_reverse_coverage(uint32_t allele) const {
        return Maths::median(allele_to_reverse_coverages[allele].begin(), allele_to_reverse_coverages[allele].end());
    }

    inline uint32_t get_sum_reverse_coverage(uint32_t allele) const {
        return Maths::sum(allele_to_reverse_coverages[allele].begin(), allele_to_reverse_coverages[allele].end());
    }

    inline uint32_t get_total_mean_coverage(uint32_t allele) const {
        return this->get_mean_forward_coverage(allele) + this->get_mean_reverse_coverage(allele);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // non-trivial getters over format fields
    double get_gaps (uint32_t allele) const;

    std::vector<double> get_likelihoods_for_all_alleles () const;


    using IndexAndConfidenceAndMaxLikelihood = std::tuple<size_t, double, double>;
    boost::optional<IndexAndConfidenceAndMaxLikelihood> get_confidence () const;
    std::string get_confidence_to_string () const;

    using GenotypeAndMaxLikelihood = std::pair<uint32_t, double>;
    boost::optional<GenotypeAndMaxLikelihood> get_genotype_from_coverage ();

    bool check_if_coverage_information_is_correct() const;

    /////////////////////////////////////////////
    // get_likelihoods_for_all_alleles() helpers
    uint32_t get_total_mean_coverage_given_a_mininum_threshold(uint32_t allele, uint32_t minimum_threshold) const;
    uint32_t get_total_mean_coverage_over_all_alleles_given_a_mininum_threshold(uint32_t minimum_threshold) const;
    static double compute_likelihood(bool min_coverage_threshold_is_satisfied, uint32_t expected_depth_covg, uint32_t total_mean_coverage_of_allele_above_threshold,
                                                 uint32_t total_mean_coverage_of_all_other_alleles_above_threshold, double error_rate,
                                                 double gaps);
    /////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // merge_other_sample_info_into_this() helpers
    void merge_other_sample_gt_from_max_likelihood_path_into_this(const SampleInfo &other, uint32_t allele_offset);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};


class SamplesInfos : public std::vector<SampleInfo> {
public:
    SamplesInfos(){}

    inline void push_back_several_empty_sample_infos (size_t amount, GenotypingOptions const * genotyping_options) {
        for (size_t index = 0; index < amount; ++index)
            push_back(SampleInfo(index, genotyping_options));
    }

    void merge_other_samples_infos_into_this(const SamplesInfos &other, uint32_t allele_offset);

    std::string to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const;

    inline void genotype_from_coverage () {
        for (uint32_t sample_index = 0; sample_index < size(); ++sample_index) {
            (*this)[sample_index].genotype_from_coverage();
        }
    }
};

#endif //PANDORA_SAMPLEINFO_H
