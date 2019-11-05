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
#include <boost/bind.hpp>

// TODO: use memoization to speed up everything here
// TODO: this class is doing too much. There is the concept of an allele info which can be factored out to another class
// TODO: also there is SampleInfo hierarchy hidden here, where two subclasses could be derived, one dealing with genotyping from max likelihood path
// TODO: and the other with genotyping from coverage, and maybe a third one to deal with compatible genotypes
class SampleInfo {
public:
    SampleInfo(uint32_t sample_index, GenotypingOptions const *  genotyping_options) :
    sample_index(sample_index),
    genotyping_options(genotyping_options),
    exp_depth_covg_for_this_sample(genotyping_options->get_sample_index_to_exp_depth_covg()[sample_index])
    {}

    virtual ~SampleInfo(){}



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // methods related to genotyping from the max likelihood path
    virtual inline void set_gt_from_max_likelihood_path (const boost::optional<uint32_t> &gt) {
        this->GT_from_maximum_likelihood_path = gt;
    }

    virtual inline bool is_gt_from_max_likelihood_path_valid () const {
        return (bool)(this->GT_from_maximum_likelihood_path);
    }

    virtual inline uint32_t get_gt_from_max_likelihood_path () const {
        if (not is_gt_from_max_likelihood_path_valid()) {
            throw std::runtime_error("GT_from_maximum_likelihood_path is not valid at SampleInfo::get_gt_from_max_likelihood_path_valid()");
        }
        return *(this->GT_from_maximum_likelihood_path);
    }

    virtual inline std::string gt_from_max_likelihood_path_to_string() const {
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
    virtual void add_coverage_information (const std::vector< std::vector<uint32_t> > &allele_to_forward_coverages,
                                   const std::vector< std::vector<uint32_t> > &allele_to_reverse_coverages);

    virtual inline const std::vector< std::vector<uint32_t> > & get_allele_to_forward_coverages() const {
        return allele_to_forward_coverages;
    }

    virtual inline const std::vector< std::vector<uint32_t> > & get_allele_to_reverse_coverages() const {
        return allele_to_reverse_coverages;
    }

    virtual void genotype_from_coverage ();

    virtual inline bool is_gt_from_coverages_valid () const {
        return (bool)(this->GT_from_coverages);
    }

    virtual inline uint32_t get_gt_from_coverages () const {
        if (not is_gt_from_coverages_valid()) {
            throw std::runtime_error("GT_from_coverages is not valid at SampleInfo::get_gt_from_coverages()");
        }
        return *(this->GT_from_coverages);
    }

    virtual inline double get_likelihood_of_gt_from_coverages () const {
        if (not is_gt_from_coverages_valid()) {
            throw std::runtime_error("GT_from_coverages is not valid at SampleInfo::get_gt_from_coverages()");
        }
        return *(this->likelihood_of_GT_from_coverages);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // methods related to genotyping from the coverage - compatible
    virtual inline bool is_gt_from_coverages_compatible_valid () const {
        return (bool)(this->GT_from_coverages_compatible);
    }

    virtual inline void set_gt_coverages_compatible (const boost::optional<uint32_t> &gt) {
        this->GT_from_coverages_compatible = gt;
    }

    virtual inline uint32_t get_gt_coverages_compatible () const {
        if (not is_gt_from_coverages_compatible_valid()) {
            throw std::runtime_error("GT_from_coverages_compatible is not valid at SampleInfo::get_gt_coverages_compatible()");
        }
        return *(this->GT_from_coverages_compatible);
    }

    virtual inline std::string gt_from_coverages_compatible_to_string() const {
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
    virtual inline void merge_other_sample_info_into_this (const SampleInfo &other);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // trivial getters over format fields
    virtual inline size_t get_number_of_alleles() const {
        return allele_to_forward_coverages.size();
    }

    virtual inline uint32_t get_mean_forward_coverage(uint32_t allele) const {
        return Maths::mean(allele_to_forward_coverages[allele].begin(), allele_to_forward_coverages[allele].end());
    }

    virtual inline uint32_t get_median_forward_coverage(uint32_t allele) const {
        return Maths::median(allele_to_forward_coverages[allele].begin(), allele_to_forward_coverages[allele].end());
    }

    virtual inline uint32_t get_sum_forward_coverage(uint32_t allele) const {
        return Maths::sum(allele_to_forward_coverages[allele].begin(), allele_to_forward_coverages[allele].end());
    }

    virtual inline uint32_t get_mean_reverse_coverage(uint32_t allele) const {
        return Maths::mean(allele_to_reverse_coverages[allele].begin(), allele_to_reverse_coverages[allele].end());
    }

    virtual inline uint32_t get_median_reverse_coverage(uint32_t allele) const {
        return Maths::median(allele_to_reverse_coverages[allele].begin(), allele_to_reverse_coverages[allele].end());
    }

    virtual inline uint32_t get_sum_reverse_coverage(uint32_t allele) const {
        return Maths::sum(allele_to_reverse_coverages[allele].begin(), allele_to_reverse_coverages[allele].end());
    }

    virtual inline uint32_t get_mean_coverage_both_alleles (uint32_t allele) const {
        return this->get_mean_forward_coverage(allele) + this->get_mean_reverse_coverage(allele);
    }

    virtual inline uint32_t get_min_coverage_threshold_for_this_sample() const {
        return std::max(genotyping_options->get_min_allele_covg(),
                uint(genotyping_options->get_min_fraction_allele_covg() * exp_depth_covg_for_this_sample));
    }


    //other trivial getter
    virtual uint32_t get_sample_index() const {
        return sample_index;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // non-trivial getters over format fields
    virtual double get_gaps (uint32_t allele) const;

    virtual std::vector<double> get_likelihoods_for_all_alleles () const;


    using IndexAndConfidenceAndMaxLikelihood = std::tuple<size_t, double, double>;
    virtual boost::optional<IndexAndConfidenceAndMaxLikelihood> get_confidence () const;
    virtual std::string get_confidence_to_string () const;

    using GenotypeAndMaxLikelihood = std::pair<uint32_t, double>;
    virtual boost::optional<GenotypeAndMaxLikelihood> get_genotype_from_coverage () const;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    virtual std::string to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_compatible_coverage) const;
protected:
    uint32_t sample_index;
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


    virtual bool check_if_coverage_information_is_correct() const;


    /////////////////////////////////////////////
    // get_likelihoods_for_all_alleles() helpers
    virtual uint32_t get_total_mean_coverage_given_a_minimum_threshold(uint32_t allele, uint32_t minimum_threshold) const;
    virtual uint32_t get_total_mean_coverage_over_all_alleles_given_a_minimum_threshold(uint32_t minimum_threshold) const;
    virtual double compute_likelihood(bool min_coverage_threshold_is_satisfied, double expected_depth_covg, double total_mean_coverage_of_allele_above_threshold,
                                      double total_mean_coverage_of_all_other_alleles_above_threshold, double error_rate,
                                                 double gaps) const;
    /////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // merge_other_sample_info_into_this() helpers
    virtual void merge_other_sample_gt_from_max_likelihood_path_into_this(const SampleInfo &other, uint32_t allele_offset);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // to_string() helpers - TODO: this could be simplified with generic lambdas, but these are just available in C++14
    template <class METHOD_TYPE>
    static void print_the_output_of_method_for_each_allele (std::stringstream &out, const SampleInfo* const sample_info,
            const METHOD_TYPE &method) {
        out << ":";

        for (uint32_t allele = 0; allele < sample_info->get_number_of_alleles(); ++allele) {
            if (allele > 0) {
                out << ",";
            }
            out << boost::bind(method, sample_info, _1)(allele);
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};

//NB: we use template here so that we can change what the container holds
//This allows us to hold whatever thing (e.g. a mock) to be able to test this class
template<class SAMPLE_TYPE>
class SampleIndexToSampleInfoTemplate : public std::vector<SAMPLE_TYPE> {
public:
    SampleIndexToSampleInfoTemplate(){}
    virtual ~SampleIndexToSampleInfoTemplate(){}

    virtual inline void emplace_back_several_empty_sample_infos (size_t amount, GenotypingOptions const * genotyping_options) {
        size_t initial_size = this->size();
        for (size_t index = initial_size; index < initial_size + amount; ++index) {
            this->emplace_back(index, genotyping_options);
        }
    }


    virtual inline void merge_other_samples_infos_into_this(const SampleIndexToSampleInfoTemplate<SAMPLE_TYPE> &other) {
        bool same_number_of_samples = this->size() == other.size();
        assert(same_number_of_samples);

        for (size_t sample_index = 0; sample_index < this->size(); ++sample_index) {
            (*this)[sample_index].merge_other_sample_info_into_this(other[sample_index]);
        }
    }


    virtual std::string to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage) const {
        std::stringstream out;

        for (uint32_t sample_info_index = 0; sample_info_index < this->size(); ++sample_info_index) {
            const SAMPLE_TYPE &sample_info = (*this)[sample_info_index];
            out << sample_info.to_string(genotyping_from_maximum_likelihood, genotyping_from_coverage);

            bool is_the_last_sample_info = sample_info_index == this->size()-1;
            if (not is_the_last_sample_info)
                out << "\t";
        }

        return out.str();
    }

    virtual inline void genotype_from_coverage () {
        for (SAMPLE_TYPE &sample_info : (*this)) {
            sample_info.genotype_from_coverage();
        }
    }
};

// the SampleIndexToSampleInfoTemplate that is used in production
using SampleIndexToSampleInfo = SampleIndexToSampleInfoTemplate<SampleInfo>;

#endif //PANDORA_SAMPLEINFO_H
