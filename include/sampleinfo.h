#ifndef PANDORA_SAMPLEINFO_H
#define PANDORA_SAMPLEINFO_H

#include <unordered_map>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <boost/optional.hpp>
#include "Maths.h"


// TODO: use memoization to speed up everything here
class SampleInfo {
private:
    std::vector< std::vector<uint32_t> > allele_to_forward_coverages;
    std::vector< std::vector<uint32_t> > allele_to_reverse_coverages;

public:
    SampleInfo(const std::vector< std::vector<uint32_t> > &allele_to_forward_coverages,
            const std::vector< std::vector<uint32_t> > &allele_to_reverse_coverages) :
            allele_to_forward_coverages(allele_to_forward_coverages),
            allele_to_reverse_coverages(allele_to_reverse_coverages) {
        // at least two alleles because a VCF record has a ref and at least one alt
        bool there_are_at_least_two_alleles = allele_to_forward_coverages.size() >= 2 and allele_to_reverse_coverages.size() >= 2;
        bool same_length_vectors = allele_to_forward_coverages.size() == allele_to_reverse_coverages.empty();
        assert(there_are_at_least_two_alleles and same_length_vectors);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // map-like getter
//    template <class RETURN_TYPE>
//    RETURN_TYPE operator[] (const std::string &format_field) const {
//        template <class STORAGE_TYPE>
//        const std::vector<std::string> SampleInfo<STORAGE_TYPE>::keys_that_can_be_merged = {"MEAN_FWD_COVG", "MEAN_REV_COVG",
//                                                                                            "MED_FWD_COVG", "MED_REV_COVG",
//                                                                                            "SUM_FWD_COVG", "SUM_REV_COVG",
//                                                                                            "LIKELIHOOD", "GT_CONF", "GAPS"};
//    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



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

    uint32_t get_total_mean_coverage(uint32_t allele) const {
        return this->get_mean_forward_coverage(allele) + this->get_mean_reverse_coverage(allele);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // non-trivial getters over format fields
    double get_gaps (uint32_t allele, uint32_t threshold_to_consider_a_gap) const;

    std::vector<double> get_likelihoods_for_all_alleles (uint32_t expected_depth_covg, double error_rate, uint32_t min_allele_covg,
                                                         double min_fraction_allele_covg, uint32_t min_kmer_covg) const;


    using IndexAndConfidence = std::pair<size_t, double>;
    boost::optional<IndexAndConfidence> get_confidence (uint32_t expected_depth_covg, double error_rate, uint32_t min_allele_covg,
                                                         double min_fraction_allele_covg, uint32_t min_kmer_covg,
                                                         const uint32_t &min_total_covg, const uint32_t &min_diff_covg) const;

    boost::optional<uint32_t> get_genotype(uint32_t expected_depth_covg, double error_rate, uint32_t min_allele_covg,
                                           double min_fraction_allele_covg, uint32_t min_kmer_covg,
                                           const uint32_t &min_total_covg, const uint32_t &min_diff_covg, const uint16_t confidence_threshold) {
        auto index_and_confidence_optional = get_confidence(expected_depth_covg, error_rate, min_allele_covg,
                min_fraction_allele_covg, min_kmer_covg, min_total_covg, min_diff_covg);

        if (index_and_confidence_optional) {
            size_t index_of_max_likelihood;
            double confidence;
            std::tie(index_of_max_likelihood, confidence) = *index_and_confidence_optional;

            bool satisfy_confidence_threshold = confidence  > confidence_threshold;
            if (satisfy_confidence_threshold) {
                return index_of_max_likelihood;
            }
        }
        return boost::none;
    }

private:
    /////////////////////////////////////////////
    // get_likelihoods_for_all_alleles() helpers
    uint32_t get_total_mean_coverage_given_a_mininum_threshold(uint32_t allele, uint32_t minimum_threshold) const;
    uint32_t get_total_mean_coverage_over_all_alleles_given_a_mininum_threshold(uint32_t minimum_threshold) const;
    static double compute_likelihood(bool min_coverage_threshold_is_satisfied, uint32_t expected_depth_covg, uint32_t total_mean_coverage_of_allele_above_threshold,
                                                 uint32_t total_mean_coverage_of_all_other_alleles_above_threshold, double error_rate,
                                                 double gaps);
    /////////////////////////////////////////////



public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // merge-related methods
    void merge_other_sample_info_into_this (const SampleInfo &other) {
        allele_to_forward_coverages.insert(allele_to_forward_coverages.end(),
                other.allele_to_forward_coverages.begin(), other.allele_to_forward_coverages.end());
        allele_to_reverse_coverages.insert(allele_to_reverse_coverages.end(),
                                           other.allele_to_reverse_coverages.begin(), other.allele_to_reverse_coverages.end());
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};


class SamplesInfos : public std::vector<SampleInfo> {
public:
    inline void push_back_several_empty_sample_infos (size_t amount) {
        for (size_t index = 0; index < amount; ++index)
            this->push_back(SampleInfo({{0,0}, {0,0}}, {{0,0}, {0,0}}));
    }

    void merge_other_samples_infos_into_this(const SamplesInfos &other);
};

#endif //PANDORA_SAMPLEINFO_H
