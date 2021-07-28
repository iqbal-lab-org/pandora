#ifndef __KMERGRAPHWITHCOVERAGE_H_INCLUDED__ // if kmergraphwithcoverage.h hasn't been
                                             // included yet...
#define __KMERGRAPHWITHCOVERAGE_H_INCLUDED__

class LocalPRG;

#include <cstdint>
#include <vector>
#include <iostream>
#include "prg/path.h"
#include "kmernode.h"
#include "kmergraph.h"
#include "pangenome/ns.cpp"
#include "utils.h"
#include "fatal_error.h"

/**
 * Represents an annotated KmerGraph, where the annotation is only the coverage on the
 * nodes This is more lightweight than copying a whole k-mer graph and changing the
 * coverage of each node
 * TODO: generalize it to general annotations?
 * NOTE: representing as composition instead of inheritance because one KmerGraph can be
 * referenced by several KmerGraphWithCoverage, thus KmerGraph's life cycle does not
 * depend on KmerGraphWithCoverage's life cycle
 */
class KmerGraphWithCoverage {
private:
    // each coverage is represented as a std::vector, each position of the vector
    // representing the coverage of a sample each coverage of a sample is represented by
    // a std::pair, for the forward and reverse reads coverage of that sample
    std::vector<std::vector<std::pair<uint16_t, uint16_t>>>
        node_index_to_sample_coverage;
    uint32_t exp_depth_covg;
    float binomial_parameter_p;
    float negative_binomial_parameter_p;
    float negative_binomial_parameter_r;
    int thresh;
    uint32_t total_number_samples;
    uint32_t num_reads;
    uint32_t get_covg(
        uint32_t node_id, pandora::Strand strand, uint32_t sample_id) const;
    void increment_covg(uint32_t node_id, pandora::Strand strand, uint32_t sample_id);
    void set_covg(
        uint32_t node_id, uint16_t value, pandora::Strand strand, uint32_t sample_id);

public:
    KmerGraph* kmer_prg; // the underlying KmerGraph - TODO: it is dangerous to leave
                         // this public, make it private? - this should be const

    // constructor, destructors, etc
    KmerGraphWithCoverage(KmerGraph* kmer_prg, uint32_t total_number_samples = 1)
        : node_index_to_sample_coverage(kmer_prg->nodes.size())
        , exp_depth_covg { 0 }
        , binomial_parameter_p { 1 }
        , negative_binomial_parameter_p { 0.015 }
        , negative_binomial_parameter_r { 2 }
        , thresh { -25 }
        , total_number_samples { total_number_samples }
        , num_reads { 0 }
        , kmer_prg { kmer_prg }
    {
        const bool kmer_prg_is_invalid = kmer_prg == nullptr;
        if (kmer_prg_is_invalid) {
            fatal_error("Error building Kmer Graph With Coverage: kmer PRG is invalid");
        }
        zeroCoverages();
    }
    KmerGraphWithCoverage(const KmerGraphWithCoverage& other)
        = default; // copy default constructor
    KmerGraphWithCoverage(KmerGraphWithCoverage&& other)
        = default; // move default constructor
    KmerGraphWithCoverage& operator=(const KmerGraphWithCoverage& other)
        = default; // copy assignment operator
    KmerGraphWithCoverage& operator=(KmerGraphWithCoverage&& other)
        = default; // move assignment operator
    virtual ~KmerGraphWithCoverage() = default;

    // getters
    /*
    TODO: we should return a uint16_t instead of a uint32_t, but I did not want to
    change the interface, as it can be dangerous WARNING: some parts of the code
    accumulates under the type of this return value, so it might be dangerous to change
    to uint16_t without careful analysis. e.g.: I can see summing a bunch of coverage
    and potentially getting a value larger than 65535, which would incur overflow and
    bugs. NOTE: converting uint16_t to uint32_t is safe anyway, so it is fine to leave
    like this NOTE: the advantage of returning a uint16_t here is a negligible
    improvement in RAM usage, so I don't think it is worth the risk now
    TODO: leaving return type as uint32_t to be safe for now, need a recheck later
     */
    uint32_t get_forward_covg(uint32_t node_id, uint32_t sample_id) const
    {
        return this->get_covg(node_id, pandora::Strand::Forward, sample_id);
    }
    uint32_t get_reverse_covg(uint32_t node_id, uint32_t sample_id) const
    {
        return this->get_covg(node_id, pandora::Strand::Reverse, sample_id);
    }
    uint32_t get_num_reads() const { return num_reads; }
    uint32_t get_total_number_samples() const { return total_number_samples; }

    // setters
    void increment_forward_covg(uint32_t node_id, uint32_t sample_id)
    {
        this->increment_covg(node_id, pandora::Strand::Forward, sample_id);
    }
    void increment_reverse_covg(uint32_t node_id, uint32_t sample_id)
    {
        this->increment_covg(node_id, pandora::Strand::Reverse, sample_id);
    }
    void set_forward_covg(uint32_t node_id, uint16_t value, uint32_t sample_id)
    {
        this->set_covg(node_id, value, pandora::Strand::Forward, sample_id);
    }
    void set_reverse_covg(uint32_t node_id, uint16_t value, uint32_t sample_id)
    {
        this->set_covg(node_id, value, pandora::Strand::Reverse, sample_id);
    }
    void set_exp_depth_covg(uint32_t);
    void set_binomial_parameter_p(float);
    void set_negative_binomial_parameters(const float&, const float&);
    void set_thresh(int thresh) { this->thresh = thresh; }
    void set_num_reads(uint32_t num_reads) { this->num_reads = num_reads; }

    void zeroCoverages()
    {
        for (auto& sampleCoverage : node_index_to_sample_coverage) {
            sampleCoverage
                = std::vector<std::pair<uint16_t, uint16_t>>(total_number_samples);
            sampleCoverage.shrink_to_fit(); // tries to make this information as compact
                                            // as possible (TODO: use sdsl?)
        }
    }

    float nbin_prob(uint32_t, const uint32_t& sample_id);

    float lin_prob(uint32_t, const uint32_t& sample_id);

    float bin_prob(uint32_t, const uint32_t& sample_id);

    float bin_prob(const uint32_t&, const uint32_t&, const uint32_t& sample_id);

    float get_prob(const std::string& prob_model, const uint32_t& node_id,
        const uint32_t& sample_id);

    bool coverage_is_zeroes(const uint32_t&);

    float find_max_path(std::vector<KmerNodePtr>& maxpath,
        const std::string& prob_model, const uint32_t& max_num_kmers_to_average,
        const uint32_t& sample_id, const pangenome::Node *pangenome_node = NULL);

    std::vector<std::vector<KmerNodePtr>> find_max_paths(
        uint32_t, const uint32_t& sample_id);

    void save_covg_dist(const std::string&);

    std::vector<std::vector<KmerNodePtr>> get_random_paths(uint32_t);

    float prob_path(const std::vector<KmerNodePtr>& kpath, const uint32_t& sample_id,
        const std::string& prob_model);

    float prob_paths(const std::vector<std::vector<KmerNodePtr>>&);

    void save(
        const fs::path& filepath, std::shared_ptr<LocalPRG> localprg = nullptr) const;
    void load(const std::string&);

    // test friends
    friend class KmerGraphWithCoverageTest_set_exp_depth_covg_Test;
    friend class KmerGraphWithCoverageTest_set_p_Test;
    friend class KmerGraphWithCoverageTest_set_nb_Test;
    friend class KmerGraphWithCoverageTest_prob_failNoNumReads_Test;
    friend class KmerGraphWithCoverageTest_prob_simple_Test;
    friend class KmerGraphWithCoverageTest_prob_realNodeCovgs_Test;
    friend class KmerGraphWithCoverageTest_findMaxPath_InvalidProbModel_Test;
    friend class KmerGraphWithCoverageTest_findMaxPathSimple_Test;
    friend class KmerGraphWithCoverageTest_findMaxPathSimple_WithMaxKmersInAvg_Test;
    friend class KmerGraphWithCoverageTest_findMaxPath2Level_bin_Test;
    friend class KmerGraphWithCoverageTest_findMaxPath2Level_nbin_Test;
    friend class KmerGraphWithCoverageTest_findMaxPath2Level_lin_Test;

    friend class KmerGraphWithCoverageTest_find_max_paths_2Level_Test;
    friend class KmerGraphWithCoverageTest_path_probs_Test;
};

#endif
