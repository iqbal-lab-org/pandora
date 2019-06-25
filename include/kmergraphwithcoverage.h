#ifndef __KMERGRAPHWITHCOVERAGE_H_INCLUDED__   // if kmergraphwithcoverage.h hasn't been included yet...
#define __KMERGRAPHWITHCOVERAGE_H_INCLUDED__

class LocalPRG;

#include <cstdint>
#include <vector>
#include <iostream>
#include "prg/path.h"
#include "kmernode.h"
#include "kmergraph.h"
#include "pangenome/ns.cpp"

/**
 * Represents an annotated KmerGraph, where the annotation is only the coverage on the nodes
 * This is more lightweight than copying a whole k-mer graph and changing the coverage of each node
 * TODO: generalize it to general annotations?
 * NOTE: representing as composition instead of inheritance because one KmerGraph can be referenced by several KmerGraphWithCoverage,
 * thus KmerGraph's life cycle does not depend on KmerGraphWithCoverage's life cycle
 */
class KmerGraphWithCoverage {
private:
    //each coverage is represented as a std::vector, each position of the vector representing the coverage of a sample
    //each coverage of a sample is represented by a std::pair, for the forward and reverse reads coverage of that sample
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> nodeIndex2SampleCoverage;
    uint32_t exp_depth_covg;
    float p;
    float nb_p;
    float nb_r;
    int thresh;
    uint32_t total_number_samples;
    uint32_t num_reads;

public:
    KmerGraph * kmer_prg; //the underlying KmerGraph - TODO: it is dangerous to leave this public, make it private? - this should be const

    //constructor, destructors, etc
    KmerGraphWithCoverage(KmerGraph * kmer_prg, uint32_t total_number_samples=1) :
            nodeIndex2SampleCoverage(kmer_prg->nodes.size()),
            exp_depth_covg{0}, p{1}, nb_p{0.015}, nb_r{2}, thresh{-25}, total_number_samples{total_number_samples}, num_reads{0}, kmer_prg{kmer_prg} {
        assert(kmer_prg != nullptr);
        zeroCoverages();
    }
    KmerGraphWithCoverage(const KmerGraphWithCoverage &other) = default; //copy default constructor
    KmerGraphWithCoverage(KmerGraphWithCoverage &&other) = default; //move default constructor
    KmerGraphWithCoverage& operator=(const KmerGraphWithCoverage& other) = default; //copy assignment operator
    KmerGraphWithCoverage& operator=(KmerGraphWithCoverage&& other) = default; //move assignment operator
    virtual ~KmerGraphWithCoverage() = default;

    //getter
    uint32_t get_covg(uint32_t node_id, bool strand, uint32_t sample_id) const;
    uint32_t get_num_reads() const { return num_reads; }
    uint32_t get_total_number_samples() const {return total_number_samples; }

    //setters
    void increment_covg(uint32_t node_id, bool strand, uint32_t sample_id);
    void set_covg(uint32_t node_id, uint32_t value, bool strand, uint32_t sample_id);
    void set_exp_depth_covg(const uint32_t);
    void set_p(const float);
    void set_nb(const float &, const float &);
    void set_thresh (int thresh) { this->thresh = thresh; }
    void set_num_reads(uint32_t num_reads) { this->num_reads = num_reads; }

    void zeroCoverages() {
        for (auto &sampleCoverage: nodeIndex2SampleCoverage)
            sampleCoverage = std::vector<std::pair<uint32_t, uint32_t>>(total_number_samples);
    }

    float nb_prob(uint32_t, const uint32_t &sample_id);

    float lin_prob(uint32_t, const uint32_t &sample_id);

    float prob(uint32_t, const uint32_t &sample_id);

    float prob(const uint32_t &, const uint32_t &, const uint32_t &sample_id);

    float get_prob(const std::string& prob_model, const uint32_t &node_id, const uint32_t &sample_id);

    bool coverage_is_zeroes(const uint32_t&);

    float find_max_path(std::vector<KmerNodePtr> &maxpath, const std::string& prob_model,
                        const uint32_t &max_num_kmers_to_average, const uint32_t &sample_id);

    std::vector<std::vector<KmerNodePtr>> find_max_paths(uint32_t, const uint32_t &sample_id);

    void save_covg_dist(const std::string &);

    std::vector<std::vector<KmerNodePtr>> get_random_paths(uint32_t);

    float prob_path(const std::vector<KmerNodePtr> &kpath, const uint32_t &sample_id, const std::string& prob_model);

    float prob_paths(const std::vector<std::vector<KmerNodePtr>> &);


    void save(const std::string &, const std::shared_ptr<LocalPRG> = nullptr);
    void load(const std::string &);


    //test friends
    friend class KmerGraphWithCoverageTest_set_p_Test;
    friend class KmerGraphWithCoverageTest_prob_Test;
    friend class KmerGraphWithCoverageTest_findMaxPathSimple_Test;
    friend class KmerGraphWithCoverageTest_findMaxPath2Level_Test;
    friend class KmerGraphWithCoverageTest_find_max_paths_2Level_Test;
    friend class KmerGraphWithCoverageTest_path_prob_Test;
    friend class KmerGraphWithCoverageTest_path_probs_Test;
};

#endif
