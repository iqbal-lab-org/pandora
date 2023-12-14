#ifndef __LOCALPRG_H_INCLUDED__ // if localPRG.h hasn't been included yet...
#define __LOCALPRG_H_INCLUDED__

#include <cstring>
#include <cstdint>
#include <vector>
#include <iostream>
#include <memory>
#include "interval.h"
#include "index.h"
#include "localgraph.h"
#include "prg/path.h"
#include "pangenome/pannode.h"
#include "kmergraph.h"
#include "kmergraphwithcoverage.h"
#include "vcf.h"
#include "fastaq.h"
#include <boost/filesystem.hpp>
#include "globals.h"

class IndexingLimitReached : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

/**
 * Represents a PRG of the many given as input to pandora
 */
class LocalPRG {
private:
    uint32_t next_id; // TODO: this should be definitely removed - it is a variable that
                      // works only in a method, not an object variable
    std::string buff; // TODO: this should be definitely removed - it is a variable that
                      // works only in a method, not an object variable
    std::vector<LocalNodePtr> nodes_along_path_core(const prg::Path&) const;

    static void check_if_vector_of_subintervals_is_consistent_with_envelopping_interval(
        const std::vector<Interval>& subintervals,
        const Interval& envelopping_interval);

    void check_if_we_already_indexed_too_many_kmers(const uint32_t num_kmers_added,
        const uint32_t indexing_upper_bound) const;

    void add_node_to_current_leaves(const KmerNodePtr &kn,
        std::deque<KmerNodePtr> &current_leaves, uint32_t indexing_upper_bound) const;

public:
    uint32_t next_site; // denotes the id of the next variant site to be processed -
                        // TODO: maybe this should not be an object variable
    uint32_t id; // id of this LocalPRG in the full graph (first gene is 0, second is 1,
                 // and so on...)
    std::string name; // name (fasta comment)
    std::string seq; // seq of LocalPRG (the PRG as string itself)
    LocalGraph prg; // the graph that represents this LocalPRG
    KmerGraph kmer_prg; // the kmer sketch graph
    // VCF vcf;
    std::vector<uint32_t> num_hits;

    static bool do_path_memoization_in_nodes_along_path_method;

    LocalPRG(uint32_t id, const std::string& name, const std::string& seq);

    bool operator==(const LocalPRG &other) const {
        if (this->id != other.id) return false;
        if (this->seq != other.seq) return false;
        if (this->prg != other.prg) return false;
        if (this->kmer_prg != other.kmer_prg) return false;
        return true;
    }

    // functions used to create LocalGraph from PRG string, and to sketch graph
    bool isalpha_string(const std::string&) const;

    std::string string_along_path(const prg::Path&) const;

    static std::string string_along_path(const std::vector<LocalNodePtr>&);

    std::vector<LocalNodePtr> nodes_along_path(prg::Path&) const;

    std::vector<Interval> split_by_site(const Interval&) const;

    std::vector<uint32_t> build_graph(
        const Interval&, const std::vector<uint32_t>&, uint32_t current_level = 0);

    std::vector<PathPtr> shift(prg::Path) const;

    void minimizer_sketch(Index* index, const uint32_t w,
        const uint32_t k,
        const uint32_t indexing_upper_bound=INDEXING_UPPER_BOUND_DEFAULT,
        double percentageDone = -1.0);

    // functions used once hits have been collected against the PRG
    std::vector<KmerNodePtr> kmernode_path_from_localnode_path(
        const std::vector<LocalNodePtr>&) const;

    std::vector<LocalNodePtr> localnode_path_from_kmernode_path(
        const std::vector<KmerNodePtr>&, const uint32_t w = 0) const;

    void write_covgs_to_file(
        const boost::filesystem::path&, const std::vector<uint32_t>&) const;

    void write_path_to_fasta(const boost::filesystem::path&,
        const std::vector<LocalNodePtr>&, const float&) const;

    void append_path_to_fasta(const boost::filesystem::path&,
        const std::vector<LocalNodePtr>&, const float&) const;

    void write_aligned_path_to_fasta(const boost::filesystem::path&,
        const std::vector<LocalNodePtr>&, const float&) const;

    void add_new_records_and_genotype_to_vcf_using_max_likelihood_path_of_the_sample(
        VCF& vcf, const std::vector<LocalNodePtr>& rpath,
        const std::vector<LocalNodePtr>& sample_path,
        const std::string& sample_name = "sample") const;

    std::vector<LocalNodePtr> find_alt_path(const std::vector<LocalNodePtr>&,
        const uint32_t, const std::string&, const std::string&) const;

    template<typename RNG>
    std::string random_path(RNG &&rng) {
        std::vector<LocalNodePtr> npath;
        npath.push_back(prg.nodes.at(0));
        while (not npath.back()->outNodes.empty()) {
            uint32_t random_number = rng();
            size_t random_neighbour = random_number % npath.back()->outNodes.size();
            npath.push_back(npath.back()->outNodes[random_neighbour]);
        }
        return string_along_path(npath);
    }

    // TODO: I really feel like these methods are not responsability of a LocalPRG
    // TODO: many of them should be in VCF class, or in the KmerGraphWithCoverage or
    // Fastaq
    void build_vcf_from_reference_path(
        VCF& vcf, const std::vector<LocalNodePtr>& ref) const;

    virtual std::pair<std::vector<uint32_t>, std::vector<uint32_t>>
    get_forward_and_reverse_kmer_coverages_in_range(
        const KmerGraphWithCoverage& kmer_graph_with_coverage,
        const std::vector<KmerNodePtr>& kmer_path,
        const std::vector<LocalNodePtr>& local_path, const uint32_t& range_pos_start,
        const uint32_t& range_pos_end, const uint32_t& sample_id) const;

protected: // helper methods of get_forward_and_reverse_kmer_coverages_in_range():
    virtual uint32_t get_number_of_bases_in_local_path_before_a_given_position(
        const std::vector<LocalNodePtr>& local_path, uint32_t position) const;

    virtual uint32_t get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node(
        const KmerNodePtr& previous_kmer_node,
        const KmerNodePtr& current_kmer_node) const;

public:
    void add_sample_covgs_to_vcf(VCF& vcf, const KmerGraphWithCoverage& kg,
        const std::vector<LocalNodePtr>& ref_path, const std::string& sample_name,
        const uint32_t& sample_id) const;

    void add_consensus_path_to_fastaq(Fastaq&, pangenome::NodePtr, std::vector<KmerNodePtr>&,
        std::vector<LocalNodePtr>&, const uint32_t, const bool, const uint32_t,
        const uint32_t& max_num_kmers_to_average, const uint32_t& sample_id,
        float min_absolute_gene_coverage, float min_relative_gene_coverage,
        float max_relative_gene_coverage, float min_gene_coverage_proportion, bool no_gene_coverage_filtering) const;
    std::vector<LocalNodePtr> get_valid_vcf_reference(const std::string&) const;

    void add_variants_to_vcf(VCF&, pangenome::NodePtr, const std::string&,
        const std::vector<KmerNodePtr>&, const std::vector<LocalNodePtr>&,
        const uint32_t& sample_id = 0, const std::string& sample_name = "sample");

    // friends definitions
    friend std::ostream& operator<<(std::ostream& out, const LocalPRG& data);
    friend class prg::Path; // for memoization
};

bool operator<(const std::pair<std::vector<LocalNodePtr>, float>& p1,
    const std::pair<std::vector<LocalNodePtr>, float>& p2);

bool operator!=(
    const std::vector<KmerNodePtr>& lhs, const std::vector<KmerNodePtr>& rhs);

std::vector<uint32_t> get_covgs_along_localnode_path(const pangenome::NodePtr,
    const std::vector<LocalNodePtr>&, const std::vector<KmerNodePtr>&,
    const uint32_t& sample_id);

#endif
