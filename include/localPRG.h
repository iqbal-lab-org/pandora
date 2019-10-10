#ifndef __LOCALPRG_H_INCLUDED__   // if localPRG.h hasn't been included yet...
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


using PanNodePtr = std::shared_ptr<pangenome::Node>;
namespace fs = boost::filesystem;

/**
 * Represents a PRG of the many given as input to pandora
 */
class LocalPRG {
private:
    uint32_t next_id; //TODO: this should be definitely removed - it is a variable that works only in a method, not an object variable
    std::string buff; //TODO: this should be definitely removed - it is a variable that works only in a method, not an object variable
    std::vector<LocalNodePtr> nodes_along_path_core(const prg::Path &) const;

public:
    uint32_t next_site; //denotes the id of the next variant site to be processed - TODO: maybe this should not be an object variable
    uint32_t id; //id of this LocalPRG in the full graph (first gene is 0, second is 1, and so on...)
    std::string name; //name (fasta comment)
    std::string seq; //seq of LocalPRG (the PRG as string itself)
    LocalGraph prg; //the graph that represents this LocalPRG
    KmerGraph kmer_prg; //the kmer sketch graph
    //VCF vcf;
    std::vector<uint32_t> num_hits;

    static bool do_path_memoization_in_nodes_along_path_method;

    LocalPRG(uint32_t id, const std::string &name, const std::string &seq);

    // functions used to create LocalGraph from PRG string, and to sketch graph
    bool isalpha_string(const std::string &) const;

    std::string string_along_path(const prg::Path &) const;

    static std::string string_along_path(const std::vector<LocalNodePtr> &);

    std::vector<LocalNodePtr> nodes_along_path(prg::Path &) const;

    std::vector<Interval> split_by_site(const Interval &) const;

    std::vector<uint32_t> build_graph(const Interval &,
                                      const std::vector<uint32_t> &,
                                      uint32_t current_level = 0);

    std::vector<PathPtr> shift(prg::Path) const;

    void minimizer_sketch(const std::shared_ptr<Index> &index, const uint32_t w, const uint32_t k, double percentageDone=-1.0);

    // functions used once hits have been collected against the PRG
    std::vector<KmerNodePtr> kmernode_path_from_localnode_path(const std::vector<LocalNodePtr> &) const;

    std::vector<LocalNodePtr>
    localnode_path_from_kmernode_path(const std::vector<KmerNodePtr> &, const uint32_t w = 0) const;

    void write_covgs_to_file(const boost::filesystem::path &, const std::vector<uint32_t> &) const;

    void write_path_to_fasta(const boost::filesystem::path &, const std::vector<LocalNodePtr> &, const float &) const;

    void append_path_to_fasta(const boost::filesystem::path &, const std::vector<LocalNodePtr> &, const float &) const;

    void write_aligned_path_to_fasta(const boost::filesystem::path &,
                                     const std::vector<LocalNodePtr> &,
                                     const float &) const;

    void add_sample_gt_to_vcf(VCF &,
                              const std::vector<LocalNodePtr> &,
                              const std::vector<LocalNodePtr> &,
                              const std::string &sample_name = "sample") const;

    std::vector<LocalNodePtr> find_alt_path(const std::vector<LocalNodePtr> &,
                                            const uint32_t,
                                            const std::string &,
                                            const std::string &) const;

    std::string random_path();

    //TODO: I really feel like these methods are not responsability of a LocalPRG
    //TODO: many of them should be in VCF class, or in the KmerGraphWithCoverage or Fastaq
    void build_vcf(VCF &, const std::vector<LocalNodePtr> &) const;


    virtual std::pair<std::vector<uint32_t>, std::vector<uint32_t>>
    append_kmer_covgs_in_range(const KmerGraphWithCoverage &kmer_graph_with_coverage,
                               const std::vector<KmerNodePtr> &kmer_path,
                               const std::vector<LocalNodePtr> &local_path, const uint32_t &range_pos_start,
                               const uint32_t &range_pos_end, const uint32_t &sample_id) const;

protected: //helper methods of append_kmer_covgs_in_range():
    virtual uint32_t
    get_number_of_bases_in_local_path_before_a_given_position(const std::vector<LocalNodePtr> &local_path,
                                                              uint32_t position) const;

    virtual uint32_t
    get_number_of_bases_that_are_exclusively_in_the_previous_kmer_node(const KmerNodePtr &previous_kmer_node,
                                                                       const KmerNodePtr &current_kmer_node) const;

public:

    void add_sample_covgs_to_vcf(VCF &, const KmerGraphWithCoverage &, const std::vector<LocalNodePtr> &,
                                     const uint32_t &min_kmer_covg, const std::string &sample_name="sample",
                                     const uint32_t &sample_id=0) const;

    void add_consensus_path_to_fastaq(Fastaq &, PanNodePtr, std::vector<KmerNodePtr> &, std::vector<LocalNodePtr> &,
                                          const uint32_t, const bool, const uint32_t,
                                          const uint32_t &max_num_kmers_to_average, const uint32_t &sample_id) const;
    std::vector<LocalNodePtr> get_valid_vcf_reference(const std::string &) const;

    void add_variants_to_vcf(VCF &, PanNodePtr, const std::string &, const std::vector<KmerNodePtr> &,
                                 const std::vector<LocalNodePtr> &, const uint32_t &min_kmer_covg,
                                 const uint32_t &sample_id=0, const std::string &sample_name="sample");


    //friends definitions
    friend std::ostream &operator<<(std::ostream &out, const LocalPRG &data);
    friend class prg::Path; //for memoization
};

bool operator<(const std::pair<std::vector<LocalNodePtr>, float> &p1,
               const std::pair<std::vector<LocalNodePtr>, float> &p2);

bool operator!=(const std::vector<KmerNodePtr> &lhs, const std::vector<KmerNodePtr> &rhs);

std::vector<uint32_t>
get_covgs_along_localnode_path(const PanNodePtr, const std::vector<LocalNodePtr> &, const std::vector<KmerNodePtr> &,
                               const uint32_t &sample_id);

#endif
