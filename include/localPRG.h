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
#include "vcf.h"
#include "fastaq.h"

typedef std::shared_ptr<pangenome::Node> PanNodePtr;

class LocalPRG {
    uint32_t next_id;
    std::string buff;
public:
    uint32_t next_site;
    uint32_t id;
    std::string name;
    std::string seq;
    LocalGraph prg;
    KmerGraph kmer_prg;
    //VCF vcf;
    std::vector<uint32_t> num_hits;

    LocalPRG(uint32_t, std::string, std::string);

    // functions used to create LocalGraph from PRG string, and to sketch graph
    bool isalpha_string(const std::string &) const;

    std::string string_along_path(const Path &) const;

    std::string string_along_path(const std::vector<LocalNodePtr> &) const;

    std::vector<LocalNodePtr> nodes_along_path(const Path &) const;

    std::vector<Interval> split_by_site(const Interval &) const;

    std::vector<uint32_t> build_graph(const Interval &,
                                      const std::vector<uint32_t> &,
                                      uint32_t current_level = 0);

    std::vector<Path> shift(Path) const;

    void minimizer_sketch(Index *idx, const uint32_t w, const uint32_t k);

    // functions used once hits have been collected against the PRG
    std::vector<KmerNodePtr> kmernode_path_from_localnode_path(const std::vector<LocalNodePtr> &) const;

    std::vector<LocalNodePtr>
    localnode_path_from_kmernode_path(const std::vector<KmerNodePtr> &, const uint32_t w = 0) const;

    void write_covgs_to_file(const string &, const std::vector<uint32_t> &) const;

    void write_path_to_fasta(const std::string &, const std::vector<LocalNodePtr> &, const float &) const;

    void append_path_to_fasta(const std::string &, const std::vector<LocalNodePtr> &, const float &) const;

    void write_aligned_path_to_fasta(const std::string &,
                                     const std::vector<LocalNodePtr> &,
                                     const float &) const;

    void build_vcf(VCF &, const std::vector<LocalNodePtr> &) const;

    void add_sample_gt_to_vcf(VCF &,
                              const std::vector<LocalNodePtr> &,
                              const std::vector<LocalNodePtr> &,
                              const std::string &sample_name = "sample") const;

    std::vector<LocalNodePtr> find_alt_path(const std::vector<LocalNodePtr> &,
                                            const uint32_t,
                                            const std::string &,
                                            const std::string &) const;

    void append_kmer_covgs_in_range(const KmerGraph &,
                               const std::vector<KmerNodePtr> &,
                               const std::vector<LocalNodePtr> &,
                               const uint32_t &,
                               const uint32_t &,
                               std::vector<uint32_t> &,
                               std::vector<uint32_t> &) const;

    void add_sample_covgs_to_vcf(VCF &, const KmerGraph &,
                                 const std::vector<LocalNodePtr> &,
                                 const std::vector<KmerNodePtr> &sample_kmer_path,
                                 const std::string &sample_name = "sample") const;

    void add_consensus_path_to_fastaq (Fastaq&,
                                       PanNodePtr,
                                       std::vector<KmerNodePtr>&,
                                       std::vector<LocalNodePtr>&,
                                       const uint32_t w,
                                       const bool bin = false,
                                       const uint32_t global_covg = 1);

    void add_variants_to_vcf(VCF&,
                             PanNodePtr,
                             const std::string&,
                             const std::vector<KmerNodePtr>&,
                             const std::vector<LocalNodePtr>&,
                             const std::string& sample_name="sample");

    std::vector<KmerNodePtr> find_path_and_variants(PanNodePtr,
                                                    const std::string &,
                                                    const uint32_t w,
                                                    const std::string &vcf_ref,
                                                    const bool output_vcf = false,
                                                    const bool output_comparison_paths = false,
                                                    const bool output_covgs = false,
                                                    const bool bin = false,
                                                    const uint32_t global_covg = 1,
                                                    const bool regenotype = false) const;

    friend std::ostream &operator<<(std::ostream &out, const LocalPRG &data);
};

bool operator<(const std::pair<std::vector<LocalNodePtr>, float> &p1,
               const std::pair<std::vector<LocalNodePtr>, float> &p2);

bool operator!=(const std::vector<KmerNodePtr> &lhs, const std::vector<KmerNodePtr> &rhs);

std::vector<uint32_t> get_covgs_along_localnode_path(const PanNodePtr,
                                                     const std::vector<LocalNodePtr> &,
                                                     const std::vector<KmerNodePtr> &);

#endif
