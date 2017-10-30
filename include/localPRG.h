#ifndef __LOCALPRG_H_INCLUDED__   // if localPRG.h hasn't been included yet...
#define __LOCALPRG_H_INCLUDED__

#include <string>
#include <vector>
#include <tuple>
#include <set>
#include <ostream>
#include "minimizer.h"
#include "interval.h"
#include "localgraph.h"
#include "path.h"
#include "index.h"
#include "minihits.h"
#include "pangenome/pannode.h"
#include "kmergraph.h"
#include "vcf.h"

class LocalNode;
class KmerNode;

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
    bool isalpha_string(const std::string&) const;
    std::string string_along_path(const Path&) const;
    std::vector<LocalNodePtr> nodes_along_path(const Path&) const;
    std::vector<Interval> split_by_site(const Interval&) const;	
    std::vector<uint32_t> build_graph(const Interval&, const std::vector<uint32_t>&, uint32_t current_level=0);
    std::vector<Path> shift(Path) const;
    void minimizer_sketch (Index* idx, const uint32_t w, const uint32_t k);

    // functions used once hits have been collected against the PRG
    //void update_covg_with_hit(MinimizerHit*);
    std::vector<LocalNodePtr> localnode_path_from_kmernode_path(const std::vector<KmerNodePtr>&, const uint w=0) const;
    void write_path_to_fasta(const std::string&, const std::vector<LocalNodePtr>&, const float&) const;
    void append_path_to_fasta(const std::string&, const std::vector<LocalNodePtr>&, const float&) const;
    void write_aligned_path_to_fasta(const std::string&, const std::vector<LocalNodePtr>&, const float&) const;
    //void build_vcf();
    void build_vcf(VCF&, const std::vector<LocalNodePtr>&) const;
    //void add_sample_to_vcf(const std::vector<LocalNodePtr>&);
    void add_sample_to_vcf(VCF&, const std::vector<LocalNodePtr>&, const std::vector<LocalNodePtr>&, const std::string& sample_name="sample") const;
    std::vector<KmerNodePtr> find_path_and_variants(PanNodePtr, const std::string&, uint w=0, bool max_path=true, bool min_path=false, bool output_vcf = false, bool output_comparison_paths = false) const;

  friend std::ostream& operator<< (std::ostream& out, const LocalPRG& data);  
};

bool operator < (const std::pair<std::vector<LocalNodePtr>, float> &p1, const std::pair<std::vector<LocalNodePtr>, float> &p2);
bool operator != (const std::vector<KmerNodePtr>& lhs, const std::vector<KmerNodePtr>& rhs);
#endif
