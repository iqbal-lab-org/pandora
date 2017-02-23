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
#include "pannode.h"
#include "maxpath.h"
#include "kmergraph.h"

class LocalNode;
class PanNode;

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

    std::vector<Path> kmer_paths; // added during index construction when PRG is sketched

    std::vector<uint32_t> kmer_path_clashes; // number of places elsewhere in PRG where a kmer occurs TO DO!!
    std::vector<std::vector<uint32_t>> kmer_path_hit_counts; //[forward vector, reverse vector, both vector] 
    std::vector<std::vector<float>> kmer_path_probs;         //[forward vector, reverse vector, both vector]
    std::map<uint32_t, std::vector<std::vector<MaxPath>>> max_path_index;          // maps from pre_site_ids seen before to a vector
										   // of size three representing (forward, reverse, both), 
										   // each vector containing the maximally probable node_path, 
										   // which mini kmers overlap
										   // the path and its confidence
										   // used when reads suggest this PRG is present in a sample
										   // in the process of working out which path through the PRG is present
										   // Also used to store partial paths (hence vector) during this process
								   
    LocalPRG(uint32_t, std::string, std::string);

    // functions used to create LocalGraph from PRG string, and to sketch graph
    bool isalpha_string(const std::string&);
    std::string string_along_path(const Path&);
    std::vector<LocalNode*> nodes_along_path(const Path&);
    std::vector<Interval> split_by_site(const Interval&);	
    std::vector<uint32_t> build_graph(const Interval&, const std::vector<uint32_t>&, uint32_t current_level=0);
    void minimizer_sketch (Index* idx, const uint32_t w, const uint32_t k);

    // functions used once hits have been collected against the PRG
    void update_covg_with_hit(MinimizerHit*);
    std::vector<uint32_t> find_overlapping_kmer_paths(MaxPath&);
    void filter_branching_kmer_paths(MaxPath&, const std::vector<float>&, const std::vector<uint32_t>&);
    void update_kmers_on_node_path(MaxPath&, const std::vector<float>&);
    void update_kmers_on_node_paths(std::vector<MaxPath>&);
    void get_kmer_path_hit_counts(const PanNode*);
    void get_kmer_path_probs(const PanNode*, uint32_t, float);
    void infer_most_likely_prg_paths_for_corresponding_pannode(const PanNode*, uint32_t, float);
    std::vector<float> get_covered_maxpath_log_probs(uint, uint);
    void write_max_paths_to_fasta(const std::string&);

  friend std::ostream& operator<< (std::ostream& out, const LocalPRG& data);  
};

bool operator < (const std::pair<std::vector<LocalNode*>, float> &p1, const std::pair<std::vector<LocalNode*>, float> &p2);

#endif
