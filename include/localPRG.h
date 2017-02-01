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

class LocalNode;
class PanNode;

using namespace std;

class LocalPRG {
    uint32_t next_id;
    string buff;
  public:
    uint32_t next_site;

    uint32_t id;
    string name;
    string seq;
    LocalGraph prg;

    vector<Path> kmer_paths; // added during index construction when PRG is sketched

    vector<uint32_t> kmer_path_clashes; // number of places elsewhere in PRG where a kmer occurs TO DO!!
    vector<vector<uint32_t>> kmer_path_hit_counts; //[forward vector, reverse vector, both vector] 
    vector<vector<float>> kmer_path_probs; //[forward vector, reverse vector, both vector]
    map<uint32_t, vector<vector<MaxPath>>> max_path_index;         // maps from pre_site_ids seen before to a vector
								   // of size three representing (forward, reverse, both), 
								   // each vector containing the maximally probable node_path, 
								   // which mini kmers overlap
								   // the path and its confidence
								   // used when reads suggest this PRG is present in a sample
								   // in the process of working out which path through the PRG is present
								   // Also used to store partial paths (hence vector) during this process
    LocalPRG(uint32_t, string, string);
    //~LocalPRG();
    bool isalpha_string(const string&);
    string string_along_path(const Path&);
    vector<LocalNode*> path_corresponding_to_string(const string&);
    vector<LocalNode*> nodes_along_path(const Path&);
    vector<Interval> splitBySite(const Interval&);	
    vector<uint32_t> build_graph(const Interval&, const vector<uint32_t>&, uint32_t current_level=0);
    void minimizer_sketch (Index* idx, const uint32_t w, const uint32_t k);
    void update_covg_with_hit(MinimizerHit*);
    void update_kmers_on_node_path(MaxPath&, const vector<float>&);
    void update_kmers_on_node_paths(vector<MaxPath>&);
    void get_kmer_path_hit_counts(const PanNode*);
    void get_kmer_path_probs(const PanNode*, uint32_t, float);
    void infer_most_likely_prg_paths_for_corresponding_pannode(const PanNode*, uint32_t, float);
    void write_path_vs_found_path(const string&, const string&);
    void write_max_paths_to_fasta(const string&);

  friend ostream& operator<< (ostream& out, const LocalPRG& data);  
};

bool operator < (const pair<vector<LocalNode*>, float> &p1, const pair<vector<LocalNode*>, float> &p2);

#endif
