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
#include "kmergraph.h"

class LocalNode;
class PanNode;
class KmerNode;

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

    LocalPRG(uint32_t, std::string, std::string);

    // functions used to create LocalGraph from PRG string, and to sketch graph
    bool isalpha_string(const std::string&);
    std::string string_along_path(const Path&);
    std::vector<LocalNode*> nodes_along_path(const Path&);
    std::vector<Interval> split_by_site(const Interval&);	
    std::vector<uint32_t> build_graph(const Interval&, const std::vector<uint32_t>&, uint32_t current_level=0);
    std::vector<Path> shift(Path);
    void minimizer_sketch (Index* idx, const uint32_t w, const uint32_t k);

    // functions used once hits have been collected against the PRG
    std::vector<KmerNode*> find_kmernodes_on_localnode_path(std::vector<LocalNode*>&);
    void update_covg_with_hit(MinimizerHit*);
    std::vector<LocalNode*> localnode_path_from_kmernode_path(std::vector<KmerNode*>);
    void write_kmer_max_paths_to_fasta(const std::string&, float);

  friend std::ostream& operator<< (std::ostream& out, const LocalPRG& data);  
};

bool operator < (const std::pair<std::vector<LocalNode*>, float> &p1, const std::pair<std::vector<LocalNode*>, float> &p2);

#endif
