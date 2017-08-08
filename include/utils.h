#ifndef __UTILS_H_INCLUDED__   // if utils.h hasn't been included yet...
#define __UTILS_H_INCLUDED__

#include <vector>
#include "localPRG.h"

class Index;
class MinimizerHits;
class PanGraph;
class PanNode;
class LocalPRG;

template <typename T>
struct pointer_values_equal
{
    const T* to_find;
    bool operator()(const T* other) const
    {
        return *to_find == *other;
    }
};

// utility functions
std::string now();
std::vector<std::string> split(const std::string&, const std::string&);
char complement(char);
std::string rev_complement(std::string);
float lognchoosek2 (uint32_t, uint32_t, uint32_t);

//probably should be moved to map_main.cpp
void read_prg_file(std::vector<LocalPRG*>&, const std::string&);
void load_PRG_kmergraphs(std::vector<LocalPRG*>&, const uint&, const uint&, const std::string&);
//void add_read_hits(uint32_t, const std::string&, const std::string&, MinimizerHits*, Index*, const uint32_t, const uint32_t);
void add_read_hits(Seq*, MinimizerHits*, Index*);
void infer_localPRG_order_for_reads(const std::vector<LocalPRG*>& prgs, MinimizerHits*, PanGraph*, const int, const uint min_cluster_size=5);
//void pangraph_from_read_file_new(const std::string&, MinimizerHits*, PanGraph*, Index*, const std::vector<LocalPRG*>&, const uint32_t, const uint32_t, const int);
void pangraph_from_read_file(const std::string&, MinimizerHits*, PanGraph*, Index*, const std::vector<LocalPRG*>&, const uint32_t, const uint32_t, const int, const uint min_cluster_size=5, const uint genome_size=5000000);
void update_localPRGs_with_hits(PanGraph*, const std::vector<LocalPRG*>&);//, const uint32_t, const float&, bool);
//float p_null(const std::vector<LocalPRG*>&, std::set<MinimizerHit*, pComp>&, uint32_t);
void infer_most_likely_prg_path_for_pannode(const std::vector<LocalPRG*>&, PanNode*, uint32_t, float);
//uint32_t nchoosek(const uint32_t , const uint32_t );
#endif
