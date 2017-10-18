#ifndef __UTILS_H_INCLUDED__   // if utils.h hasn't been included yet...
#define __UTILS_H_INCLUDED__

#include <vector>
#include <set>
#include <memory>
#include "localPRG.h"
#include "minihits.h"

class Index;
struct MinimizerHit;
class PanGraph;
class PanNode;
class LocalPRG;
class Seq;

template <typename T>
struct pointer_values_equal
{
    const T* to_find;
    bool operator()(const T* other) const
    {
        return *to_find == *other;
    }
};

template <typename T>
struct spointer_values_equal
{
    const std::shared_ptr<T> to_find;
    bool operator()(const std::shared_ptr<T> other) const
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
void define_clusters(std::set<std::set<MinimizerHitPtr, pComp>,clusterComp>&, const std::vector<LocalPRG*>&,
                     const MinimizerHits*, const int, const float&, const uint, const uint);
void filter_clusters(std::set<std::set<MinimizerHitPtr, pComp>,clusterComp>&);
void filter_clusters2(std::set<std::set<MinimizerHitPtr, pComp>,clusterComp>&, const uint&);
void infer_localPRG_order_for_reads(const std::vector<LocalPRG*>& prgs, MinimizerHits*, PanGraph*, const int,
                                    const uint&, const float&, const uint min_cluster_size=10,
                                    const uint short_read_length=std::numeric_limits<uint>::max());
void pangraph_from_read_file(const std::string&, MinimizerHits*, PanGraph*, Index*, const std::vector<LocalPRG*>&,
                             const uint32_t, const uint32_t, const int, const float&, const uint min_cluster_size=10,
                             const uint genome_size=5000000, const bool illumina=false);
void update_localPRGs_with_hits(PanGraph*, const std::vector<LocalPRG*>&);//, const uint32_t, const float&, bool);
void infer_most_likely_prg_path_for_pannode(const std::vector<LocalPRG*>&, PanNode*, uint32_t, float);
#endif
