#ifndef __NOISEFILTERING_H_INCLUDED__ // if noise_filtering.h hasn't been included
                                      // yet...
#define __NOISEFILTERING_H_INCLUDED__

#include <string>
#include <cstdint>
#include "pangenome/pangraph.h"
#include "de_bruijn/graph.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

uint_least32_t node_plus_orientation_to_num(const uint_least32_t, const bool);

void num_to_node_plus_orientation(uint_least32_t&, bool&, const uint_least32_t);

uint_least32_t rc_num(const uint_least32_t&);

void hashed_node_ids_to_ids_and_orientations(const std::deque<uint_least32_t>&,
    std::vector<uint_least32_t>&, std::vector<bool>&);

bool overlap_forwards(
    const std::deque<uint_least32_t>&, const std::deque<uint_least32_t>&);

bool overlap_backwards(
    const std::deque<uint_least32_t>&, const std::deque<uint_least32_t>&);

std::deque<uint_least32_t> rc_hashed_node_ids(const std::deque<uint_least32_t>&);

void dbg_node_ids_to_ids_and_orientations(const debruijn::Graph&,
    const std::deque<uint32_t>&, std::vector<uint_least32_t>&, std::vector<bool>&);

void construct_debruijn_graph(
    std::shared_ptr<pangenome::Graph> pangraph, debruijn::Graph& dbg);

void remove_leaves(std::shared_ptr<pangenome::Graph>, debruijn::Graph&,
    uint_least32_t covg_thresh = 1);

void filter_unitigs(
    std::shared_ptr<pangenome::Graph>, debruijn::Graph&, const uint_least32_t&);

void detangle_pangraph_with_debruijn_graph(
    std::shared_ptr<pangenome::Graph>, debruijn::Graph&);

void clean_pangraph_with_debruijn_graph(std::shared_ptr<pangenome::Graph>,
    const uint_least32_t, const uint_least32_t, const bool illumina = false);

void write_pangraph_gfa(const std::string&, std::shared_ptr<pangenome::Graph>);

#endif
