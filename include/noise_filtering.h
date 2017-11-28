
#include "pangenome/pangraph.h"
#include "de_bruijn/graph.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

uint16_t node_plus_orientation_to_num (const uint16_t, const bool);
void num_to_node_plus_orientation (uint16_t&, bool&, const uint16_t);
void hashed_node_ids_to_ids_and_orientations(const deque<uint16_t>&, std::vector<uint16_t>&, std::vector<bool>&);
debruijn::Graph construct_debruijn_graph_from_pangraph(uint8_t, const pangenome::Graph &);
void remove_leaves(pangenome::Graph &, debruijn::Graph &);
void filter_unitigs(pangenome::Graph &, debruijn::Graph &, const uint16_t&);
void detangle_pangraph_with_debruijn_graph(pangenome::Graph &, debruijn::Graph &);
void clean_pangraph_with_debruijn_graph(pangenome::Graph &, debruijn::Graph &, const uint16_t&);
