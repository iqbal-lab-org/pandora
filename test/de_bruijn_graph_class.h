#include "de_bruijn/graph.h"
#include "de_bruijn/ns.cpp"

class GraphTester : public debruijn::Graph
{
public:
    GraphTester(uint8_t i) : Graph(i) {};
    friend class DeBruijnGraphTest_create_Test;
    friend class DeBruijnGraphTest_add_node_Test;
    friend class DeBruijnGraphTest_add_edge_Test;
    friend class DeBruijnGraphTest_remove_node_Test;
    friend class DeBruijnGraphTest_get_leaves_Test;
    friend class DeBruijnGraphTest_get_unitigs_Test;
    friend class DeBruijnGraphTest_extend_unitig_Test;
    friend class DeBruijnGraphTest_equals_Test;
};
