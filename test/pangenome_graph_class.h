#include "pangenome/pangraph.h"

class PGraphTester : public pangenome::Graph
{
public:
    PGraphTester() : Graph() {};
    friend class PangenomeGraphTest_add_node_Test;
    friend class PangenomeGraphTest_add_node_sample_Test;
    friend class PangenomeGraphTest_clear_Test;
    friend class PangenomeReadTest_remove_node_Test;
    friend class PangenomeReadTest_find_position_Test;

};
