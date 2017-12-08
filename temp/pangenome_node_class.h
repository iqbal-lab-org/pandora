#include "pangenome/pannode.h"

class PNodeTester : public pangenome::Node
{
public:
    PNodeTester(const uint32_t i, const uint32_t j, const string n) : Node(i,j,n) {};
    friend class PangenomeNodeTest_create_Test;
    friend class PangenomeNodeTest_output_samples_Test;

};
