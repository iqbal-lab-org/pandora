#include "pangenome/panread.h"

class PReadTester : public pangenome::Read
{
public:
    PReadTester(const uint32_t i) : Read(i) {};
    friend class PangenomeReadTest_create_Test;
    friend class PangenomeReadTest_add_hits_Test;
    friend class PangenomeReadTest_remove_node_Test;


};
