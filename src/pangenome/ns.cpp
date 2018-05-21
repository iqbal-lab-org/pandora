#include <iostream>
#include <memory>


using namespace std;

namespace pangenome {
    class Node;

    class Read;

    class Sample;

    class Graph;

    typedef std::shared_ptr<pangenome::Node> NodePtr;
    typedef std::shared_ptr<pangenome::Read> ReadPtr;
    typedef std::shared_ptr<pangenome::Sample> SamplePtr;
}
