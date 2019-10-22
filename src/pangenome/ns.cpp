#include <iostream>
#include <memory>

namespace pangenome {
class Node;

class Read;

class Sample;
struct SamplePtrSorterBySampleId;

class Graph;

typedef std::shared_ptr<pangenome::Node> NodePtr;
typedef std::weak_ptr<pangenome::Node> WeakNodePtr;
typedef std::shared_ptr<pangenome::Read> ReadPtr;
typedef std::shared_ptr<pangenome::Sample> SamplePtr;
}
