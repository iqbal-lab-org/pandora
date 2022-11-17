#ifndef PANDORA_FORWARD_DECLARATIONS_H
#define PANDORA_FORWARD_DECLARATIONS_H

#include <memory>
#include <set>
#include <tuple>
#include <boost/filesystem.hpp>


struct MinimizerHit;
typedef std::shared_ptr<MinimizerHit> MinimizerHitPtr;
class MinimizerHits;
typedef std::set<MinimizerHits> MinimizerHitClusters;
class LocalPRG;
class KmerNode;
typedef std::shared_ptr<KmerNode> KmerNodePtr;
using SampleIdText = std::string;
using SampleFpath = std::string;
using SampleData = std::pair<SampleIdText, SampleFpath>;
namespace fs = boost::filesystem;

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

#endif // PANDORA_FORWARD_DECLARATIONS_H
