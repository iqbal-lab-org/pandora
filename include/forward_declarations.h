#ifndef PANDORA_FORWARD_DECLARATIONS_H
#define PANDORA_FORWARD_DECLARATIONS_H

#include <memory>
#include <set>
#include <tuple>

struct MinimizerHit;
class MinimizerHits;
typedef std::shared_ptr<MinimizerHit> MinimizerHitPtr;
struct pComp;
typedef std::set<MinimizerHitPtr, pComp> Hits;
struct clusterComp;
typedef std::set<Hits, clusterComp> MinimizerHitClusters;
class LocalPRG;
class KmerNode;
using SampleIdText = std::string;
using SampleFpath = std::string;
using SampleData = std::pair<SampleIdText, SampleFpath>;

#endif // PANDORA_FORWARD_DECLARATIONS_H
