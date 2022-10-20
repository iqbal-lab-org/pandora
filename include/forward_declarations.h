#ifndef PANDORA_FORWARD_DECLARATIONS_H
#define PANDORA_FORWARD_DECLARATIONS_H

#include <memory>
#include <set>
#include <tuple>

struct MinimizerHit;
typedef std::shared_ptr<MinimizerHit> MinimizerHitPtr;
class MinimizerHits;
typedef std::set<MinimizerHits> MinimizerHitClusters;
class LocalPRG;
class KmerNode;
using SampleIdText = std::string;
using SampleFpath = std::string;
using SampleData = std::pair<SampleIdText, SampleFpath>;

#endif // PANDORA_FORWARD_DECLARATIONS_H
