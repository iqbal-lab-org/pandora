#include <cmath>
#include "pathAbundanceEstimator.h"

PathAbundanceEstimator::PathAbundanceEstimator(std::vector<std::vector<std::pair<uint16_t, uint16_t>>> hitCntPerRead4Paths,
                                               std::vector<std::deque<KmerNodePtr>> paths,
                                               double epsilonIn,
                                               uint16_t maxItrCntIn):
  epsilon(epsilonIn), maxItrCnt(maxItrCntIn) {
  pathCnts.resize(paths.size());
  std::fill_n(pathCnts.begin(), pathCnts.size(), 1);
  readProbs.resize(hitCntPerRead4Paths.size());

  // calculate probability of each read coming from any of the paths in its compatible list of paths
  uint16_t i = 0;
  for (auto readIt = hitCntPerRead4Paths.begin(); readIt != hitCntPerRead4Paths.end(); readIt++) {
    std::vector<std::pair<uint16_t, double>> readProb;
    // calculate the probability of each read coming from a path : (# of hits) / (path total # of minimizers)
    for (auto hitIt = readIt->begin(); hitIt != readIt->end(); hitIt++) {
      // hitIt->first : path ID
      // hitIt->second : # of hits
      readProb.emplace_back(hitIt->first, static_cast<double>(hitIt->second)/static_cast<double>(paths[hitIt->first].size()));
    }
    readProbs[i++] = readProb;
  }

};

std::vector<double>& PathAbundanceEstimator::runEM() {
  uint16_t iterationCntr = 0;
  bool converged = false;

  uint32_t i = 0;
  std::vector<double> tmpPathCnts(pathCnts.size());
  while (iterationCntr++ < maxItrCnt && !converged) {
    converged = true;
    // E Step
    for (auto readIt = readProbs.begin(); readIt != readProbs.end(); readIt++) {
      // first go over all compatible paths for a read and calculate the denomerator to normalize the probabilities for this read condition on each paths
      double denom = 0;
      // keep track of the new probabilities in newReadProbs just to skip multiplication again
      // NOTE the value calculated in this step is not actually probability but a value proportional to the probability which will be normalized in next step
      std::vector<double> newReadProbs(readIt->size());
      i = 0;
      for (auto probIt = readIt->begin(); probIt != readIt->end(); probIt++) {
        newReadProbs[i] = pathCnts[probIt->first] * probIt->second;
        denom += newReadProbs[i++];
      }
      // normalize probabilities for each read

      i = 0;
      for (auto probIt = readIt->begin(); probIt != readIt->end(); probIt++) {
        tmpPathCnts[probIt->first] += newReadProbs[i++] / denom;
      }
    }

    // M Step
    // update path counts
    for (i = 0; i < tmpPathCnts.size(); i++) {
      if (std::abs(tmpPathCnts[i] - pathCnts[i])/pathCnts[i] > epsilon) {
        converged = false;
      }
      pathCnts[i] = tmpPathCnts[i];
      tmpPathCnts[i] = 0;
    }
  }
  if (iterationCntr == maxItrCnt)
    std::cerr << "INFO: Didn't converge.\n";
  // TODO question: do we need to replace p(p|gama)p(r|p) instead of p(r|p) in the readProbs vector at the end?

  return pathCnts;
};

std::vector<double>& PathAbundanceEstimator::getPathCnts() {
  return pathCnts;
}
