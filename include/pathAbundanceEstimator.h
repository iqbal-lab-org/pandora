#ifndef __Path_ABUNDANCE_ESTIMATOR_H
#define __Path_ABUNDANCE_ESTIMATOR_H

#include "kmergraph.h"

class PathAbundanceEstimator {
private:
  double epsilon;
  uint16_t maxItrCnt;
  std::vector<double> pathCnts;
  std::vector<std::vector<std::pair<uint16_t, double>>> readProbs;	
  bool updatePathCnts();

public:
  // default value for epsilon based on EM termination criteria adopted from Bray et al. 2016
  PathAbundanceEstimator(std::vector<std::vector<std::pair<uint16_t, uint16_t>>> hitCntPerRead4Paths,
                         std::vector<std::deque<KmerNodePtr>> paths,
                         double epsilonIn = 1e-8,
                         uint16_t maxItrCntIn = 100);
  std::vector<double>& runEM();
  std::vector<double>& getPathCnts();
  std::vector<uint16_t> getMostLikelyAssignments();

friend class PathAbundanceEstimatorTest_constructor_Test;
};

#endif
