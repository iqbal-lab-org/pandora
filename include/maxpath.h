#ifndef __MAXPATH_H_INCLUDED__   // if maxpath.h hasn't been included yet...
#define __MAXPATH_H_INCLUDED__

#include <ostream>
#include <vector>
#include "localnode.h"

struct MaxPath
{
    std::vector<LocalNode*> npath; //node path
    std::vector<int> kmers_on_path; // indicates which of the localPRG hits lie along node path
    uint32_t num_equivalent_paths;
    float prob;
    float mean_prob;
    float median_prob;

    std::string direction;

    MaxPath();
    MaxPath(std::vector<LocalNode*>, std::vector<int>, uint32_t);

    void extend(const MaxPath);
    float get_prob(const std::vector<float>&);
    float get_mean_prob(const std::vector<float>&);
    float get_median_prob(const std::vector<float>&);
};

struct VMPgreater
{
    bool operator()( const std::vector<MaxPath>& lx, const std::vector<MaxPath>& rx );
};

#endif
