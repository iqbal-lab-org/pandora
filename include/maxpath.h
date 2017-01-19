#ifndef __MAXPATH_H_INCLUDED__   // if maxpath.h hasn't been included yet...
#define __MAXPATH_H_INCLUDED__

#include <ostream>
#include <vector>
#include "localnode.h"

using namespace std;

struct MaxPath
{
    vector<LocalNode*> npath; //node path
    vector<bool> kmers_on_path; // indicates which of the localPRG hits lie along node path
    uint32_t num_equivalent_paths;
    float prob;
    float mean_prob;

    MaxPath();
    MaxPath(vector<LocalNode*>, vector<bool>, uint32_t);

    void extend(const MaxPath);
    float get_prob(const vector<float>&);
    float get_mean_prob(const vector<float>&);
};

#endif
