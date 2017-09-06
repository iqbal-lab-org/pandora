#ifndef __ESTIMATEPARAMETERS_H_INCLUDED__   // if estimate_parameters.h hasn't been included yet...
#define __ESTIMATEPARAMETERS_H_INCLUDED__

class Pangraph;
class LocalPRG;

uint find_mean_covg(std::vector<uint>&);

int find_prob_thresh(std::vector<uint>&);

void estimate_parameters(PanGraph*, std::string&, uint32_t, float&);

#endif
