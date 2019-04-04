#ifndef __ESTIMATEPARAMETERS_H_INCLUDED__   // if estimate_parameters.h hasn't been included yet...
#define __ESTIMATEPARAMETERS_H_INCLUDED__

#include <cstring>
#include <cstdint>
#include <boost/filesystem.hpp>


namespace fs = boost::filesystem;


uint32_t find_mean_covg(std::vector<uint32_t> &);

int find_prob_thresh(std::vector<uint32_t> &);

uint32_t estimate_parameters(std::shared_ptr<pangenome::Graph>, const std::string &, const uint32_t, float &,
                         const uint32_t,
                         bool &bin, const uint32_t &sample_id);

#endif
