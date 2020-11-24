#ifndef __ESTIMATEPARAMETERS_H_INCLUDED__ // if estimate_parameters.h hasn't been
                                          // included yet...
#define __ESTIMATEPARAMETERS_H_INCLUDED__

#include <cstring>
#include <cstdint>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

double fit_mean_covg(const std::vector<uint32_t>&, const uint8_t);

double fit_variance_covg(const std::vector<uint32_t>&, double&, const uint8_t);

void fit_negative_binomial(double&, double&, float&, float&);

uint32_t find_mean_covg(std::vector<uint32_t>&);

int find_prob_thresh(std::vector<uint32_t>&);

uint32_t estimate_parameters(std::shared_ptr<pangenome::Graph> pangraph,
    const fs::path& outdir, const uint32_t k, float& e_rate, const uint32_t covg,
    bool& bin, const uint32_t& sample_id);

#endif
