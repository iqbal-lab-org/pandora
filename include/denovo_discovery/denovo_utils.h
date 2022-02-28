#ifndef __EXTRACT_READS_H_INCLUDED__
#define __EXTRACT_READS_H_INCLUDED__

#include <map>
#include <string>
#include <boost/filesystem/path.hpp>
#include <forward_declarations.h>
#include <iostream>
#include "utils.h"
#include <boost/algorithm/string.hpp>

namespace fs = boost::filesystem;


std::map<std::string, std::string> get_locus_to_reads(
    const fs::path &sample_outdir,
    const std::string &sample_name
);

void concatenate_all_denovo_files(const std::vector<SampleData> &samples,
    const fs::path &outdir);

#endif
