#ifndef PANDORA_MERGE_INDEX_MAIN_H
#define PANDORA_MERGE_INDEX_MAIN_H

#include <vector>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>

#include "utils.h"
#include "localPRG.h"
#include "CLI11.hpp"

namespace fs = boost::filesystem;

struct MergeIndexOptions {
    std::vector<std::string> indicies;
    fs::path outfile { "merged_index.idx" };
    int verbosity { 0 };
};

void setup_merge_index_subcommand(CLI::App& app);
int pandora_merge_index(MergeIndexOptions const& opt);
#endif // PANDORA_MERGE_INDEX_MAIN_H
