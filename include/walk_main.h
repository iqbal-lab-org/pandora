#ifndef PANDORA_WALK_MAIN_H
#define PANDORA_WALK_MAIN_H
#include <cstring>
#include <vector>
#include <iostream>

#include "localPRG.h"
#include "utils.h"
#include "fastaq_handler.h"
#include "CLI11.hpp"

struct WalkOptions {
    fs::path pandora_index_file;
    fs::path seqfile;
    bool top { false };
    bool bottom { false };
    uint8_t verbosity { 0 };
};

void setup_walk_subcommand(CLI::App& app);
int pandora_walk(WalkOptions const& opt);
#endif // PANDORA_WALK_MAIN_H
