#ifndef PANDORA_RANDOM_MAIN_H
#define PANDORA_RANDOM_MAIN_H

#include <vector>
#include <iostream>

#include "localPRG.h"
#include "utils.h"
#include "fastaq.h"
#include "CLI11.hpp"

struct RandomOptions {
    std::string prgfile;
    bool compress { false };
    uint32_t num_paths { 1 };
    uint8_t verbosity { 0 };
};

void setup_random_subcommand(CLI::App& app);
int pandora_random_path(RandomOptions const& opt);
#endif // PANDORA_RANDOM_MAIN_H
