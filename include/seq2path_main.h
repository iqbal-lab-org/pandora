#ifndef PANDORA_SEQ2PATH_MAIN_H
#define PANDORA_SEQ2PATH_MAIN_H
#include <cstring>
#include <cassert>
#include <vector>
#include <iostream>

#include <boost/log/trivial.hpp>

#include "utils.h"
#include "localPRG.h"
#include "fastaq_handler.h"
#include "CLI11.hpp"

struct Seq2PathOptions {
    std::string prgfile;
    std::string seqfile;
    uint32_t kmer_size { 15 };
    uint32_t window_size { 14 };
    bool top { false };
    bool bottom { false };
    bool flag { false };
    uint8_t verbosity { 0 };
};

void setup_seq2path_subcommand(CLI::App& app);
int pandora_seq2path(Seq2PathOptions const& opt);

#endif // PANDORA_SEQ2PATH_MAIN_H
