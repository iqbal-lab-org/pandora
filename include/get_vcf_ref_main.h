#ifndef PANDORA_GET_VCF_REF_MAIN_H
#define PANDORA_GET_VCF_REF_MAIN_H
#include <cassert>
#include <vector>
#include <iostream>
#include "localPRG.h"
#include "utils.h"
#include "fastaq_handler.h"
#include "fastaq.h"
#include "CLI11.hpp"

struct GetVcfRefOptions {
    std::string prgfile;
    std::string seqfile;
    bool compress { false };
    uint8_t verbosity { 0 };
};

void setup_get_vcf_ref_subcommand(CLI::App& app);
int pandora_get_vcf_ref(GetVcfRefOptions const& opt);
#endif // PANDORA_GET_VCF_REF_MAIN_H
