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

inline size_t get_flank_field_pos(const std::vector<std::string> &words, const std::string &flank_prefix);

std::map<std::string, std::string> get_locus_to_reads(
    const fs::path &sample_outdir,
    const std::string &sample_name
);

void concatenate_all_denovo_files(const std::vector<SampleData> &samples,
    const fs::path &outdir);

void write_denovo_header_file(const std::string &sample_name,
    const fs::path &denovo_outdir, uint32_t number_of_loci_with_denovo_variants);

#endif
