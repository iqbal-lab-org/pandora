#ifndef PANDORA_RACON_H
#define PANDORA_RACON_H

#include <string>
#include "utils.h"

class Racon {
private:
    const bool illumina;
    const uint32_t kmer_size;
    const std::string &locus;
    std::string consensus_seq;
    const fs::path &denovo_outdir;
    const std::string &locus_reads_filepath;
    const uint32_t max_number_of_rounds_to_run;
    const bool keep_extra_debugging_files;
    uint32_t number_of_rounds_executed;
    bool previous_run_improved_consensus_seq;
    std::vector<std::string> consensus_seq_already_seen;

    // Runs minimap2 with the given parameters
    // @return : PAF string
    static std::string run_minimap2_core(
        const std::string &query_filepath,
        const std::string &ref_filepath,
        const bool illumina, const uint32_t kmer_size);


    // Runs racon with the given parameters
    // @return : a polished sequence
    static std::string run_racon_core(
        const std::string &sequence_to_be_corrected,
        const std::string &reads_filepath,
        const std::string &paf_filepath,
        const std::string &target_seq_filepath);


    // Runs a round of minimap2+racon
    // @return : a bool meaning if racon improved the consensus seq or not
    // @side-effects:
    //      1. modifies consensus_seq;
    //      2. increases the number_of_rounds_executed;
    //      3. adds the new polished sequence to consensus_seq_already_seen;
    bool run_another_round();


    // Runs run_another_round() at most max_number_of_rounds_to_run or until it does
    // not improve the consensus seq anymore
    void run();

public:
    Racon(bool illumina, uint32_t kmer_size, const std::string &locus,
          const std::string &consensus_seq, const fs::path &denovo_outdir,
          const std::string &locus_reads_filepath,
          uint32_t max_number_of_rounds_to_run = 10,
          bool keep_extra_debugging_files = false)
        : illumina(illumina), kmer_size(kmer_size), locus(locus), consensus_seq(consensus_seq),
        denovo_outdir(denovo_outdir), locus_reads_filepath(locus_reads_filepath),
        max_number_of_rounds_to_run(max_number_of_rounds_to_run),
        keep_extra_debugging_files(keep_extra_debugging_files),
        number_of_rounds_executed(0), previous_run_improved_consensus_seq(true),
        consensus_seq_already_seen() {}

    const std::string & get_polished_sequence() {
        run();
        return consensus_seq;
    }

private:
    // some predefined constants
    const static uint32_t racon_window_padding = 100;
};

#endif // PANDORA_RACON_H
