#ifndef PANDORA_RACON_H
#define PANDORA_RACON_H

#include <string>
#include "utils.h"

class Racon {
private:
    const bool illumina;
    const std::string &locus;
    std::string consensus_seq;
    const fs::path &denovo_outdir;
    const std::string &locus_reads_filepath;
    uint32_t max_number_of_rounds_to_run;
    uint32_t number_of_rounds_executed;
    bool previous_run_improved_consensus_seq;
    std::vector<std::string> consensus_seq_already_seen;

    // Runs racon with the given parameters and returns a polished sequence
    static std::string run_racon_core(
        const std::string &consensus_seq,
        const std::string &reads_filepath,
        const std::string &paf_filepath,
        const std::string &polished_filepath);


    // Runs a round of minimap2+racon, effectively modifying consensus_seq and
    // increasing the number_of_rounds_executed.
    // Returns a bool, which means if racon improved the consensus seq or not.
    // If this bool is false, then it means racon cannot improve the consensus seq
    // anymore, and we should stop running it
    bool run_another_round();


    // Run racon for at most max_number_of_rounds_to_run (or until it does not improve
    // the consensus seq anymore)
    void run();

public:
    Racon(bool illumina, const std::string &locus,
          const std::string &consensus_seq, const fs::path &denovo_outdir,
          const std::string &locus_reads_filepath,
          uint32_t max_number_of_rounds_to_run = 10)
        : illumina(illumina), locus(locus), consensus_seq(consensus_seq),
        denovo_outdir(denovo_outdir), locus_reads_filepath(locus_reads_filepath),
        max_number_of_rounds_to_run(max_number_of_rounds_to_run),
        number_of_rounds_executed(0), previous_run_improved_consensus_seq(true),
        consensus_seq_already_seen() {}

    const std::string & get_polished_sequence() {
        run();
        return consensus_seq;
    }
};

#endif // PANDORA_RACON_H
