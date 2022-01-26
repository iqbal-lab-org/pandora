#include "denovo_discovery/racon.h"
#include <sstream>
#include "sequence.hpp"
#include "polisher.hpp"

std::string Racon::run_racon_core(
    const std::string &consensus_seq,
    const std::string &reads_filepath,
    const std::string &paf_filepath,
    const std::string &polished_filepath
    ) {
    uint32_t window_length = consensus_seq.size();
    double quality_threshold = 10.0;
    double error_threshold = 0.3;
    bool trim = false;

    int8_t match = 3;
    int8_t mismatch = -5;
    int8_t gap = -4;
    uint32_t type = 0;

    bool drop_unpolished_sequences = true;
    uint32_t num_threads = 1;

    uint32_t cudapoa_batches = 0;
    uint32_t cudaaligner_batches = 0;
    uint32_t cudaaligner_band_width = 0;
    bool cuda_banded_alignment = false;

    std::unique_ptr<racon::Polisher> polisher = racon::createPolisher(
        reads_filepath, paf_filepath, polished_filepath,
        type == 0 ? racon::PolisherType::kC : racon::PolisherType::kF,
        window_length, quality_threshold, error_threshold, trim, match, mismatch, gap,
        num_threads, cudapoa_batches, cuda_banded_alignment, cudaaligner_batches,
        cudaaligner_band_width);

    polisher->initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;

    polisher->polish(polished_sequences, drop_unpolished_sequences);

    std::string polished_consensus_seq;
    const bool polishing_was_successful = polished_sequences.size()>=1;
    if (polishing_was_successful) {
        polished_consensus_seq = polished_sequences[0]->data();
    } else {
        polished_consensus_seq = consensus_seq;
    }

    return polished_consensus_seq;
}

bool Racon::run_another_round() {
    // builds a mem_fd with the consensus seq
    const std::string consensus_record = ">consensus\n" + consensus_seq;
    // TODO: use mem_fd back
    // std::pair<int, std::string> consensus_fd_and_filepath = build_memfd(consensus_record);
    std::string locus_consensus_filepath;
    {
        std::stringstream ss;
        ss << locus << ".consensus.fa";
        locus_consensus_filepath = ss.str();
    }
    locus_consensus_filepath = (denovo_outdir / locus_consensus_filepath).string();
    build_file(locus_consensus_filepath, consensus_record);


    // run minimap to get overlaps
    const std::string paf_filepath {
        (denovo_outdir / (locus + ".minimap2.out.paf")).string() };
    std::stringstream minimap_ss;
    minimap_ss << "minimap2 -t 1 -x " << (illumina ? "sr" : "map-ont")
               << " -o " << paf_filepath << " " << locus_consensus_filepath << " "
               << locus_reads_filepath;
    const std::string minimap_str = minimap_ss.str();
    exec(minimap_str.c_str());

    // run racon to correct the ML seq
    const std::string polished_consensus_seq = run_racon_core(
        consensus_seq, locus_reads_filepath, paf_filepath,locus_consensus_filepath);

    const bool we_already_saw_this_consensus_seq =
        std::find(consensus_seq_already_seen.begin(), consensus_seq_already_seen.end(),
            polished_consensus_seq) != consensus_seq_already_seen.end();
    const bool racon_improved_the_consensus_seq = not we_already_saw_this_consensus_seq;

    if (not we_already_saw_this_consensus_seq) {
        consensus_seq_already_seen.push_back(polished_consensus_seq);
    }

    consensus_seq = polished_consensus_seq;
    number_of_rounds_executed++;

    fs::remove(locus_consensus_filepath);
    fs::remove(paf_filepath);

    BOOST_LOG_TRIVIAL(debug) << "Ran racon " << number_of_rounds_executed << " times";

    return racon_improved_the_consensus_seq;
}

void Racon::run() {
    while (number_of_rounds_executed < max_number_of_rounds_to_run &&
           previous_run_improved_consensus_seq) {
        previous_run_improved_consensus_seq = run_another_round();
    }
}

