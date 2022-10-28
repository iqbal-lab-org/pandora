#include "denovo_discovery/racon.h"
#include <sstream>
#include "sequence.hpp"
#include "polisher.hpp"
#include "minimap.h"
#include "kseq.h"
#include "fastaq_handler.h"
#include <zlib.h>


std::string Racon::run_racon_core(
    const std::string &sequence_to_be_corrected,
    const std::string &reads_filepath,
    const std::string &paf_filepath,
    const std::string &target_seq_filepath) {
    uint32_t window_length = sequence_to_be_corrected.size() + Racon::racon_window_padding;
    const static double quality_threshold = 10.0;
    const static double error_threshold = 0.3;
    const static bool trim = false;

    const static int8_t match = 3;
    const static int8_t mismatch = -5;
    const static int8_t gap = -4;
    const static uint32_t type = 0;

    const static bool drop_unpolished_sequences = true;
    const static uint32_t num_threads = 1;

    const static uint32_t cudapoa_batches = 0;
    const static uint32_t cudaaligner_batches = 0;
    const static uint32_t cudaaligner_band_width = 0;
    const static bool cuda_banded_alignment = false;

    bool polishing_was_successful;
    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    try {
        std::unique_ptr<racon::Polisher> polisher
            = racon::createPolisher(reads_filepath, paf_filepath, target_seq_filepath,
                type == 0 ? racon::PolisherType::kC : racon::PolisherType::kF,
                window_length, quality_threshold, error_threshold, trim, match,
                mismatch, gap, num_threads, cudapoa_batches, cuda_banded_alignment,
                cudaaligner_batches, cudaaligner_band_width);

        polisher->initialize();

        polisher->polish(polished_sequences, drop_unpolished_sequences);

        polishing_was_successful = polished_sequences.size()>=1;
    }catch (const racon::EmptyOverlapSetError &empty_overlap_set_error) {
        polishing_was_successful = false;
    }

    std::string polished_consensus_seq;
    if (polishing_was_successful) {
        polished_consensus_seq = polished_sequences[0]->data();
    } else {
        polished_consensus_seq = sequence_to_be_corrected;
    }

    return polished_consensus_seq;
}

// Runs minimap2
// Returns PAF string
std::string Racon::run_minimap2_core(
    const std::string &query_filepath,
    const std::string &ref_filepath,
    const bool illumina, const uint32_t kmer_size
    ) {
    // setup minimap2 options
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    const static int n_threads = 1;
    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);

    if (illumina) {
        mm_set_opt("sr", &iopt, &mopt);
    } else {
        mm_set_opt("map-ont", &iopt, &mopt);
    }

    iopt.k = kmer_size;
    mopt.flag |= MM_F_CIGAR; // perform alignment - this is actually required

    // open index reader
    mm_idx_reader_t *r = mm_idx_reader_open(ref_filepath.c_str(), &iopt, 0);
    if (!r) {
        fatal_error("Could not open mm_idx_reader_t for ", ref_filepath);
    }

    // open query seq file
    gzFile query_seq_fd = gzopen(query_filepath.c_str(), "r");
    if (!query_seq_fd) {
        fatal_error(
            "Could not open minimap2 query sequence file for ", query_filepath);
    }

    // minimap2 the reads
    kseq_t *ks = kseq_init(query_seq_fd);
    mm_idx_t *mi;

    std::stringstream paf_ss;
    while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
            mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
            mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
            gzrewind(query_seq_fd);
            kseq_rewind(ks);
            while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
                mm_reg1_t *reg;
                int j, i, n_reg;
                reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
                for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                    mm_reg1_t *r = &reg[j];
                    assert(r->p); // with MM_F_CIGAR, this should not be NULL
                    paf_ss << ks->name.s << "\t" << ks->seq.l << "\t" << r->qs << "\t"
                           << r->qe << "\t" << "+-"[r->rev] << "\t"
                           << mi->seq[r->rid].name << "\t" << mi->seq[r->rid].len << "\t"
                           << r->rs << "\t" << r->re << "\t" << r->mlen << "\t"
                           << r->blen << "\t" << r->mapq << "\n";
                    free(r->p);
                }
                free(reg);
            }
            mm_tbuf_destroy(tbuf);
            mm_idx_destroy(mi);
    }
    mm_idx_reader_close(r); // close the index reader
    kseq_destroy(ks); // close the query file
    gzclose(query_seq_fd);

    return paf_ss.str();
}

bool Racon::run_another_round() {
    const std::string consensus_record = ">consensus_" + locus + "\n" + consensus_seq;

    // builds a mem_fd with the consensus seq
    const std::pair<int, std::string> locus_consensus_fd_and_filepath = build_memfd(consensus_record);

    // run minimap to get overlaps
    std::string paf = run_minimap2_core(
        this->locus_reads_filepath,
        locus_consensus_fd_and_filepath.second,
        this->illumina, this->kmer_size);
    const bool paf_file_is_empty = paf.length() == 0;

    // polish the sequence
    std::string polished_consensus_seq;
    if (paf_file_is_empty) {
        // racon errors out in this case, and it would not polish anything anyway
        polished_consensus_seq = consensus_seq;
    } else {
        // builds memfd with paf file
        const std::pair<int, std::string> paf_fd_and_filepath = build_memfd(paf);

        // run racon to correct the ML seq
        polished_consensus_seq = run_racon_core(consensus_seq,
            locus_reads_filepath, paf_fd_and_filepath.second,
            locus_consensus_fd_and_filepath.second);

        close(paf_fd_and_filepath.first);
    }
    close(locus_consensus_fd_and_filepath.first);

    const bool we_already_saw_this_consensus_seq =
        std::find(consensus_seq_already_seen.begin(), consensus_seq_already_seen.end(),
            polished_consensus_seq) != consensus_seq_already_seen.end();
    if (not we_already_saw_this_consensus_seq) {
        consensus_seq_already_seen.push_back(polished_consensus_seq);
    }

    consensus_seq = polished_consensus_seq;
    number_of_rounds_executed++;

    BOOST_LOG_TRIVIAL(debug) << "Ran racon " << number_of_rounds_executed << " times";

    const bool racon_improved_the_consensus_seq = not we_already_saw_this_consensus_seq;
    return racon_improved_the_consensus_seq;
}

void Racon::run() {
    if (keep_extra_debugging_files) {
        // build a file with Racon input
        const std::string consensus_record = ">consensus_" + locus + "\n" + consensus_seq;
        build_file((denovo_outdir / (locus+".consensus.fa")).string(), consensus_record);
    }

    while (number_of_rounds_executed < max_number_of_rounds_to_run &&
           previous_run_improved_consensus_seq) {
        previous_run_improved_consensus_seq = run_another_round();
    }

    if (keep_extra_debugging_files) {
        // build a file on disk with the polished consensus
        std::string polished_consensus_filepath;
        {
            std::stringstream ss;
            ss << locus << ".consensus.racon." << number_of_rounds_executed << "_rounds.fa";
            polished_consensus_filepath = ss.str();
        }
        polished_consensus_filepath = (denovo_outdir / polished_consensus_filepath).string();
        build_file(polished_consensus_filepath,
            std::string(">polished_consensus_" + locus + "\n" + consensus_seq));
    }
}
