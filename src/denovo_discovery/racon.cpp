#include "denovo_discovery/racon.h"
#include <sstream>
#include "sequence.hpp"
#include "polisher.hpp"
#include "minimap.h"
#include "kseq.h"
#include "fastaq_handler.h"
#include <zlib.h>


std::string Racon::run_racon_core(
    const std::string &consensus_seq,
    const std::string &reads_filepath,
    const std::string &paf_filepath,
    const std::string &polished_filepath
    ) {
    BOOST_LOG_TRIVIAL(info) << "Start Running Racon";
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

    BOOST_LOG_TRIVIAL(info) << "End Running Racon";
    return polished_consensus_seq;
}

// Runs minimap2
// @return the number of reads mapped (entries in the PAF file)
uint32_t run_minimap2_core(bool illumina, const std::string &locus_consensus_filepath,
    const std::string &locus_reads_filepath, const std::string &paf_filepath) {
    BOOST_LOG_TRIVIAL(info) << "Start Running Minimap2";

    // setup minimap2 options
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    int n_threads = 1;
    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);

    if (illumina) {
        mm_set_opt("sr", &iopt, &mopt);
    } else {
        mm_set_opt("map-ont", &iopt, &mopt);
    }

    // TODO: set k?
    // iopt.k = kmer_prg->k;
    mopt.flag |= MM_F_CIGAR; // perform alignment

    // open index reader
    mm_idx_reader_t *r = mm_idx_reader_open(locus_consensus_filepath.c_str(), &iopt, 0);
    if (!r) {
        fatal_error("Could not open mm_idx_reader_t for ", locus_consensus_filepath);
    }

    // open query seq file
    gzFile query_seq_fd = gzopen(locus_reads_filepath.c_str(), "r");
    if (!query_seq_fd) {
        fatal_error(
            "Could not open minimap2 query sequence file for ", locus_reads_filepath);
    }

    // minimap2 the reads
    kseq_t *ks = kseq_init(query_seq_fd);
    mm_idx_t *mi;
    FILE* paf_filehandler = fopen(paf_filepath.c_str(), "w");
    if (!paf_filehandler) {
        fatal_error("Could not open PAF file: ", paf_filepath);
    }

    uint32_t number_of_reads_mapped = 0;
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
                    fprintf(paf_filehandler, "%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
                    fprintf(paf_filehandler, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", mi->seq[r->rid].name, mi->seq[r->rid].len,
                        r->rs, r->re, r->mlen, r->blen, r->mapq);
                    free(r->p);
                    ++number_of_reads_mapped;
                }
                free(reg);
            }
            mm_tbuf_destroy(tbuf);
            mm_idx_destroy(mi);
    }
    fclose(paf_filehandler);
    mm_idx_reader_close(r); // close the index reader
    kseq_destroy(ks); // close the query file
    gzclose(query_seq_fd);

    BOOST_LOG_TRIVIAL(info) << "End Running Minimap2";

    return number_of_reads_mapped;
}

bool Racon::run_another_round() {
    // builds a mem_fd with the consensus seq
    const std::string consensus_record = ">consensus\n" + consensus_seq;
    // TODO: use mem_fd when we can
//    std::pair<int, std::string> consensus_fd_and_filepath =
//        build_memfd(consensus_record, ".fa");
//    std::string locus_consensus_filepath = consensus_fd_and_filepath.second;
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
    uint32_t number_of_reads_mapped = run_minimap2_core(illumina,
        locus_consensus_filepath, locus_reads_filepath, paf_filepath);

    const bool paf_file_is_empty = number_of_reads_mapped == 0;

    // polish the sequence
    std::string polished_consensus_seq;
    if (paf_file_is_empty) {
        // racon errors out in this case, and it would not polish anything anyway
        polished_consensus_seq = consensus_seq;
    } else {
        // run racon to correct the ML seq
        polished_consensus_seq = run_racon_core(consensus_seq, locus_reads_filepath,
            paf_filepath, locus_consensus_filepath);
    }

    const bool we_already_saw_this_consensus_seq =
        std::find(consensus_seq_already_seen.begin(), consensus_seq_already_seen.end(),
            polished_consensus_seq) != consensus_seq_already_seen.end();
    const bool racon_improved_the_consensus_seq = not we_already_saw_this_consensus_seq;

    if (not we_already_saw_this_consensus_seq) {
        consensus_seq_already_seen.push_back(polished_consensus_seq);
    }

    consensus_seq = polished_consensus_seq;
    number_of_rounds_executed++;

    // fs::remove(locus_consensus_filepath);
    // fs::remove(paf_filepath);

    BOOST_LOG_TRIVIAL(debug) << "Ran racon " << number_of_rounds_executed << " times";

    return racon_improved_the_consensus_seq;
}

void Racon::run() {
    while (number_of_rounds_executed < max_number_of_rounds_to_run &&
           previous_run_improved_consensus_seq) {
        previous_run_improved_consensus_seq = run_another_round();
    }
}

