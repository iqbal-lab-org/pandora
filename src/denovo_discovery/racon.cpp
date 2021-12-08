#include "denovo_discovery/racon.h"
#include <sstream>


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
    minimap_ss << "minimap2 -t 1 -x map-ont -o " << paf_filepath << " " <<
        locus_consensus_filepath << " " << locus_reads_filepath;
    const std::string minimap_str = minimap_ss.str();
    exec(minimap_str.c_str());

    // run racon to correct the ML seq
    const std::string racon_out {
        (denovo_outdir / (locus + ".racon.out")).string() };
    std::stringstream racon_ss;
    racon_ss << "racon -u --no-trimming -t 1 " << locus_reads_filepath << " " << paf_filepath <<
        " " << locus_consensus_filepath << " > " << racon_out;
    const std::string racon_str = racon_ss.str();
    exec(racon_str.c_str());

    // get the polished consensus seq
    const std::vector<std::string> polished_record
        = get_vector_of_strings_from_file(racon_out);

    const bool only_one_denovo_sequence = polished_record.size()==2;
    if (!only_one_denovo_sequence) {
        fatal_error("Racon outputted no sequences or more than one sequence");
    }

    const std::string &polished_consensus_seq = polished_record[1];

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
    fs::remove(racon_out);

    BOOST_LOG_TRIVIAL(debug) << "Ran racon " << number_of_rounds_executed << " times";

    return racon_improved_the_consensus_seq;
}

void Racon::run() {
    while (number_of_rounds_executed < max_number_of_rounds_to_run &&
           previous_run_improved_consensus_seq) {
        previous_run_improved_consensus_seq = run_another_round();
    }
}

