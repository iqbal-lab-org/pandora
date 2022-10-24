#include "sam_file.h"
#include "minihits.h"
#include "minihit.h"
#include "version.h"
#include "globals.h"

std::string SAMFile::get_header() const {
    std::stringstream ss;

    for (const auto &prg_name_and_length : prg_name_to_length) {
        ss << "@SQ\t"
           << "SN:" << prg_name_and_length.first << "\t"
           << "LN:" << prg_name_and_length.second << "\n";
    }
    ss << "@PG\tID:pandora\tPN:pandora\tVN:" << PANDORA_VERSION
       << "\tCL: " << PandoraGlobals::command_line << "\n";
    ss << "@CO\tThe reference length (in @SQ header lines) and the POS field refer to "
          "the string representation of the PRGs\n";
    ss << "@CO\tLF: left flank sequence, the sequence before the first "
          "mapped kmer, soft-clipped, max " << flank_size << " bps\n";
    ss << "@CO\tRF: right flank sequence, the sequence after the last "
          "mapped kmer, soft-clipped, max " << flank_size << " bps\n";
    ss << "@CO\tMP: number of minimizer matches on the plus strand\n";
    ss << "@CO\tMM: number of minimizer matches on the minus strand\n";
    ss << "@CO\tPP: Prg Paths of the cluster of hits: the PRG path of each "
          "hit in considered cluster of hits\n";
    ss << "@CO\tNM: Total number of mismatches in the quasi-alignment\n";
    ss << "@CO\tAS: Alignment score (number of matches)\n";
    ss << "@CO\tnn: Number of ambiguous bases in the quasi-alignment\n";
    ss << "@CO\tcm: Number of minimizers in the quasi-alignment\n";

    return ss.str();
}

void SAMFile::create_final_sam_file() {
    GenericFile sam_file(filepath);
    sam_file << get_header();
    std::ifstream sam_records_fh;
    open_file_for_reading(tmp_filepath.string(), sam_records_fh);
    sam_file << sam_records_fh.rdbuf();
    sam_records_fh.close();
}

SAMFile::~SAMFile() {
    try {
        delete tmp_sam_file;  // closes tmp_sam_file
        create_final_sam_file();
        fs::remove(tmp_filepath);
    } catch (...) {
        fatal_error("Error creating SAM file: ", filepath.string());
    };
}

std::vector<bool> SAMFile::get_mapped_positions_bitset(const Seq &seq, const MinimizerHits &cluster) const {
    std::vector<bool> mapped_positions_bitset(seq.full_seq.size(), false);
    for (const MinimizerHitPtr &hit : cluster) {
        const uint32_t read_start_position = hit->get_read_start_position();
        const uint32_t read_end_position = hit->get_read_start_position()+hit->get_prg_path().length();
        for (uint32_t pos = read_start_position; pos < read_end_position; ++pos) {
            mapped_positions_bitset[pos] = true;
        }
    }
    return mapped_positions_bitset;
}

uint32_t get_first_mapped_position(const std::vector<bool> &mapped_positions_bitset) {
    return std::find(mapped_positions_bitset.begin(), mapped_positions_bitset.end(), 1) -
           mapped_positions_bitset.begin();
}
uint32_t get_last_mapped_position(const std::vector<bool> &mapped_positions_bitset) {
    return mapped_positions_bitset.rend() - std::find(
        mapped_positions_bitset.rbegin(), mapped_positions_bitset.rend(), 1);
}

Cigar SAMFile::get_cigar(const std::vector<bool> &mapped_positions_bitset) const {
    Cigar cigar;
    const uint32_t first_mapped_position = get_first_mapped_position(mapped_positions_bitset);
    const uint32_t last_mapped_position = get_last_mapped_position(mapped_positions_bitset);

    if (first_mapped_position > 0) {
        cigar.add_entry('H', first_mapped_position);
    }

    uint32_t consecutive_equal_bits_length = 0;
    bool consecutive_bit_value = true;
    for (uint32_t pos = first_mapped_position; pos < last_mapped_position; ++pos) {
        const bool current_bit = mapped_positions_bitset[pos];
        const bool extend_window = current_bit == consecutive_bit_value;
        if (extend_window) {
            ++consecutive_equal_bits_length;
        }else {
            const char cigar_char = consecutive_bit_value ? '=' : 'X';
            cigar.add_entry(cigar_char, consecutive_equal_bits_length);
            consecutive_equal_bits_length = 1;
            consecutive_bit_value = !consecutive_bit_value;
        }
    }
    {
        const char cigar_char = consecutive_bit_value ? '=' : 'X';
        cigar.add_entry(cigar_char, consecutive_equal_bits_length);
    }

    const uint32_t number_of_bases_after_last = mapped_positions_bitset.size() - last_mapped_position;
    if (number_of_bases_after_last > 0) {
        cigar.add_entry('H', number_of_bases_after_last);
    }
    return cigar;
}

std::string SAMFile::get_segment_sequence(const Seq &seq,
                                          const std::vector<bool> &mapped_positions_bitset) const {
    const uint32_t first_mapped_position = get_first_mapped_position(mapped_positions_bitset);
    const uint32_t last_mapped_position = get_last_mapped_position(mapped_positions_bitset);
    const uint32_t number_of_bases_between_first_and_last = last_mapped_position - first_mapped_position;
    return seq.full_seq.substr(first_mapped_position, number_of_bases_between_first_and_last);
}

void SAMFile::write_sam_record_from_hit_cluster(
    const Seq &seq, const MinimizerHitClusters &clusters) {
    bool at_least_a_single_mapping_was_output = false;
    for (const MinimizerHits &cluster : clusters) {
        if (cluster.empty()) {
            continue;
        }
        at_least_a_single_mapping_was_output = true;

        uint32_t plus_strand_count, minus_strand_count;
        std::tie(plus_strand_count, minus_strand_count) = cluster.get_strand_counts();
        const bool is_plus_strand = plus_strand_count >= minus_strand_count;
        const uint32_t flag = is_plus_strand ? 0 : 16;

        const MinimizerHitPtr first_hit = *(cluster.begin());
        const std::string &prg_name = prgs[first_hit->get_prg_id()]->name;

        const uint32_t alignment_start = first_hit->get_prg_path().begin()->start + 1;

        const std::vector<bool> mapped_positions_bitset =
            get_mapped_positions_bitset(seq, cluster);
        const Cigar cigar = get_cigar(mapped_positions_bitset);

        std::string segment_sequence = get_segment_sequence(seq, mapped_positions_bitset);

        const uint32_t first_mapped_pos = get_first_mapped_position(mapped_positions_bitset);
        const uint32_t left_flank_start = first_mapped_pos <= flank_size ? 0 :
                                      first_mapped_pos - flank_size;
        const uint32_t left_flank_length = first_mapped_pos - left_flank_start;
        std::string left_flank = seq.full_seq.substr(left_flank_start, left_flank_length);

        const uint32_t last_mapped_pos = get_last_mapped_position(mapped_positions_bitset);
        const uint32_t right_flank_start = last_mapped_pos;
        const uint32_t right_flank_length =
            (right_flank_start + flank_size) <= seq.full_seq.size() ? flank_size :
            seq.full_seq.size() - right_flank_start;
        std::string right_flank = seq.full_seq.substr(right_flank_start, right_flank_length);

        const bool is_reverse_complemented = not is_plus_strand;
        if (is_reverse_complemented) {
            segment_sequence = rev_complement(segment_sequence);
            left_flank = rev_complement(left_flank);
            right_flank = rev_complement(right_flank);
            std::swap(left_flank, right_flank);
        }

        std::stringstream cluster_of_hits_prg_paths_ss;
        for (const MinimizerHitPtr &hit : cluster) {
            cluster_of_hits_prg_paths_ss << hit->get_prg_path() << "->";
        }

        auto number_ambiguous_bases = std::count_if(
            segment_sequence.begin(), segment_sequence.end(),
            [](char base) -> bool {
                switch (base) {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                    return true;
                default:
                    return false;
                }
            }
        );

        uint32_t number_of_mismatches = cigar.number_of_mismatches();
        uint32_t alignment_score = segment_sequence.size() - number_of_mismatches;

        (*tmp_sam_file)  << seq.name << "\t"
                         << flag << "\t"
                         << prg_name << "\t"
                         << alignment_start << "\t"
                         << "255\t"
                         << cigar << "\t"
                         << "*\t0\t0\t"
                         << segment_sequence << "\t"
                         << "*\t"
                         << "LF:Z:" << left_flank << "\t"
                         << "RF:Z:" << right_flank << "\t"
                         << "MP:i:" << plus_strand_count << "\t"
                         << "MM:i:" << minus_strand_count << "\t"
                         << "PP:Z:" << cluster_of_hits_prg_paths_ss.str() << "\t"
                         << "NM:i:" << number_of_mismatches << "\t"
                         << "AS:i:" << alignment_score << "\t"
                         << "nn:i:" << number_ambiguous_bases << "\t"
                         << "cm:i:" << cluster.size() << "\n";

        const uint32_t prg_length = prgs[first_hit->get_prg_id()]->seq.size();
        prg_name_to_length[prg_name] = prg_length;
    }

    if (!at_least_a_single_mapping_was_output) {
        (*tmp_sam_file) << seq.name << "\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
    }
}

