#include "sam_file.h"
#include "minihits.h"
#include "minihit.h"
#include "version.h"
#include "globals.h"

SAMFile::SAMFile(const fs::path &filepath,
                 const std::vector<std::shared_ptr<LocalPRG>>& prgs,
                 const uint32_t flank_size) :
    GenericFile(filepath), prgs(prgs), flank_size(flank_size) {
    (*this) << "@PG\tID:pandora\tPN:pandora\tVN:" << PANDORA_VERSION
            << "\tCL: " << PandoraGlobals::command_line << "\n";
    (*this) << "@CO\tLF: left flank sequence, the sequence before the first "
               "mapped kmer, soft-clipped, max " << flank_size << " bps\n";
    (*this) << "@CO\tRF: right flank sequence, the sequence after the last "
               "mapped kmer, soft-clipped, max " << flank_size << " bps\n";
    (*this) << "@CO\tPPCH: Prg Paths of the Cluster of Hits: the PRG path of each "
               "hit in considered cluster of hits\n";
    (*this) << "@CO\tNM: Total number of mismatches in the quasi-alignment\n";
    (*this) << "@CO\tAS: Alignment score (number of matches)\n";
    (*this) << "@CO\tnn: Number of ambiguous bases in the quasi-alignment\n";
    (*this) << "@CO\tcm: Number of minimizers in the quasi-alignment\n";
}

std::vector<bool> SAMFile::get_mapped_positions_bitset(const Seq &seq, const Hits &cluster) const {
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
        cigar.add_entry('S', first_mapped_position);
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
        cigar.add_entry('S', number_of_bases_after_last);
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
    for (const Hits &cluster : clusters) {
        if (cluster.empty()) {
            continue;
        }
        at_least_a_single_mapping_was_output = true;

        const MinimizerHitPtr first_hit = *(cluster.begin());
        const std::vector<bool> mapped_positions_bitset =
            get_mapped_positions_bitset(seq, cluster);
        const Cigar cigar = get_cigar(mapped_positions_bitset);
        const std::string segment_sequence = get_segment_sequence(seq, mapped_positions_bitset);

        const uint32_t first_mapped_pos = get_first_mapped_position(mapped_positions_bitset);
        const uint32_t left_flank_start = first_mapped_pos <= flank_size ? 0 :
                                      first_mapped_pos - flank_size;
        const uint32_t left_flank_length = first_mapped_pos - left_flank_start;
        const std::string left_flank = seq.full_seq.substr(left_flank_start, left_flank_length);

        const uint32_t last_mapped_pos = get_last_mapped_position(mapped_positions_bitset);
        const uint32_t right_flank_start = last_mapped_pos;
        const uint32_t right_flank_length =
            (right_flank_start + flank_size) <= seq.full_seq.size() ? flank_size :
            seq.full_seq.size() - right_flank_start;
        const std::string right_flank = seq.full_seq.substr(right_flank_start, right_flank_length);

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

        (*this)  << seq.name << "[" << first_mapped_pos << ":" << last_mapped_pos << "]\t"
                 << "0\t"
                 << prgs[first_hit->get_prg_id()]->name << "\t"
                 << first_hit->get_prg_path() << "\t"
                 << "255\t"
                 << cigar << "\t"
                 << "*\t0\t0\t"
                 << segment_sequence << "\t"
                 << "*\t"
                 << "LF:Z:" << left_flank << "\t"
                 << "RF:Z:" << right_flank << "\t"
                 << "PPCH:Z:" << cluster_of_hits_prg_paths_ss.str() << "\t"
                 << "NM:i:" << number_of_mismatches << "\t"
                 << "AS:i:" << alignment_score << "\t"
                 << "nn:i:" << number_ambiguous_bases << "\t"
                 << "cm:i:" << cluster.size() << "\n";
    }

    if (!at_least_a_single_mapping_was_output) {
        (*this) << seq.name << "\t" << "4\t*\t*\t0\t*\t*\t0\t0\t*\t*\t" << "\n";
    }
}

