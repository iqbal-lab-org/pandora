#include "sam_file.h"
#include "minihits.h"
#include "minihit.h"
#include "version.h"

SAMFile::SAMFile(const fs::path &filepath,
                 const std::vector<std::shared_ptr<LocalPRG>>& prgs) :
    GenericFile(filepath), prgs(prgs) {
    file_handler << "@PG\tID:pandora\tPN:pandora\tVN:" << PANDORA_VERSION
                 << "\tCL:TODO" << std::endl;
}

std::vector<bool> SAMFile::get_mapped_positions_bitset(const Seq &seq, const Hits &cluster) const {
    std::vector<bool> mapped_positions_bitset(seq.seq.size(), 0);
    for (const MinimizerHitPtr &hit : cluster) {
        const uint32_t read_start_position = hit->get_read_start_position();
        const uint32_t read_end_position = hit->get_read_start_position()+hit->get_prg_path().length();
        for (uint32_t pos = read_start_position; pos < read_end_position; ++pos) {
            mapped_positions_bitset[pos] = 1;
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

std::string SAMFile::get_cigar(const std::vector<bool> &mapped_positions_bitset) const {
    std::stringstream cigar_ss;

    const uint32_t first_mapped_position = get_first_mapped_position(mapped_positions_bitset);
    const uint32_t last_mapped_position = get_last_mapped_position(mapped_positions_bitset);

    if (first_mapped_position > 0) {
        cigar_ss << first_mapped_position << "S";
    }

    uint32_t consecutive_equal_bits_length = 0;
    bool consecutive_bit_value = true;
    for (uint32_t pos = first_mapped_position; pos < last_mapped_position; ++pos) {
        const bool current_bit = mapped_positions_bitset[pos];
        const bool extend_window = current_bit == consecutive_bit_value;
        if (extend_window) {
            ++consecutive_equal_bits_length;
        }else {
            const std::string cigar_char = consecutive_bit_value ? "=" : "X";
            cigar_ss << consecutive_equal_bits_length << cigar_char;
            consecutive_equal_bits_length = 1;
            consecutive_bit_value = !consecutive_bit_value;
        }
    }
    {
        const std::string cigar_char = consecutive_bit_value ? "=" : "X";
        cigar_ss << consecutive_equal_bits_length << cigar_char;
    }

    const uint32_t number_of_bases_after_last = mapped_positions_bitset.size() - last_mapped_position;
    if (number_of_bases_after_last > 0) {
        cigar_ss << number_of_bases_after_last << "S";
    }
    return cigar_ss.str();
}

std::string SAMFile::get_segment_sequence(const Seq &seq,
                                          const std::vector<bool> &mapped_positions_bitset) const {
    const uint32_t first_mapped_position = get_first_mapped_position(mapped_positions_bitset);
    const uint32_t last_mapped_position = get_last_mapped_position(mapped_positions_bitset);
    const uint32_t number_of_bases_between_first_and_last = last_mapped_position - first_mapped_position;
    return seq.seq.substr(first_mapped_position, number_of_bases_between_first_and_last);
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
        file_handler << seq.name  << "\t"
                     << "0\t"
                     << prgs[first_hit->get_prg_id()]->name << "\t"
                     << first_hit->get_prg_path() << "\t"
                     << "254\t" << "\t"
                     << "*\t0\t0\t"
                     << seq.seq << "\t"
                     << "*\t" << std::endl;
    }

    if (!at_least_a_single_mapping_was_output) {
        file_handler << seq.name << "\t" << "4\t*\t*\t0\t*\t*\t0\t0\t*\t*\t" << std::endl;
    }
}

