#include "denovo_discovery/denovo_discovery.h"


fs::path get_discovered_paths_fname(const GeneIntervalInfo &gene_slice, const uint32_t &local_assembly_kmer_size) {
    return fs::path(
            std::string(gene_slice.pnode->get_name()).append(".").append(std::to_string(gene_slice.interval.start))
                                                     .append("-").append(std::to_string(gene_slice.interval.get_end()))
                                                     .append("_local_assembly_K")
                                                     .append(std::to_string(local_assembly_kmer_size)).append(".fa"));
}

void
denovo_discovery::find_candidates(const std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &candidate_coordinates,
                                  const std::string &reads_filepath, const fs::path &output_directory,
                                  const double &error_rate, const uint32_t &local_assembly_kmer_size) {

    const auto pileups { denovo_discovery::collect_read_pileups(candidate_coordinates, reads_filepath) };
    for (const auto &pileup : pileups) {
        const auto &sequences { pileup.second };
        const auto &coordinate_info { pileup.first };
        const auto &interval_sequence { coordinate_info.seq };

        const auto current_denovo_paths_filename {
                get_discovered_paths_fname(coordinate_info, local_assembly_kmer_size) };
        const auto denovo_paths_dir { output_directory / "denovo_paths" };
        fs::create_directories(denovo_paths_dir);
        const auto discovered_paths_filepath { denovo_paths_dir / current_denovo_paths_filename };

        const auto read_covg { sequences.size() };
        const auto length_of_interval_slice { interval_sequence.length() };
        const double expected_kmer_covg {
                denovo_discovery::calculate_kmer_coverage(read_covg, length_of_interval_slice, local_assembly_kmer_size,
                                                          error_rate) };

        BOOST_LOG_TRIVIAL(debug) << "Running local assembly for: " << coordinate_info.pnode->get_name()
                                 << " - interval [" << coordinate_info.interval.start << ", "
                                 << coordinate_info.interval.get_end() << "]";

        local_assembly(sequences, interval_sequence, coordinate_info.flank_seq_left, coordinate_info.flank_seq_right,
                       discovered_paths_filepath, local_assembly_kmer_size, interval_sequence.length() * 2,
                       expected_kmer_covg);
    }
}

double
denovo_discovery::calculate_kmer_coverage(const uint32_t &read_covg, const uint32_t &ref_length, const uint32_t &k,
                                          const double &error_rate) {
    if (ref_length == 0) {
        throw std::invalid_argument("ref_length should be greater than 0.");
    } else if (k == 0) {
        throw std::invalid_argument("K should be greater than 0.");
    } else if (error_rate < 0) {
        throw std::invalid_argument("error_rate should not be a negative value.");
    }

    const auto numerator { read_covg * (ref_length - k + 1) * std::pow(1 - error_rate, k) };
    const auto denominator { ref_length };
    return numerator / denominator;
}