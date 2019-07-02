#include "denovo_discovery/denovo_discovery.h"


void DenovoDiscovery::find_paths_through_candidate_region(CandidateRegion &candidate_region) {
    const auto read_covg { candidate_region.pileup.size() };
    const auto length_of_candidate_sequence { candidate_region.max_likelihood_sequence.length() };
    const double expected_kmer_covg { calculate_kmer_coverage(read_covg, length_of_candidate_sequence) };

    BOOST_LOG_TRIVIAL(debug) << "Running local assembly for: " << candidate_region.get_name() << " - interval ["
                             << candidate_region.get_interval() << "]";

    if (candidate_region.pileup.empty()) {
        BOOST_LOG_TRIVIAL(debug) << "No sequences to assemble. Skipping local assembly.";
        return;
    }

    const auto max_path_length { length_of_candidate_sequence + 50 };

    if (kmer_size > max_path_length) {
        BOOST_LOG_TRIVIAL(debug) << "Kmer size " << kmer_size << " is greater than the maximum path length "
                                 << max_path_length << ". Skipping local assembly.";
        return;
    }

    LocalAssemblyGraph graph;

    try {
        Graph gatb_graph = LocalAssemblyGraph::create(new BankStrings(candidate_region.pileup),
                                                      "-kmer-size %d -abundance-min %d -verbose 0", kmer_size,
                                                      min_covg_for_node_in_assembly_graph);
        if (clean_assembly_graph) {
            clean(gatb_graph);
        }
        graph = gatb_graph; //TODO: use move constructor

    } catch (gatb::core::system::Exception &error) {
        BOOST_LOG_TRIVIAL(debug) << "Couldn't create GATB graph." << "\n\tEXCEPTION: " << error.getMessage();
        remove_graph_file();
        return;
    }

    Node start_node, end_node;
    bool start_found { false };
    bool end_found { false };
    const auto start_kmers { generate_start_kmers(candidate_region.max_likelihood_sequence, kmer_size, kmer_size) };
    const auto end_kmers { generate_end_kmers(candidate_region.max_likelihood_sequence, kmer_size, kmer_size) };

    for (uint_least16_t start_idx = 0; start_idx < start_kmers.size(); start_idx++) {
        const auto &current_start_kmer { start_kmers[start_idx] };
        std::tie(start_node, start_found) = graph.get_node(current_start_kmer);

        if (not start_found) {
            continue;
        }

        for (uint_least16_t end_idx = 0; end_idx < end_kmers.size(); end_idx++) {
            const auto &current_end_kmer { end_kmers[end_idx] };

            std::tie(end_node, end_found) = graph.get_node(current_end_kmer);

            if (end_found) {
                auto tree_of_nodes_visited_during_dfs { graph.depth_first_search_from(start_node) };
                auto denovo_paths {
                        graph.get_paths_between(current_start_kmer, current_end_kmer, tree_of_nodes_visited_during_dfs,
                                                max_path_length, expected_kmer_covg) };
                candidate_region.denovo_paths.insert(candidate_region.denovo_paths.begin(), denovo_paths.begin(),
                                                     denovo_paths.end());
                if (not candidate_region.denovo_paths.empty()) {
                    // add flank to each sequence - the whole sequence from the start to the end of the gene
                    for (auto &current_path : candidate_region.denovo_paths) {
                        const auto start_kmer_offset { candidate_region.max_likelihood_sequence.substr(0, start_idx) };
                        const auto end_kmer_offset { candidate_region.max_likelihood_sequence
                                                                     .substr(candidate_region.max_likelihood_sequence
                                                                                             .length() - end_idx) };
                        current_path = std::string(candidate_region.left_flanking_sequence).append(start_kmer_offset)
                                                                                           .append(current_path)
                                                                                           .append(end_kmer_offset)
                                                                                           .append(candidate_region
                                                                                                           .right_flanking_sequence);
                    }
                }

                remove_graph_file();
                return;
            }
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "Could not find any combination of start and end k-mers. Skipping local assembly for "
                             << candidate_region.get_name();
    remove_graph_file();

}


DenovoDiscovery::DenovoDiscovery(const uint_least8_t &kmer_size, const double &read_error_rate)
        : kmer_size { kmer_size }, read_error_rate { read_error_rate } {}


double DenovoDiscovery::calculate_kmer_coverage(const uint32_t &read_covg, const uint32_t &ref_length) const {
    if (ref_length == 0) {
        throw std::invalid_argument("ref_length should be greater than 0.");
    } else if (kmer_size == 0) {
        throw std::invalid_argument("K should be greater than 0.");
    } else if (read_error_rate < 0) {
        throw std::invalid_argument("error_rate should not be a negative value.");
    }

    const auto numerator { read_covg * (ref_length - kmer_size + 1) * std::pow(1 - read_error_rate, kmer_size) };
    const auto denominator { ref_length };
    return numerator / denominator;
}

