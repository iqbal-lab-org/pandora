#include "denovo_discovery/denovo_discovery.h"

void DenovoDiscovery::find_paths_through_candidate_region(
    CandidateRegion& candidate_region, const fs::path& temp_dir) const
{
    const auto read_covg { candidate_region.pileup.size() };
    const auto length_of_candidate_sequence {
        candidate_region.max_likelihood_sequence.length()
    };
    const double expected_kmer_covg { calculate_kmer_coverage(
        read_covg, length_of_candidate_sequence) };
    const fs::path GATB_graph_filepath(temp_dir / "GATB_graph");

    BOOST_LOG_TRIVIAL(debug) << "Running local assembly for: "
                             << candidate_region.get_name() << " - interval ["
                             << candidate_region.get_interval() << "]";

    if (candidate_region.pileup.empty()) {
        BOOST_LOG_TRIVIAL(debug)
            << "No sequences to assemble. Skipping local assembly.";
        return;
    }

    const auto max_path_length { length_of_candidate_sequence + max_insertion_size };

    if (kmer_size > max_path_length) {
        BOOST_LOG_TRIVIAL(debug)
            << "Kmer size " << kmer_size << " is greater than the maximum path length "
            << max_path_length << ". Skipping local assembly.";
        return;
    }

    LocalAssemblyGraph graph;

    try {
        const std::string GATB_graph_filepath_as_string = GATB_graph_filepath.string();
        Graph gatb_graph = LocalAssemblyGraph::create(
            new BankStrings(candidate_region.pileup),
            "-kmer-size %d -abundance-min %d -verbose 0 -nb-cores 1 -out %s", kmer_size,
            min_covg_for_node_in_assembly_graph, GATB_graph_filepath_as_string.c_str());
        if (this->clean_assembly_graph) {
            BOOST_LOG_TRIVIAL(debug) << "Cleaning local assembly graph...";
            clean(gatb_graph);
            BOOST_LOG_TRIVIAL(debug) << "Local assembly graph cleaned";
        }
        graph = gatb_graph; // TODO: use move constructor

    } catch (gatb::core::system::Exception& error) {
        BOOST_LOG_TRIVIAL(debug) << "Couldn't create GATB graph."
                                 << "\n\tEXCEPTION: " << error.getMessage();
        remove_graph_file(GATB_graph_filepath);
        return;
    }
    graph.set_max_nb_paths(this->max_nb_paths);

    Node start_node, end_node;
    bool start_found { false };
    bool end_found { false };
    const auto start_kmers { generate_start_kmers(
        candidate_region.max_likelihood_sequence, kmer_size, kmer_size) };
    const auto end_kmers { generate_end_kmers(
        candidate_region.max_likelihood_sequence, kmer_size, kmer_size) };

    for (uint_least16_t start_idx = 0; start_idx < start_kmers.size(); start_idx++) {
        const auto& current_start_kmer { start_kmers[start_idx] };
        BOOST_LOG_TRIVIAL(trace)
            << "Searching for start anchor kmer " << current_start_kmer;
        std::tie(start_node, start_found) = graph.get_node(current_start_kmer);

        if (not start_found) {
            BOOST_LOG_TRIVIAL(trace)
                << "Did not find start anchor kmer " << current_start_kmer;
            continue;
        }
        BOOST_LOG_TRIVIAL(trace) << "Found start anchor kmer " << current_start_kmer;

        for (uint_least16_t end_idx = 0; end_idx < end_kmers.size(); end_idx++) {
            const auto& current_end_kmer { end_kmers[end_idx] };

            BOOST_LOG_TRIVIAL(trace)
                << "Searching for end anchor kmer " << current_end_kmer;
            std::tie(end_node, end_found) = graph.get_node(current_end_kmer);

            if (end_found) {
                BOOST_LOG_TRIVIAL(trace)
                    << "Found end anchor kmer " << current_end_kmer;
                DenovoPaths denovo_paths;
                FoundPaths abandoned;
                std::tie(denovo_paths, abandoned) = graph.get_paths_between(
                    start_node, end_node, max_path_length, expected_kmer_covg);

                if (abandoned) {
                    remove_graph_file(GATB_graph_filepath);
                    return;
                }

                if (not denovo_paths.empty()) {
                    // add flank to each sequence - the whole sequence from the start to
                    // the end of the gene
                    for (const std::string &denovo_path : denovo_paths) {
                        const auto start_kmer_offset {
                            candidate_region.max_likelihood_sequence.substr(
                                0, start_idx)
                        };
                        const auto end_kmer_offset {
                            candidate_region.max_likelihood_sequence.substr(
                                candidate_region.max_likelihood_sequence.length()
                                - end_idx)
                        };
                        std::string denovo_path_with_flanks
                            = std::string(candidate_region.left_flanking_sequence)
                                  .append(start_kmer_offset)
                                  .append(denovo_path)
                                  .append(end_kmer_offset)
                                  .append(candidate_region.right_flanking_sequence);

                        bool denovo_path_is_indeed_a_new_sequence = denovo_path_with_flanks != candidate_region.local_node_max_likelihood_sequence;
                        if (denovo_path_is_indeed_a_new_sequence) {
                            candidate_region.denovo_paths.push_back(denovo_path_with_flanks);
                        }
                    }

                    remove_graph_file(GATB_graph_filepath);
                    return;
                }
            } else {
                BOOST_LOG_TRIVIAL(trace)
                    << "Did not find end anchor kmer " << current_end_kmer;
            }
        }
    }
    BOOST_LOG_TRIVIAL(debug)
        << "Could not find any combination of start and end "
           "k-mers with a path between them. Skipping local assembly for "
        << candidate_region.get_name();
    remove_graph_file(GATB_graph_filepath);
}

DenovoDiscovery::DenovoDiscovery(const uint16_t& kmer_size,
    const double& read_error_rate, const int max_nb_paths,
    const uint16_t max_insertion_size,
    const uint16_t min_covg_for_node_in_assembly_graph, const bool clean)
    : kmer_size { kmer_size }
    , read_error_rate { read_error_rate }
    , max_nb_paths { max_nb_paths }
    , max_insertion_size { max_insertion_size }
    , min_covg_for_node_in_assembly_graph { min_covg_for_node_in_assembly_graph }
    , clean_assembly_graph { clean }
{
}

double DenovoDiscovery::calculate_kmer_coverage(
    const uint32_t& read_covg, const uint32_t& ref_length) const
{
    if (ref_length == 0) {
        throw std::invalid_argument("ref_length should be greater than 0.");
    } else if (kmer_size == 0) {
        throw std::invalid_argument("K should be greater than 0.");
    } else if (read_error_rate < 0) {
        throw std::invalid_argument("error_rate should not be a negative value.");
    }

    const auto numerator { read_covg * (ref_length - kmer_size + 1)
        * std::pow(1 - read_error_rate, kmer_size) };
    const auto denominator { ref_length };
    return numerator / denominator;
}
