#include "denovo_discovery/local_assembly.h"


LocalAssemblyGraph& LocalAssemblyGraph::operator=(const Graph &graph) {
    if (this != &graph) {
        _kmerSize = graph._kmerSize;
        _storageMode = graph._storageMode;
        _name = graph._name;
        _info = graph._info;
        _bloomKind = graph._bloomKind;
        _debloomKind = graph._debloomKind;
        _debloomImpl = graph._debloomImpl;
        _branchingKind = graph._branchingKind;
        _state = graph._state;

        setStorage(graph._storage);

        if (graph._variant) { *((GraphDataVariant *) _variant) = *((GraphDataVariant *) graph._variant); }
    }
    return *this;
}


bool string_ends_with(std::string const &query, std::string const &ending) {
    if (query.length() < ending.length()) {
        return false;
    }
    return 0 == query.compare(query.length() - ending.length(), ending.length(), ending);
}


std::pair<Node, bool> LocalAssemblyGraph::get_node(const std::string &query_kmer) {
    const auto query_node { buildNode(query_kmer.c_str()) };
    Node requested_node;
    bool node_found { false };
    auto nodes_in_graph { iterator() };

    for (nodes_in_graph.first(); !nodes_in_graph.isDone(); nodes_in_graph.next()) {
        const auto &current_node { nodes_in_graph.item() };
        node_found = (current_node == query_node);

        if (node_found) {
            const bool strands_are_opposite { toString(current_node) != query_kmer };

            if (strands_are_opposite) {
                requested_node = reverse(current_node);
            } else {
                requested_node = current_node;
            }
            break;
        }
    }
    return std::make_pair(requested_node, node_found);
}


// Non-recursive implementation of DFS from "Algorithm Design" - Kleinberg and Tardos (First Edition)
DfsTree LocalAssemblyGraph::depth_first_search_from(const Node &start_node, bool reverse) {
    BOOST_LOG_TRIVIAL(debug) << "Starting DFS...";
    std::stack <Node> nodes_to_explore({ start_node });

    std::unordered_set <std::string> explored_nodes;
    DfsTree tree_of_nodes_visited;

    while (not nodes_to_explore.empty()) {
        auto &current_node { nodes_to_explore.top() };
        nodes_to_explore.pop();

        bool previously_explored { explored_nodes.find(toString(current_node)) != explored_nodes.end() };
        if (previously_explored) {
            continue;
        }

        explored_nodes.insert(toString(current_node));
        auto children_of_current_node { (reverse ? predecessors(current_node) : successors(current_node)) };
        tree_of_nodes_visited[toString(current_node)] = children_of_current_node;

        for (unsigned int i = 0; i < children_of_current_node.size(); ++i) {
            nodes_to_explore.push(children_of_current_node[i]);
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "DFS finished.";
    return tree_of_nodes_visited;
}


/* The aim of this function is to take a DFS tree, and return all paths within this tree that start at start_kmer
 * and end at end_kmer. Allowing for the different combinations in the number of cycles if the path contains any.
 *
 * The associated util function is a recursive function that generates a path down to a "leaf" of the tree and
 * then comes back up to the next unexplored branching point.
 */
DenovoPaths
LocalAssemblyGraph::get_paths_between(const Node &start_node, const Node &end_node,
                                      const uint32_t &max_path_length, const double &expected_coverage) {
    DenovoPaths paths_between_queries;
    const std::string start_kmer = toString(start_node);
    const std::string end_kmer = toString(end_node);

    auto tree { depth_first_search_from(start_node) };

    //check if end node is in forward tree, if not just return
    bool is_end_kmer_reachable_from_start_kmer =
            tree.find(end_kmer) == tree.end();
    if (is_end_kmer_reachable_from_start_kmer) {
        return paths_between_queries;
    }

    auto reverse_tree { depth_first_search_from(end_node, true) };


    BOOST_LOG_TRIVIAL(debug) << "Enumerating all paths in DFS tree between " << start_kmer << " and " << end_kmer;
    std::string path_accumulator { start_kmer.substr(0, start_kmer.length() - 1) };
    uint8_t retries { 1 };

    do {
        paths_between_queries.clear();

        const float required_percent_of_expected_covg { retries * COVG_SCALING_FACTOR };
        if (required_percent_of_expected_covg > 1.0) {
            BOOST_LOG_TRIVIAL(debug) << "Abandoning local assembly for slice as too many paths.";
            break;
        }
        build_paths_between(start_kmer, end_kmer, path_accumulator, tree, reverse_tree, paths_between_queries, max_path_length,
                            expected_coverage, required_percent_of_expected_covg);
        retries++;
    } while (paths_between_queries.size() > MAX_NUMBER_CANDIDATE_PATHS);

    BOOST_LOG_TRIVIAL(debug) << "Path enumeration complete. There were " << std::to_string(paths_between_queries.size())
                             << " paths found.";
    return paths_between_queries;
}


void LocalAssemblyGraph::build_paths_between(const std::string &start_kmer, const std::string &end_kmer,
                                             std::string path_accumulator,
                                             DfsTree &tree,
                                             DfsTree &reverse_tree,
                                             DenovoPaths &paths_between_queries, const uint32_t &max_path_length,
                                             const double &expected_kmer_covg,
                                             const float &required_percent_of_expected_covg,
                                             uint32_t num_kmers_below_threshold) {
    if (path_accumulator.length() > max_path_length or paths_between_queries.size() > MAX_NUMBER_CANDIDATE_PATHS) {
        return;
    }

    bool start_kmer_cannot_reach_end_kmer = reverse_tree.find(start_kmer) == reverse_tree.end();
    if (start_kmer_cannot_reach_end_kmer) {
        return;
    }


    auto start_node { buildNode(start_kmer.c_str()) };
    const auto kmer_coverage { queryAbundance(start_node) };
    const auto max_num_kmers_allowed_below_covg_threshold { start_kmer.length() };

    if (kmer_coverage < (expected_kmer_covg * required_percent_of_expected_covg)) {
        num_kmers_below_threshold++;
        if (num_kmers_below_threshold >= max_num_kmers_allowed_below_covg_threshold) {
            return;
        }
    }

    path_accumulator.push_back(start_kmer.back());

    if (string_ends_with(path_accumulator, end_kmer) and path_accumulator.length() > end_kmer.length()) {
        paths_between_queries.push_back(path_accumulator);
        // we dont return here so as to make sure we get all possible cycle repetitions up to the maximum length
    }

    auto &children_of_start_node { tree[start_kmer] };
    const auto num_children { children_of_start_node.size() };

    for (unsigned int i = 0; i < num_children; ++i) {
        const auto next_start_kmer { toString(children_of_start_node[i]) };
        build_paths_between(next_start_kmer, end_kmer, path_accumulator, tree, reverse_tree, paths_between_queries, max_path_length,
                            expected_kmer_covg, required_percent_of_expected_covg, num_kmers_below_threshold);
    }
}


void remove_graph_file() {
    const fs::path p { "dummy.h5" };
    fs::remove(p);
}


void clean(Graph &graph, const uint16_t &num_cores) {
    Simplifications <Graph, Node, Edge> graph_simplifications(graph, num_cores);
    graph_simplifications._doTipRemoval = true;
    graph_simplifications._doBulgeRemoval = false;
    graph_simplifications._doECRemoval = false;

    graph_simplifications
            ._tipLen_Topo_kMult = 2; // remove all tips of length <= k * X bp  [default '2.500000'] set to 0 to turn off
    graph_simplifications
            ._tipLen_RCTC_kMult = 0;  // remove tips that pass coverage criteria, of length <= k * X bp  [default '10.000000'] set to 0 to turn off
    graph_simplifications
            ._tipRCTCcutoff = 2; // tip relative coverage coefficient: mean coverage of neighbors >  X * tip coverage default 2.0
    graph_simplifications.simplify();
}


std::string reverse_complement(const std::string &forward) {
    const auto len { forward.size() };
    std::string reverse(len, ' ');
    for (size_t k = 0; k < len; k++) {
        const char base { forward[k] };
        const char magic = base & 2 ? 4 : 21;
        reverse[len - k - 1] = base ^ magic;
    }
    reverse[len] = '\0';
    return reverse;
}


std::vector <std::string>
generate_start_kmers(const std::string &sequence, const uint16_t &k, uint32_t num_to_generate) {
    const auto seq_len { sequence.length() };
    const auto more_combinations_requested_than_possible { k + (num_to_generate - 1) > seq_len };

    if (more_combinations_requested_than_possible) {
        num_to_generate = seq_len - k + 1;
    }

    const auto generate_kmers_from { sequence.substr(0, num_to_generate + (k - 1)) };

    return all_kmers_in(generate_kmers_from, k);
}


std::vector <std::string> generate_end_kmers(const std::string &sequence, const uint32_t &k, uint32_t num_to_generate) {
    const auto seq_len { sequence.length() };
    const auto more_combinations_requested_than_possible { k + (num_to_generate - 1) > seq_len };

    if (more_combinations_requested_than_possible) {
        num_to_generate = seq_len - k + 1;
    }

    const auto generate_kmers_from { sequence.substr(seq_len - (num_to_generate + (k - 1))) };
    auto end_kmers { all_kmers_in(generate_kmers_from, k) };

    std::reverse(end_kmers.begin(), end_kmers.end());

    return end_kmers;
}


std::vector <std::string> all_kmers_in(const std::string &query, const uint_least8_t k_size) {
    std::vector <std::string> kmers;

    if (k_size > query.length()) {
        return kmers;
    }

    for (size_t i = 0; i < query.length() - (k_size - 1); i++) {
        kmers.emplace_back(query.substr(i, k_size));
    }

    return kmers;
}