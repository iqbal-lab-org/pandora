#include <boost/filesystem/path.hpp>
#include "denovo_discovery/local_assembly.h"


bool string_ends_with(std::string const &query, std::string const &ending) {
    if (query.length() < ending.length()) {
        return false;
    }
    return 0 == query.compare(query.length() - ending.length(), ending.length(), ending);
}


std::pair<Node, bool> get_node(const std::string &query_kmer, const Graph &graph) {
    const auto query_node { graph.buildNode(query_kmer.c_str()) };
    Node requested_node;
    bool node_found { false };
    auto nodes_in_graph { graph.iterator() };

    for (nodes_in_graph.first(); !nodes_in_graph.isDone(); nodes_in_graph.next()) {
        const auto &current_node { nodes_in_graph.item() };

        // does not take strand into account
        node_found = (current_node == query_node);
        if (node_found) {
            // now test whether the strands are the same
            // todo make if_strand_same function
            if (graph.toString(current_node) != query_kmer) {
                requested_node = graph.reverse(current_node);
            } else {
                requested_node = current_node;
            }
            break;
        }
    }
    return std::make_pair(requested_node, node_found);
}


// Non-recursive implementation of DFS from "Algorithm Design" - Kleinberg and Tardos (First Edition)
DfsTree depth_first_search_from(const Node &start_node, const Graph &graph) {
    BOOST_LOG_TRIVIAL(debug) << "Starting DFS...";
    std::stack<Node> nodes_to_explore({ start_node });

    std::unordered_set<std::string> explored_nodes;
    DfsTree tree_of_nodes_visited;

    while (not nodes_to_explore.empty()) {
        auto &current_node { nodes_to_explore.top() };
        nodes_to_explore.pop();

        bool previously_explored { explored_nodes.find(graph.toString(current_node)) != explored_nodes.end() };
        if (previously_explored) {
            continue;
        }

        explored_nodes.insert(graph.toString(current_node));
        auto children_of_current_node { graph.successors(current_node) };
        tree_of_nodes_visited[graph.toString(current_node)] = children_of_current_node;

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
get_paths_between(const std::string &start_kmer, const std::string &end_kmer, DfsTree &tree, const Graph &graph,
                  const uint32_t &max_path_length, const double &expected_coverage) {
    BOOST_LOG_TRIVIAL(debug) << "Enumerating all paths in DFS tree between " << start_kmer << " and " << end_kmer;
    std::string path_accumulator { start_kmer.substr(0, start_kmer.length() - 1) };
    DenovoPaths paths_between_queries;
    uint8_t retries { 1 };

    do {
        paths_between_queries.clear();

        const float required_percent_of_expected_covg { retries * COVG_SCALING_FACTOR };
        if (required_percent_of_expected_covg > 1.0) {
            BOOST_LOG_TRIVIAL(debug) << "Abandoning local assembly for slice as too many paths.";
            break;
        }
        build_paths_between(start_kmer, end_kmer, path_accumulator, graph, tree, paths_between_queries, max_path_length,
                            expected_coverage, required_percent_of_expected_covg);
        retries++;
    } while (paths_between_queries.size() > MAX_NUMBER_CANDIDATE_PATHS);

    BOOST_LOG_TRIVIAL(debug) << "Path enumeration complete. There were " << std::to_string(paths_between_queries.size())
                             << " paths found.";
    return paths_between_queries;
}


void build_paths_between(const std::string &start_kmer, const std::string &end_kmer, std::string path_accumulator,
                         const Graph &graph, std::unordered_map<string, GraphVector<Node>> &tree,
                         DenovoPaths &paths_between_queries, const uint32_t &max_path_length,
                         const double &expected_kmer_covg, const float &required_percent_of_expected_covg,
                         uint32_t num_kmers_below_threshold) {
    if (path_accumulator.length() > max_path_length or paths_between_queries.size() > MAX_NUMBER_CANDIDATE_PATHS) {
        return;
    }

    auto start_node { graph.buildNode(start_kmer.c_str()) };
    const auto kmer_coverage { graph.queryAbundance(start_node) };
    const auto max_num_kmers_allowed_below_covg_threshold { start_kmer.length() };

    if (kmer_coverage < (expected_kmer_covg * COVG_SCALING_FACTOR)) {
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
        const auto next_start_kmer { graph.toString(children_of_start_node[i]) };
        build_paths_between(next_start_kmer, end_kmer, path_accumulator, graph, tree, paths_between_queries,
                            max_path_length, expected_kmer_covg, required_percent_of_expected_covg,
                            num_kmers_below_threshold);
    }
}


void remove_graph_file() {
    const fs::path p { "dummy.h5" };
    fs::remove(p);
}


void do_graph_clean(Graph &graph, const uint16_t &num_cores) {
    Simplifications<Graph, Node, Edge> graph_simplifications(graph, num_cores);
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


// todo: the generate kmers functions have a lot of overlapping code - refactor
std::vector<std::string>
generate_start_kmers(const std::string &sequence, const uint16_t &k, uint32_t num_to_generate) {
    const auto seq_len { sequence.length() };
    const auto more_combinations_requested_than_possible { k + (num_to_generate - 1) > seq_len };

    if (k > seq_len) {
        std::vector<std::string> empty;
        return empty;
    } else if (more_combinations_requested_than_possible) {
        num_to_generate = seq_len - k + 1;
    }

    std::vector<std::string> start_kmers;

    for (uint32_t i = 0; i < num_to_generate; i++) {
        start_kmers.emplace_back(sequence.substr(i, k));
    }
    return start_kmers;
}


std::vector<std::string> generate_end_kmers(const std::string &sequence, const uint32_t &k, uint32_t num_to_generate) {
    const auto seq_len { sequence.length() };
    const auto more_combinations_requested_than_possible { k + (num_to_generate - 1) > seq_len };

    if (k > seq_len) {
        std::vector<std::string> empty;
        return empty;
    } else if (more_combinations_requested_than_possible) {
        num_to_generate = seq_len - k + 1;
    }
    std::vector<std::string> end_kmers;

    for (uint32_t i = 0; i < num_to_generate; i++) {
        end_kmers.emplace_back(sequence.substr(seq_len - k - i, k));
    }
    return end_kmers;
}

