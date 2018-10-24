#include <boost/filesystem/path.hpp>
#include "denovo_discovery/local_assembly.h"


bool has_ending(std::string const &fullString, std::string const &ending) {
    if (fullString.length() < ending.length()) {
        return false;
    }
    return 0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending);
}


std::pair<Node, bool> get_node(const std::string &kmer, const Graph &graph) {

    const auto query_node{graph.buildNode(kmer.c_str())};

    Node node;
    bool node_found{false};

    auto it = graph.iterator();
    for (it.first(); !it.isDone(); it.next()) {
        const auto &current_node{it.item()};

        // compare the kmer values for nodes - does not take strand into account
        node_found = (current_node == query_node);
        if (node_found) {
            // now test whether the strands are the same
            if (graph.toString(current_node) != kmer) {
                node = graph.reverse(current_node);
            }
            else {
                node = current_node;
            }
            break;
        }
    }
    return std::make_pair(node, node_found);
}


// Non-recursive implementation of DFS from "Algorithm Design" - Kleinberg and Tardos (First Edition)
DfsTree DFS(const Node &start_node, const Graph &graph) {
    BOOST_LOG_TRIVIAL(debug) << "Starting DFS...";
    std::stack<Node> nodes_to_explore({start_node});

    std::unordered_set<std::string> explored_nodes;
    DfsTree tree;

    while (not nodes_to_explore.empty()) {
        auto &current_node = nodes_to_explore.top();
        nodes_to_explore.pop();

        bool previously_explored = explored_nodes.find(graph.toString(current_node)) != explored_nodes.end();
        if (previously_explored) {
            continue;
        }

        explored_nodes.insert(graph.toString(current_node));

        auto neighbours = graph.successors(current_node);
        tree[graph.toString(current_node)] = neighbours;

        for (unsigned int i = 0; i < neighbours.size(); ++i) {
            auto child = neighbours[i];
            nodes_to_explore.push(child);
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "DFS finished.";
    return tree;
}

/* The aim of this function is to take a DFS tree, and return all paths within this tree that start at start_kmer
 * and end at end_kmer. Allowing for the different combinations in the number of cycles if the path contains any.
 *
 * The associated util function is a recursive function that generates a path down to a "leaf" of the tree and
 * then comes back up to the next unexplored branching point.
 */
Paths get_paths_between(const std::string &start_kmer, const std::string &end_kmer,
                        std::unordered_map<string, GraphVector<Node>> &tree, const Graph &graph,
                        const unsigned long max_path_length, const double &expected_coverage) {
    BOOST_LOG_TRIVIAL(debug) << "Enumerating all paths in DFS tree between " << start_kmer << " and " << end_kmer;
    std::string initial_acc = start_kmer.substr(0, start_kmer.length() - 1);

    Paths result = {};
    get_paths_between_util(start_kmer, end_kmer, initial_acc, graph, tree,
                           result, max_path_length, expected_coverage);
    BOOST_LOG_TRIVIAL(debug) << "Path enumeration complete. There were " << std::to_string(result.size())
                             << " paths found.";
    return result;
}


void get_paths_between_util(const std::string &start_kmer, const std::string &end_kmer, std::string path_accumulator,
                            const Graph &graph, std::unordered_map<string, GraphVector<Node>> &tree, Paths &full_paths,
                            const unsigned long &max_path_length, const double &expected_kmer_covg,
                            unsigned int kmers_below_threshold) {
    // gather information on kmer coverages
    auto start_node{graph.buildNode(start_kmer.c_str())};
    const auto kmer_coverage{graph.queryAbundance(start_node)};

    // do coverage check
    // if there are k k-mers with coverage <= expected_covg * coverage scaling factor - stop recursing for this path
    // todo: revisit this scaling with zam
    if (kmer_coverage < (expected_kmer_covg * COVG_SCALING_FACTOR)) {
        kmers_below_threshold++;
        if (kmers_below_threshold >= start_kmer.length()) {
            return;
        }
    }

    auto &child_nodes = tree[start_kmer];
    auto num_children = child_nodes.size();

    if (path_accumulator.length() > max_path_length) {
        return;
    }

    path_accumulator.push_back(start_kmer.back());

    // makes sure we get all possible cycle repitions up to the maximum length
    if (has_ending(path_accumulator, end_kmer)) {
        full_paths.push_back(path_accumulator);
    }

    for (unsigned int i = 0; i < num_children; ++i) {
        auto kmer = graph.toString(child_nodes[i]);
        get_paths_between_util(kmer, end_kmer, path_accumulator, graph, tree,
                               full_paths, max_path_length, expected_kmer_covg, kmers_below_threshold);
    }
}


void write_paths_to_fasta(const boost::filesystem::path &filepath,
                          const Paths &paths,
                          unsigned long line_width) {
    const auto header = ">path";
    std::ofstream out_file(filepath.string());

    for (const auto &path: paths) {
        out_file << header << "\n";

        for (unsigned long i = 0; i < path.length(); i += line_width) {
            out_file << path.substr(i, line_width) << "\n";
        }
    }

    out_file.close();
    BOOST_LOG_TRIVIAL(debug) << "Local assembly paths written to " << filepath;
}

void local_assembly(const std::vector<std::string> &sequences,
                    std::unordered_set<std::string> &start_kmers,
                    std::unordered_set<std::string> &end_kmers,
                    const fs::path &out_path,
                    const unsigned int kmer_size,
                    const unsigned long max_path_length,
                    const double &expected_coverage,
                    const bool clean_graph,
                    const unsigned int min_coverage) {
    if (sequences.empty()) {
        BOOST_LOG_TRIVIAL(warning) << "Sequences vector to assemble is empty. Skipping local assembly.";
        return;
    }

    Graph graph;  // have to predefine as actually initialisation is inside try block

    // make sure the max_path_length is actually longer than the kmer size
    if (kmer_size > max_path_length) {
        BOOST_LOG_TRIVIAL(warning) << "Kmer size " << std::to_string(kmer_size)
                                   << " is greater than the maximum path length " << std::to_string(max_path_length)
                                   << ". Skipping local assembly.";
        return;
    }

    try {
        graph = Graph::create(
                new BankStrings(sequences),
                "-kmer-size %d -abundance-min %d -verbose 0", kmer_size, min_coverage);
    }
    catch (gatb::core::system::Exception &error) {
        BOOST_LOG_TRIVIAL(warning) << "Couldn't create GATB graph." << "\n\tEXCEPTION: "
                                   << error.getMessage();
        // remove h5 file that GATB has written to file for this graph
        remove_graph_file();
        return;
    }

    if (clean_graph) {
        do_graph_clean(graph);
    }

    Node start_node;
    Node end_node;
    std::string start_kmer, end_kmer;

    bool start_found{false};
    bool end_found{false};
    for (const auto &s_kmer: start_kmers) {
        std::tie(start_node, start_found) = get_node(s_kmer, graph);

        if (not start_found) {
            continue;
        }
        for (const auto &e_kmer: end_kmers) {
            // make sure end kmer doesnt exist in the set of start kmers
            if (start_kmers.find(e_kmer) != start_kmers.end()) {
                continue;
            }

            std::tie(end_node, end_found) = get_node(e_kmer, graph);

            if (end_found) {
                start_kmer = s_kmer;
                end_kmer = e_kmer;
                goto dfs;
            }
        }
    }
    BOOST_LOG_TRIVIAL(warning) << "Could not find any combination of start and end k-mers. Skipping local assembly.";
    remove_graph_file();
    return;

    dfs:

    auto tree = DFS(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, max_path_length, expected_coverage);

    if (not result.empty()) {
        write_paths_to_fasta(out_path, result);
    }

    // remove h5 file that GATB has written to file for this graph
    remove_graph_file();
}


void remove_graph_file() {
    const fs::path p{"dummy.h5"};
    fs::remove(p);
}


void do_graph_clean(Graph &graph, const int num_cores) {
    Simplifications<Graph, Node, Edge> graph_simplifications(graph, num_cores);
    graph_simplifications._doTipRemoval = true;
    graph_simplifications._doBulgeRemoval = false;
    graph_simplifications._doECRemoval = false;

    graph_simplifications._tipLen_Topo_kMult = 2; // remove all tips of length <= k * X bp  [default '2.500000'] set to 0 to turn off
    graph_simplifications._tipLen_RCTC_kMult = 0;  // remove tips that pass coverage criteria, of length <= k * X bp  [default '10.000000'] set to 0 to turn off
    graph_simplifications._tipRCTCcutoff = 2; // tip relative coverage coefficient: mean coverage of neighbors >  X * tip coverage default 2.0
    graph_simplifications.simplify();
}


std::string reverse_complement(const std::string &forward) {
    const auto len = forward.size();
    std::string reverse(len, ' ');
    for (size_t k = 0; k < len; k++) {
        char base = forward[k];
        char magic = base & 2 ? 4 : 21;
        reverse[len - k - 1] = base ^ magic;
    }
    reverse[len] = '\0';
    return reverse;
}


bool file_exists(const std::string &name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

//todo: make start and end kmers a vector or similar. i.e would prefer to keep ordering of kmers
std::unordered_set<std::string>
generate_start_kmers(const std::string &sequence, const unsigned int k, unsigned int n) {
    const auto L{sequence.length()};
    if (k > L) {
        BOOST_LOG_TRIVIAL(error) << "Cannot generate kmers when K " << std::to_string(k)
                                 << " is greater than the length of the sequence " << std::to_string(L);
        std::unordered_set<std::string> empty;
        return empty;
    } else if (k + (n - 1) > L) {  // more combinations are requested than is possible
        n = L - k + 1;  // make n the largest value it can take on
    }

    std::unordered_set<std::string> kmers;

    for (unsigned int i = 0; i < n; i++) {
        const auto kmer{sequence.substr(i, k)};
        kmers.insert(kmer);
    }
    return kmers;
}

std::unordered_set<std::string> generate_end_kmers(const std::string &sequence, const unsigned int k, unsigned int n) {
    const auto L{sequence.length()};
    if (k > L) {
        BOOST_LOG_TRIVIAL(error) << "Cannot generate kmers when K " << std::to_string(k)
                                 << " is greater than the length of the sequence " << std::to_string(L);
        std::unordered_set<std::string> empty;
        return empty;
    } else if (k + (n - 1) > L) {  // more combinations are requested than is possible
        n = L - k + 1;  // make n the largest value it can take on
    }
    std::unordered_set<std::string> kmers;

    for (unsigned int i = 0; i < n; i++) {
        const auto kmer{sequence.substr(L - k - i, k)};
        kmers.insert(kmer);
    }
    return kmers;
}