#include "local_assembly.h"


bool has_ending(std::string const &fullString, std::string const &ending) {
    if (fullString.length() < ending.length())
        return false;
    return 0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending);
}


std::pair<Node, bool> get_node(const std::string &kmer, const Graph &graph) {
    Node node = {};
    bool found = false;

    auto it = graph.iterator();
    for (it.first(); !it.isDone(); it.next()) {
        const auto &current = it.item();

        bool node_found = graph.toString(current) == kmer;
        if (node_found) {
            node = current;
            found = true;
            break;
        }
    }
    return std::make_pair(node, found);
}


/* Non-recursive implementation of DFS from "Algorithm Design" - Kleinberg and Tardos (First Edition)
 *
 * DFS(s):
 *     Initialise S to be a stack with one element s
 *     Initialise parent to be an array
 *     Initialise dfs tree T
 *     While S is not empty
 *         Take a node u from S
 *         If Explored[u] = false then
 *             Set Explored[u] = true
 *             If u != s then
 *                 Add edge (u, parent[u]) to tree T
 *             Endif
 *             For each edge (u,v) incident to u
 *                 Add v to the stack S
 *                 Set parent[v] = u
 *             Endfor
 *         Endif
 *     Endwhile
 *     Return T
 */
DfsTree DFS(const Node &start_node, const Graph &graph) {
    std::stack<Node> nodes_to_explore({start_node});

    std::set<std::string> explored_nodes;
    DfsTree tree = {};

    while (not nodes_to_explore.empty()) {
        auto current_node = nodes_to_explore.top();
        nodes_to_explore.pop();

        bool previously_explored = explored_nodes.find(graph.toString(current_node)) != explored_nodes.end();
        if (previously_explored)
            continue;

        explored_nodes.insert(graph.toString(current_node));

        auto neighbours = graph.successors(current_node);
        tree[graph.toString(current_node)] = neighbours;

        for (auto i = 0; i < neighbours.size(); ++i) {
            auto child = neighbours[i];
            nodes_to_explore.push(child);
        }
    }
    return tree;
}

/* The aim of this function is to take a DFS tree, and return all paths within this tree that start at start_kmer
 * and end at end_kmer. Allowing for the different combinations in the number of cycles if the path contains any.
 *
 * The associated util function is a recursive function that generates a path down to a "leaf" of the tree and
 * then comes back up to the next unexplored branching point.
 */
Paths get_paths_between(const std::string &start_kmer,
                        const std::string &end_kmer,
                        DfsTree &tree,
                        const Graph &graph) {
    std::string initial_acc = start_kmer.substr(0, start_kmer.length() - 1);

    Paths result = {};
    get_paths_between_util(start_kmer, end_kmer, initial_acc, graph, tree, result);
    return result;
}


void get_paths_between_util(const std::string &start_kmer,
                            const std::string &end_kmer,
                            std::string path_accumulator,
                            const Graph &graph,
                            DfsTree &tree,
                            Paths &full_paths,
                            const long max_length) {
    auto &child_nodes = tree[start_kmer];
    auto num_children = child_nodes.size();

    if (path_accumulator.length() > max_length)
        return;

    path_accumulator.push_back(start_kmer.back());

    // makes sure we get all possible cycle repitions up to the maximum length
    if (has_ending(path_accumulator, end_kmer)) {
        full_paths.push_back(path_accumulator);
    }

    for (int i = 0; i < num_children; ++i) {
        auto kmer = graph.toString(child_nodes[i]);
        get_paths_between_util(kmer,
                               end_kmer,
                               path_accumulator,
                               graph,
                               tree,
                               full_paths);
    }
}


void write_paths_to_fasta(const std::string &filepath, Paths &paths, unsigned long line_width) {
    const auto header = ">path";
    std::ofstream out_file(filepath);

    for (auto &path: paths) {
        out_file << header << "\n";

        for (unsigned long i = 0; i < path.length(); i += line_width) {
            out_file << path.substr(i, line_width) << "\n";
        }
    }

    out_file.close();
}


void local_assembly(const std::string &filepath,
                    std::string &start_kmer,
                    std::string &end_kmer,
                    const std::string &out_path,
                    const int kmer_size) {

    Graph graph;  // have to predefine as actually initialisation is inside try block

    // check if filepath exists
    const bool exists {file_exists(filepath)};
    if (not exists) {
        std::cerr << filepath << " does not exist. Skipping local assembly.\n";
        return;
    }

    try {
        const Graph graph = Graph::create(
                Bank::open(filepath),
                "-kmer-size %d -abundance-min 1 -verbose 0", kmer_size
        );
    }
    catch (gatb::core::system::Exception &error){
        std::cerr << "Couldn't create GATB graph for " << filepath << "\n";
        std::cerr << error.getMessage() << "\n";
        return;
    }

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);
    if (not found) {
        auto tmp_copy = start_kmer;
        start_kmer = reverse_complement(end_kmer);
        end_kmer = reverse_complement(tmp_copy);
        std::tie(start_node, found) = get_node(start_kmer, graph);
        if (not found) {
            std::cerr << "Start kmer not found in " << filepath << "\n";
            return;
        }
    }

    auto tree = DFS(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph);
    write_paths_to_fasta(out_path, result);
}


std::string reverse_complement(const std::string forward) {
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


bool file_exists(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}