#include "local_assembly.h"


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
            Node child = neighbours[i];
            nodes_to_explore.push(child);
        }
    }

//    for (auto &kv: tree) {
//        std::cout << "Key: " << kv.first << "\t Successors: " << kv.second.size() << "\n";
//    }

    return tree;
}

void get_paths_between(const std::string &start_kmer,
                       const std::string &end_kmer,
                       DfsTree &tree,
                       const Graph &graph,
                       Paths &result) {
    std::string initial_acc = start_kmer.substr(0, start_kmer.length() - 1);
    Paths full_paths;

    get_paths_between_util(start_kmer, end_kmer, initial_acc, graph, tree,
                           full_paths);

    for (auto &path : full_paths) {
        // find last occurrence of end kmer in current path
        size_t found = path.rfind(end_kmer);

        if (found == std::string::npos)  // if it wasnt found, skip
            continue;

        const std::string trimmed_path = path.substr(0, found + end_kmer.length());
        result.insert(trimmed_path);
    }
}

void get_paths_between_util(const std::string &start_kmer, const std::string &end_kmer, std::string acc,
                            const Graph &graph,
                            DfsTree &tree, Paths &full_paths) {
    size_t num_children = tree[start_kmer].size();

    if (num_children == 0 || acc.length() > g_max_length) {
        full_paths.insert(acc + start_kmer.back());
    }
    else {
        acc += start_kmer.back();

        // makes sure we get all possible cycle repitions up to the maximum length
        if (has_ending(acc, end_kmer)) {
            full_paths.insert(acc);
        }

        for (int i = 0; i < num_children; ++i) {
            get_paths_between_util(graph.toString(tree[start_kmer][i]), end_kmer, acc, graph, tree,
                                   full_paths);

        }
    }
}


bool has_ending(std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}