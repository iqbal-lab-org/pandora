#include "local_assembly.h"


std::pair<Node, bool> get_node(const std::string &kmer, Graph &graph) {
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

        auto neighbors = graph.successors(current_node);
        tree[graph.toString(current_node)] = neighbors;

        for (auto i = 0; i < neighbors.size(); ++i) {
            Node child = neighbors[i];
            nodes_to_explore.push(child);
        }
    }
    return tree;
}

void print_path(DfsTree &tree,
                const std::string &start_node,
                Graph &graph,
                std::vector<std::string> &result) {
    std::string initial_acc = start_node.substr(0, start_node.size() - 1);
    helper(start_node, initial_acc, graph, tree, result);
}

void helper(const std::string &node, std::string acc, Graph &graph, DfsTree &tree,
            std::vector<std::string> &result) {
    size_t num_children = tree[node].size();
    if (num_children == 0) {
        result.push_back(acc + node.back());
    } else {
        acc += node.back();
        for (int i = 0; i < num_children; ++i) {
            helper(graph.toString(tree[node][i]), acc, graph, tree, result);
        }
    }
}
