#include <local_assembly.h>


bool get_node(Node &node, Graph &graph, Node &found_node) {
    GraphIterator<Node> it = graph.iterator ();
    std::string required_kmer = graph.toString(node);
    for (it.first(); !it.isDone(); it.next())
    {
        Node& current = it.item();
        if (graph.toString(current) == required_kmer) {
            found_node = current;
            return true;
        }
    }
    return false;
}

u_int64_t graph_size(Graph &graph) {
    // We get an iterator for all nodes of the graph.
    GraphIterator<Node> it = graph.iterator ();

    // returns the number of nodes in the graph
    return it.size();
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
std::unordered_map<std::string, std::string>& DFS(Node &start_node, Graph &graph) {
    std::stack<Node> nodes_to_explore;  // S from pseudocode
    nodes_to_explore.push(start_node);  // s from pseudocode
    std::unordered_map<std::string, std::string> parent;
    DfsTree tree { DfsTree() };
    std::set<std::string> explored;
    bool u_explored;

    while (!(nodes_to_explore.empty())) {
        // Take a node u from S
        Node &current_node { nodes_to_explore.top() };  // u from pseudocode
        nodes_to_explore.pop();
        // If Explored[u] = false then
        u_explored = explored.find(graph.toString(current_node)) != explored.end();
        if (!u_explored) {
            // Set Explored[u] = true
            explored.insert(graph.toString(current_node));
            // If u != s then
            if (current_node != start_node) {
                // Add edge (u, parent[u]) to tree T
                // todo
            }
            // For each edge (u,v) incident to u
            // We get the neighbors of this current node
            GraphVector<Node> neighbors = graph.successors(current_node);

            // We loop each node.
            for (int count = 0; count < neighbors.size(); ++count)
            {
                Node v = neighbors[count];
                // Add v to the stack S
                nodes_to_explore.push(v);
                // Set parent[v] = u
                parent[graph.toString(v)] = graph.toString(current_node);
            }
        }
    }
    return parent;
}


DfsNode::DfsNode() {
    DfsNode *parent{nullptr};
    DfsNode *left_child{nullptr};
    DfsNode *right_sibling{nullptr};
}
DfsNode *const DfsNode::get_parent() {
        return parent;
}
DfsNode *const DfsNode::get_left_child() {
        return left_child;
}
DfsNode *const DfsNode::get_right_sibling() {
        return right_sibling;
}
void DfsNode::set_parent(DfsNode *node) {
    parent = node;
}
void DfsNode::set_left_child(DfsNode *node) {
    left_child = node;
}
void DfsNode::set_right_sibling(DfsNode *node) {
    right_sibling = node;
}
bool DfsNode::has_children() {
    return left_child != nullptr;
}
bool DfsNode::is_rightmost_child() {
    return right_sibling == nullptr;
}
bool DfsNode::is_root() {
    return parent == nullptr;
}



DfsTree::DfsTree() {
    DfsNode *root { nullptr };
}
