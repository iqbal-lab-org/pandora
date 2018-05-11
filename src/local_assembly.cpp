#include <local_assembly.h>

bool kmer_in_graph(const char *kmer, Graph &graph) {
    // get the first node of the kmer
    Node node = graph.buildNode(kmer);

    return graph.contains(node);
}

