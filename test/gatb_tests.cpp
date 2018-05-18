#include "gtest/gtest.h"
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <local_assembly.h>

const std::string fastqPath = "../../test/test_cases/test.fastq";

TEST(GatbFastq, create) {

    // we define a try/catch block in case some method fails
    try {
        // we decalre an input Bank
        IBank *inputBank = Bank::open(fastqPath);

        // create an iterator over this Bank
        Iterator<Sequence> *it = inputBank->iterator();

        // loop over the sequences
        for (it->first(); !it->isDone(); it->next()) {
            // shortcut
            Sequence &seq = it->item();

            // we dump the sequence size and the comment
            std::cout << "[" << seq.getDataSize() << "]" << seq.getComment() << "\n";

            // we dump the sequence
            std::cout << seq.toString() << "\n";

            // we dump the quality score
            std::cout << seq.getQuality() << "\n";
        }
    }
    catch (Exception &e) {
        std::cerr << "EXCEPTION: " << e.getMessage() << "\n";
    }
}


TEST(GatbKmer, create) {
    const char *seq = "CATTGATAGTGGATGGT";

    // we declare a kmer model with kmer sizing 5 letters.
    // note that we want "direct" kmers, not the min(forward,revcomp)
    // default behaviour
    Kmer<>::ModelDirect model {5};
    // we declare an iterator relying on that kmer model
    Kmer<>::ModelDirect::Iterator it {model};

    // we create a data object representation of the sequence (in ASCII format)
    Data data ((char*) seq);

    // we configure the iterator with our sequence
    it.setData (data);

    // and now, we can iterate over the kmers.
    for (it.first(); !it.isDone(); it.next()) {
        std::cout << "kmer " << model.toString(it->value()) << ", value " << it->value() << "\n";
    }
}


TEST(GatbDeBruijnFromFasta, create) {
    try {
        // we load the graph from the given sequence file. in this example, we
        // create the graph only with kmers observed at least 1 times in the data
        // set. we do that using parameter '-abundance-min'.
        Graph graph = Graph::create(Bank::open(fastqPath), "-kmer-size 5 -verbose 0 -abundance-min %d", 1);

        // we dump some information about the graph
        std::cout << graph.getInfo() << "\n";

        // Note: Graph::create will take care of 'bank' object and will delete it if nobody
        // else needs it. i.e no need to call delete on the bank object
    }
    catch (Exception &e) {
        std::cerr << "EXCEPTION: " << e.getMessage() << "\n";
    }
}


TEST(GatbDeBruijnIterateNodes, create) {
    try {
    // we load the graph from the given sequence file. in this example, we
    // create the graph only with kmers observed at least 1 times in the data
    // set. we do that using parameter '-abundance-min'.
    Graph graph = Graph::create(Bank::open(fastqPath), "-kmer-size 5 -verbose 0 -abundance-min %d", 1);

    // we get an iterator for all nodes of the graph
    GraphIterator<Node> it = graph.iterator ();

    // we loop over each node
    for (it.first(); !it.isDone(); it.next()) {
        // the currently iterated node is available with it.item().
        // here, we use it just to dump an ascii representation of each node.
        std::cout << graph.toString(it.item()) << "\n";
    }

    // Note: Graph::create will take care of 'bank' object and will delete it if nobody
    // else needs it. i.e no need to call delete on the bank object
    }
    catch (Exception &e) {
    std::cerr << "EXCEPTION: " << e.getMessage() << "\n";
    }
}


TEST(GatbDeBruijnNodeNeighbours, create) {
    try {
        // We create the graph with a bank holding one sequence and use a
        // specific kmer size and kmer solid abundance set to 1.
        // Using such a string and kmer size we only have two possible nodes
        // in the De Bruijn graph: AATG and ATGC.
        // Of course, in real life we'll load DNA reads from a Fasta/Fastq file.
        Graph graph = Graph::create(
                new BankStrings("AATGC", NULL),
                "-kmer-size 4 -abundance-min 1 -verbose 0"
                );

        // we get an iterator for all nodes of the graph
        GraphIterator<Node> it = graph.iterator();

        // we check that we only have two possible nodes
        EXPECT_EQ((int)2, it.size());

        // we iterate over the nodes
        for (it.first(); !it.isDone(); it.next()) {
            // we get a node
            Node &current = it.item();

            // we get the ascii representation of the current iterated node
            std::string s = graph.toString(current);

            // analyse one of the nodes
            if (s == "AATG") {
                // we get the neighbours of the specific current
                GraphVector<Node> neighbours = graph.successors(current);

                // we check that we only got one successor
                EXPECT_EQ((int)1, neighbours.size());

                // another way to check the number of successors
                EXPECT_EQ((int)1, graph.outdegree(current));

                // check the number of predecessors
                EXPECT_EQ((int)0, graph.indegree(current));

                // we check that it is the correct neighbour
                EXPECT_EQ((std::string)"ATGC", graph.toString(neighbours[0]));
            }
        }

        std::cout << "Graph exploration OK" << "\n";
    }
    catch (Exception &e) {
        std::cerr << "EXCEPTION: " << e.getMessage() << "\n";
    }
}


TEST(GetNodeFromGraph, create) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", NULL),
            "-kmer-size 5 -abundance-min 1 -verbose 0"
    );
    const char *kmer = "AATGTCAGG";
    Node real_node = graph.buildNode(kmer);
    EXPECT_TRUE(get_node(real_node, graph));
    // We get the neighbors of this real node and make sure it has the neighbours we expect
    GraphVector<Node> neighbours = graph.successors(real_node);
    EXPECT_EQ(graph.toString(neighbours[0]), "ATGTC");

    // check the negative case also holds. we create a kmer that should not exist in the graph
    Node fake_node = graph.buildNode("FAKE");
    EXPECT_FALSE(get_node(fake_node, graph));
}

TEST(DFSTest, create) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", "AATGTAAGG", "AATGTATCGTGATG", NULL),
//            new BankStrings("AATC", "AATA", "AATG", NULL),
            "-kmer-size 5 -abundance-min 1 -verbose 0"
    );
    Node start_node = graph.buildNode("AATGT");
    bool node_found = get_node(start_node, graph);
    assert(node_found);

    std::cout << "Sequence 1: AATGTCAGG\n";
    std::cout << "Sequence 2: AATGTAAGG\n";

    std::unordered_map<std::string, GraphVector<Node>>& tree = DFS(start_node, graph);
    std::cout << "DFS Tree:\n";
    for (auto kv : tree) {
        std::cout << kv.first << " ";
        std::cout << kv.second.size() << "\n";
    }
    std::cout << "print_path function output:\n";
    print_path(tree, "AATGT", graph);

}
