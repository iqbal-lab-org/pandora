#include "gtest/gtest.h"
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <local_assembly.h>


const std::string fastqPath = "../../test/test_cases/test.fastq";
const int g_kmer_size = 5;

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
    Kmer<>::ModelDirect model{5};
    // we declare an iterator relying on that kmer model
    Kmer<>::ModelDirect::Iterator it{model};

    // we create a data object representation of the sequence (in ASCII format)
    Data data((char *) seq);

    // we configure the iterator with our sequence
    it.setData(data);

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
        Graph graph = Graph::create(Bank::open(fastqPath), "-kmer-size %d -verbose 0 -abundance-min 1", g_kmer_size);

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
        Graph graph = Graph::create(Bank::open(fastqPath), "-kmer-size %d -verbose 0 -abundance-min 1", g_kmer_size);

        // we get an iterator for all nodes of the graph
        GraphIterator<Node> it = graph.iterator();

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
        EXPECT_EQ((int) 2, it.size());

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
                EXPECT_EQ((int) 1, neighbours.size());

                // another way to check the number of successors
                EXPECT_EQ((int) 1, graph.outdegree(current));

                // check the number of predecessors
                EXPECT_EQ((int) 0, graph.indegree(current));

                // we check that it is the correct neighbour
                EXPECT_EQ((std::string) "ATGC", graph.toString(neighbours[0]));
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
            "-kmer-size %d -abundance-min 1 -verbose 0", g_kmer_size
    );
    auto kmer = "AATGT";
    Node real_node;
    bool found;
    std::tie(real_node, found) = get_node(kmer, graph);
    EXPECT_TRUE(found);

    // We get the neighbors of this real node and make sure it has the neighbours we expect
    GraphVector<Node> neighbours = graph.successors(real_node);
    EXPECT_EQ(graph.toString(neighbours[0]), "ATGTC");
}


TEST(GetNodeFromGraph, GivenGraphAndKmer_KmerFoundInGraph) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_kmer_size
    );
    auto kmer = "AATGT";

    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    auto &result = found;
    EXPECT_TRUE(result);
}


TEST(GetNodeFromGraph, GivenGraphAndMissingKmer_KmerNotFoundInGraph) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_kmer_size
    );

    auto kmer = "ACTGT";
    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    auto &result = found;
    EXPECT_FALSE(found);
}


TEST(GetNodeFromGraph, GivenGraphAndKmer_CorrectNodeReturned) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_kmer_size
    );
    auto kmer = "AATGT";

    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    auto result = graph.toString(node);
    auto &expected = kmer;

    EXPECT_EQ(expected, result);
}


TEST(GetNodeFromGraph, GivenGraphAndMissingKmer_CorrectEmptyNodeReturned) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_kmer_size
    );
    auto kmer = "ACTGT";

    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    auto result = node;
    Node expected = {};

    EXPECT_EQ(expected, result);
}


TEST(DFSTest, create) {
    const std::string s1{"AATGTAAGG"};
    const std::string s2{"AATGTCAGG"};
    const std::string s3{"AATGTTAGG"};
    std::vector<std::string> seqs = {s1, s2, s3};

    Graph graph = Graph::create(
            new BankStrings(seqs),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_kmer_size
    );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node("AATGT", graph);

    auto tree = DFS(start_node, graph);

    std::vector<std::string> result;
    const auto end_kmer{"AGG"};
    get_paths_between("AATGT", end_kmer, tree, graph, result);
    std::sort(result.begin(), result.end());

    for (int i = 0; i < result.size(); ++i) {
        EXPECT_EQ(result[i], seqs[i]);
    }
}


TEST(DFSTest, SimpleGraphTwoNodes_ReturnSeqPassedIn) {
    const auto seq{"ATGCAG"};
    const auto start_kmer{"ATGCA"};
    const auto end_kmer{"TGCAG"};

    const Graph graph = Graph::create(
            new BankStrings(seq, NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_kmer_size
            );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = DFS(start_node, graph);

    std::vector<std::string> result;
    get_paths_between(start_kmer, end_kmer, tree, graph, result);

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(result.at(0), seq);
}


TEST(DFSTest, SimpleGraphSixNodes_ReturnSeqPassedIn) {
    const auto seq{"ATGCAGTACA"};
    const auto start_kmer{"ATGCA"};
    const auto end_kmer{"GTACA"};

    const Graph graph = Graph::create(
            new BankStrings(seq, NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_kmer_size
            );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = DFS(start_node, graph);

    std::vector<std::string> result;
    get_paths_between(start_kmer, end_kmer, tree, graph, result);

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(result.at(0), seq);
}
