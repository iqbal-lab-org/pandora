#include "gtest/gtest.h"
#include <gatb/gatb_core.hpp>
#include <iostream>
#include "denovo_discovery/local_assembly.h"
#include <cstdio>
#include <unordered_set>
#include <algorithm>


const uint32_t g_test_kmer_size = 5;
const uint32_t g_test_max_path = 50;


TEST(GetNodeFromGraph, LowestKmerOfNode_KmerFoundInGraphAndNeighbour) {
    Graph graph = Graph::create(new BankStrings("AATGTCAGG", NULL), "-kmer-size %d -abundance-min 1 -verbose 0",
                                g_test_kmer_size);
    auto kmer = "AATGT";
    Node real_node;
    bool found;
    std::tie(real_node, found) = get_node(kmer, graph);
    EXPECT_TRUE(found);

    // We get the neighbors of this real node and make sure it has the neighbours we expect
    GraphVector<Node> neighbours = graph.successors(real_node);
    EXPECT_EQ(graph.toString(neighbours[0]), "ATGTC");

    auto result = graph.toString(real_node);
    auto &expected = kmer;

    EXPECT_EQ(expected, result);
    remove_graph_file();

}


TEST(GetNodeFromGraph, HighestKmerOfNode_KmerFoundInGraph) {
    Graph graph = Graph::create(new BankStrings("AATGTCAGG", NULL), "-kmer-size %d -abundance-min 1 -verbose 0",
                                g_test_kmer_size);
    auto kmer = "ACATT";

    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    auto result = graph.toString(node);
    auto &expected = kmer;

    EXPECT_EQ(expected, result);
    remove_graph_file();
}


TEST(GetNodeFromGraph, HighestKmerOfNode_KmerFoundInGraphAndNeighbour) {
    Graph graph = Graph::create(new BankStrings("TCGTTGTCACT", NULL), "-kmer-size %d -abundance-min 1 -verbose 0",
                                g_test_kmer_size);

    const auto kmer { "TCGTT" };

    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    EXPECT_TRUE(found);

    const auto result = graph.toString(node);
    const auto &expected = kmer;

    EXPECT_EQ(expected, result);

    const auto expected_neighbour { "CGTTG" };

    GraphVector<Node> neighbours = graph.successors(node);
    const auto result_neighbour { graph.toString(neighbours[0]) };

    EXPECT_EQ(result_neighbour, expected_neighbour);

    remove_graph_file();
}


TEST(GetNodeFromGraph, LowestKmerOfStartNode_KmerFoundInGraphButNotNeighbourOfStart) {
    Graph graph = Graph::create(new BankStrings("TCGTTGTCACT", NULL), "-kmer-size %d -abundance-min 1 -verbose 0",
                                g_test_kmer_size);

    const auto kmer { "AACGA" };

    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    EXPECT_TRUE(found);

    const auto result = graph.toString(node);
    const auto &expected = kmer;

    EXPECT_EQ(expected, result);

    const auto expected_neighbour { "CGTTG" };

    GraphVector<Node> neighbours = graph.successors(node);
    const auto result_neighbour { graph.toString(neighbours[0]) };

    EXPECT_NE(result_neighbour, expected_neighbour);  // NE = not equal

    remove_graph_file();
}


TEST(GetNodeFromGraph, RevcompKmerOfInitialSeq_KmerFoundInGraphAndCorrectNeighbourFound) {
    Graph graph = Graph::create(new BankStrings("TCGTTGTCACT", NULL), "-kmer-size %d -abundance-min 1 -verbose 0",
                                g_test_kmer_size);

    const auto kmer { "TGACA" };

    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    EXPECT_TRUE(found);

    const auto result = graph.toString(node);
    const auto &expected = kmer;

    EXPECT_EQ(expected, result);

    const auto expected_neighbour { "GACAA" };

    GraphVector<Node> neighbours = graph.successors(node);
    const auto result_neighbour { graph.toString(neighbours[0]) };

    EXPECT_EQ(result_neighbour, expected_neighbour);

    remove_graph_file();
}


TEST(GetNodeFromGraph, NonExistentKmer_NotFoundInGraphAndNodeEmpty) {
    Graph graph = Graph::create(new BankStrings("AATGTCAGG", NULL), "-kmer-size %d -abundance-min 1 -verbose 0",
                                g_test_kmer_size);
    auto kmer = "ACTGT";

    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    auto result = node;
    Node expected = {};

    EXPECT_EQ(expected, result);
    remove_graph_file();
}


TEST(GetPathsBetweenTest, OnlyReturnPathsBetweenStartAndEndKmers) {
    const std::string s1 { "AATGTAAGG" };
    const std::string s2 { "AATGTCAGG" };
    const std::string s3 { "AATGTTAGG" };
    std::vector<std::string> seqs = { s1, s2, s3 };

    Graph graph = Graph::create(new BankStrings(seqs), "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size);

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node("AATGT", graph);

    auto tree = depth_first_search_from(start_node, graph);

    const auto end_kmer { "AGG" };
    auto result = get_paths_between("AATGT", end_kmer, tree, graph, g_test_max_path);

    DenovoPaths expected_seqs(seqs.begin(), seqs.end());
    EXPECT_EQ(result, expected_seqs);
    remove_graph_file();
}


TEST(DepthFirstSearchFromTest, SimpleGraphTwoNodesReturnSeqPassedIn) {
    const auto seq { "ATGCAG" };
    const auto start_kmer { "ATGCA" };
    const auto end_kmer { "TGCAG" };

    const Graph graph = Graph::create(new BankStrings(seq, NULL), "-kmer-size %d -abundance-min 1 -verbose 0",
                                      g_test_kmer_size);

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = depth_first_search_from(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_test_max_path);

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(*result.begin(), seq);
    remove_graph_file();
}


TEST(DepthFirstSearchFromTest, SimpleGraphSixNodesReturnSeqPassedIn) {
    const auto seq { "ATGCAGTACA" };
    const auto start_kmer { "ATGCA" };
    const auto end_kmer { "GTACA" };

    const Graph graph = Graph::create(new BankStrings(seq, NULL), "-kmer-size %d -abundance-min 1 -verbose 0",
                                      g_test_kmer_size);

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = depth_first_search_from(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_test_max_path);

    bool original_seq_found = false;
    // make sure all paths begin and end with correct kmer
    for (auto &path: result) {
        EXPECT_EQ(path.substr(0, g_test_kmer_size), start_kmer);
        EXPECT_EQ(path.substr(path.length() - g_test_kmer_size, path.length()), end_kmer);

        if (path == seq) {
            original_seq_found = true;
        }
    }

    EXPECT_TRUE(original_seq_found);
    remove_graph_file();
}


TEST(DepthFirstSearchFromTest, TwoReadsSameSequenceReturnOneSequence) {
    const auto seq1 { "ATGCAG" };
    const auto seq2 { "ATGCAG" };
    std::vector<std::string> seqs = { seq1, seq2 };
    const auto start_kmer { "ATGCA" };
    const auto end_kmer { "TGCAG" };

    const Graph graph = Graph::create(new BankStrings(seqs), "-kmer-size %d -abundance-min 1 -verbose 0",
                                      g_test_kmer_size);

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = depth_first_search_from(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_test_max_path);

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(*result.begin(), seq1);
    remove_graph_file();
}


TEST(DepthFirstSearchFromTest, TwoReadsOneVariantReturnOriginalTwoSequences) {
    const std::string seq1 { "ATGCAGTACAA" };
    const std::string seq2 { "ATGCATTACAA" };
    std::vector<std::string> seqs = { seq1, seq2 };
    const auto start_kmer { "ATGCA" };
    const auto end_kmer { "TACAA" };

    const Graph graph = Graph::create(new BankStrings(seqs), "-kmer-size %d -abundance-min 1 -verbose 0",
                                      g_test_kmer_size);

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = depth_first_search_from(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_test_max_path);

    int original_seq_found = 0;
    for (auto &path: result) {
        EXPECT_EQ(path.substr(0, g_test_kmer_size), start_kmer);
        EXPECT_EQ(path.substr(path.length() - g_test_kmer_size, path.length()), end_kmer);

        if (path.length() == seq1.length()) {
            EXPECT_TRUE(std::find(seqs.begin(), seqs.end(), path) != seqs.end());
            ++original_seq_found;
        }
    }
    EXPECT_EQ(original_seq_found, seqs.size());
    remove_graph_file();
}


TEST(DepthFirstSearchFromTest, ThreeReadsTwoVariantsReturnOriginalSequences) {
    const std::string seq1 { "ATGCAGTACAA" };
    const std::string seq2 { "ATGCATTACAA" };
    const std::string seq3 { "ATGCACTACAA" };
    std::vector<std::string> seqs = { seq1, seq2, seq3 };
    const auto start_kmer { "ATGCA" };
    const auto end_kmer { "TACAA" };

    const Graph graph = Graph::create(new BankStrings(seqs), "-kmer-size %d -abundance-min 1 -verbose 0",
                                      g_test_kmer_size);

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = depth_first_search_from(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_test_max_path);

    int original_seq_found = 0;
    for (auto &path: result) {
        EXPECT_EQ(path.substr(0, g_test_kmer_size), start_kmer);
        EXPECT_EQ(path.substr(path.length() - g_test_kmer_size, path.length()), end_kmer);

        if (path.length() == seq1.length()) {
            bool path_in_expected = std::find(seqs.begin(), seqs.end(), path) != seqs.end();
            if (path_in_expected) {
                ++original_seq_found;
            }
        }
    }
    EXPECT_EQ(original_seq_found, seqs.size());
    remove_graph_file();
}


TEST(DepthFirstSearchFromTest, TwoReadsTwoVariantsReturnOriginalTwoSequencesPlusTwoMosaics) {
    const std::string seq1 { "TTGGTCATCCCATTATG" };
    const std::string seq2 { "TTGGTGATCCCGTTATG" };
    std::vector<std::string> seqs = { seq1, seq2 };
    const auto start_kmer { "TTGGT" };
    const auto end_kmer { "TTATG" };

    const Graph graph = Graph::create(new BankStrings(seqs), "-kmer-size %d -abundance-min 1 -verbose 0",
                                      g_test_kmer_size);

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = depth_first_search_from(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_test_max_path);

    // add other expected paths due to variants
    std::vector<std::string> expected_seqs = { seq1, seq2, "TTGGTGATCCCATTATG", "TTGGTCATCCCGTTATG" };

    std::sort(result.begin(), result.end());
    std::sort(expected_seqs.begin(), expected_seqs.end());

    EXPECT_EQ(result, expected_seqs);
    remove_graph_file();

}


TEST(DepthFirstSearchFromTest, ThreeReadsOneReverseComplimentReturnPathsForStrandOfStartAndEndKmers) {
    const std::string seq1 { "ATGTG" };
    const std::string seq2 { "TGTGC" };
    const std::string seq3 { "TGCAC" };
    std::vector<std::string> seqs = { seq1, seq2, seq3 };
    const auto start_kmer { "ATGTG" };
    const auto end_kmer { "GTGCA" };

    const Graph graph = Graph::create(new BankStrings(seqs), "-kmer-size %d -abundance-min 1 -verbose 0",
                                      g_test_kmer_size);

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = depth_first_search_from(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_test_max_path);

    // add other expected paths due to variants
    const std::string expected_seq = "ATGTGCA";

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(*result.begin(), expected_seq);
    remove_graph_file();

}


TEST(DepthFirstSearchFromTest, SimpleCycleReturnPathsOfLengthsUpToMaxPathLengthCycling) {
    const std::string seq1 { "ATATATATA" };
    const std::string seq2 { "TATAT" };
    std::vector<std::string> seqs = { seq1, seq2 };
    const auto start_kmer { "ATATA" };
    const auto end_kmer { "TATAT" };

    const Graph graph = Graph::create(new BankStrings(seqs), "-kmer-size %d -abundance-min 1 -verbose 0",
                                      g_test_kmer_size);

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = depth_first_search_from(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_test_max_path);

    const std::string min_expected_seq = "ATATAT";
    bool is_in = false;

    for (auto &path: result) {
        if (path == min_expected_seq) {
            is_in = true;
        }
    }

    EXPECT_TRUE(is_in);
    remove_graph_file();
}


TEST(StringEndsWithTest, endsWithReturnTrue) {
    std::string test = "binary";
    std::string ending = "nary";

    EXPECT_TRUE(string_ends_with(test, ending));
}


TEST(StringEndsWithTest, doesNotHaveEndingReturnFalse) {
    std::string test = "tertiary";
    std::string ending = "nary";

    EXPECT_FALSE(string_ends_with(test, ending));
}


TEST(StringEndsWithTest, endingLongerThanQueryReturnFalse) {
    std::string test = "ry";
    std::string ending = "nary";

    EXPECT_FALSE(string_ends_with(test, ending));
}


TEST(ReverseComplementTest, SingleBaseReturnCompliment) {
    const auto seq = "A";
    const auto expected = "T";
    auto result = reverse_complement(seq);
    EXPECT_EQ(expected, result);
}


TEST(ReverseComplementTest, TwoBasesReturnCompliment) {
    const auto seq = "AA";
    const auto expected = "TT";
    auto result = reverse_complement(seq);
    EXPECT_EQ(expected, result);
}


TEST(ReverseComplementTest, AllBasesReturnCompliment) {
    const auto seq = "ACTGCA";
    const auto expected = "TGCAGT";
    auto result = reverse_complement(seq);
    EXPECT_EQ(expected, result);
}


TEST(ReverseComplementTest, PalindromeReturnCompliment) {
    const auto seq = "ACGT";
    const auto expected = "ACGT";
    auto result = reverse_complement(seq);
    EXPECT_EQ(expected, result);
}


TEST(GraphCleaningTest, simpleTipRemove) {
    const int kmer_size { 21 };
    const std::vector<std::string> sequences {
            //>works well for k=21; part of genome10K.fasta
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC",
            "TGTCATCTAGTTCAACAACCAAAAAAA", //>that's the tip
            "TGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG" //>remaining part
    };
    Graph graph = Graph::create(new BankStrings(sequences), "-kmer-size %d -abundance-min 1 -verbose 0", kmer_size);
    do_graph_clean(graph);

    unsigned int num_non_deleted_nodes { 0 };
    unsigned int num_nodes { 0 };

    GraphIterator<Node> iterNodes = graph.iterator();
    for (iterNodes.first(); !iterNodes.isDone(); iterNodes.next()) {
        num_nodes++;

        if (!graph.isNodeDeleted(*iterNodes)) {
            num_non_deleted_nodes++;
        }
    }

    EXPECT_EQ(num_nodes, 624);
    EXPECT_EQ(num_non_deleted_nodes, 617);
    remove_graph_file();

}


TEST(GenerateStartKmersTest, GenerateOneKmerReturnFirstKCharacters) {
    const std::string sequence { "ACGTGCGATGCAT" };
    const auto k { 4 };
    const auto n { 1 };

    const auto result { generate_start_kmers(sequence, k, n) };
    const std::vector<std::string> expected { "ACGT" };

    EXPECT_EQ(result, expected);
}


TEST(GenerateStartKmersTest, GenerateTwoKmersReturnFirstTwoKmers) {
    const std::string sequence { "ACGTGCGATGCAT" };
    const auto k { 4 };
    const auto n { 2 };

    const auto result { generate_start_kmers(sequence, k, n) };
    const std::vector<std::string> expected { "ACGT", "CGTG" };

    EXPECT_EQ(result, expected);
}


TEST(GenerateStartKmersTest, GenerateMaxPossibleNumKmersReturnWholeSeqAsKmers) {
    const std::string sequence { "ACGTGCGA" };
    const auto k { 4 };
    const auto n { 5 };

    const auto result { generate_start_kmers(sequence, k, n) };
    const std::vector<std::string> expected { "ACGT", "CGTG", "GTGC", "TGCG", "GCGA" };

    EXPECT_EQ(result, expected);
}


TEST(GenerateStartKmersTest, GenerateTooManyKmersReturnWholeSeqAsKmers) {
    const std::string sequence { "ACGTGCGA" };
    const auto k { 4 };
    const auto n { 20 };

    const auto result { generate_start_kmers(sequence, k, n) };
    const std::vector<std::string> expected { "ACGT", "CGTG", "GTGC", "TGCG", "GCGA" };

    EXPECT_EQ(result, expected);
}


TEST(GenerateStartKmersTest, GenerateNoKmersReturnEmptySet) {
    const std::string sequence { "ACGTACGT" };
    const auto k { 4 };
    const auto n { 0 };

    const auto result { generate_start_kmers(sequence, k, n) };
    const std::vector<std::string> expected;

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmersTest, GenerateOneKmerReturnLastKCharacters) {
    const std::string sequence { "ACGTGCGATGCAT" };
    const auto k { 4 };
    const auto n { 1 };

    const auto actual { generate_end_kmers(sequence, k, n) };
    const std::vector<std::string> expected { "GCAT" };

    EXPECT_EQ(expected, actual);
}


TEST(GenerateEndKmersTest, GenerateTwoKmersReturnLastTwoKmers) {
    const std::string sequence { "ACGTGCGATGCAT" };
    const auto k { 4 };
    const auto n { 2 };

    const auto result { generate_end_kmers(sequence, k, n) };
    const std::vector<std::string> expected { "GCAT", "TGCA" };

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmersTest, GenerateMaxPossibleNumKmersReturnWholeSeqAsKmers) {
    const std::string sequence { "ACGTGCGA" };
    const auto k { 4 };
    const auto n { 5 };

    const auto result { generate_end_kmers(sequence, k, n) };
    const std::vector<std::string> expected { "GCGA", "TGCG", "GTGC", "CGTG", "ACGT" };

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmersTest, GenerateTooManyKmersReturnWholeSeqAsKmers) {
    const std::string sequence { "ACGTGCGA" };
    const auto k { 4 };
    const auto n { 20 };

    const auto result { generate_end_kmers(sequence, k, n) };
    const std::vector<std::string> expected { "GCGA", "TGCG", "GTGC", "CGTG", "ACGT" };

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmersTest, SequenceHasRepeatKmersReturnOnlyUniqueKmers) {
    const std::string sequence { "ACGTACGT" };
    const auto k { 4 };
    const auto n { 20 };

    const auto result { generate_end_kmers(sequence, k, n) };
    const std::vector<std::string> expected { "ACGT", "TACG", "GTAC", "CGTA", "ACGT" };

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmersTest, GenerateNoKmersReturnEmptySet) {
    const std::string sequence { "ACGTACGT" };
    const auto k { 4 };
    const auto n { 0 };

    const auto result { generate_end_kmers(sequence, k, n) };
    const std::vector<std::string> expected;

    EXPECT_EQ(result, expected);
}


TEST(QueryAbundance, oneKmer_ReturnOne) {
    Graph graph = Graph::create(new BankStrings("AATGT", NULL), "-kmer-size 5 -abundance-min 1 -verbose 0");
    auto kmer = "AATGT";
    auto node { graph.buildNode(kmer) };
    const auto covg { graph.queryAbundance(node) };

    // We get the neighbors of this real node and make sure it has the neighbours we expect
    EXPECT_EQ(covg, 1);
    remove_graph_file();
}


TEST(QueryAbundance, twoKmers_ReturnTwo) {
    Graph graph = Graph::create(new BankStrings("AATGTAATGT", NULL), "-kmer-size 5 -abundance-min 1 -verbose 0");
    auto kmer = "AATGT";
    auto node { graph.buildNode(kmer) };
    const auto covg { graph.queryAbundance(node) };

    // We get the neighbors of this real node and make sure it has the neighbours we expect
    EXPECT_EQ(covg, 2);
    remove_graph_file();
}


TEST(QueryAbundance, fakeKmer_ReturnZero) {
    Graph graph = Graph::create(new BankStrings("AATGT", NULL), "-kmer-size 5 -abundance-min 1 -verbose 0");
    auto kmer = "CCCCC";
    auto node { graph.buildNode(kmer) };
    const auto covg { graph.queryAbundance(node) };

    // We get the neighbors of this real node and make sure it has the neighbours we expect
    EXPECT_EQ(covg, 0);
    remove_graph_file();
}