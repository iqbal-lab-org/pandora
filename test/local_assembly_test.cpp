#include "gtest/gtest.h"
#include <gatb/gatb_core.hpp>
#include <iostream>
#include "local_assembly.h"
#include <cstdio>
#include <unordered_set>


const int g_test_kmer_size = 5;
const auto test_log_level{logging::trivial::info};


TEST(GetNodeFromGraph, create) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );
    auto kmer = "AATGT";
    Node real_node;
    bool found;
    std::tie(real_node, found) = get_node(kmer, graph);
    EXPECT_TRUE(found);

    // We get the neighbors of this real node and make sure it has the neighbours we expect
    GraphVector<Node> neighbours = graph.successors(real_node);
    EXPECT_EQ(graph.toString(neighbours[0]), "ATGTC");
    remove_graph_file();
}


TEST(GetNodeFromGraph, GivenGraphAndKmer_KmerFoundInGraph) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );
    auto kmer = "AATGT";

    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    auto &result = found;
    EXPECT_TRUE(result);
    remove_graph_file();
}


TEST(GetNodeFromGraph, GivenGraphAndMissingKmer_KmerNotFoundInGraph) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );

    auto kmer = "ACTGT";
    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    auto &result = found;
    EXPECT_FALSE(found);
    remove_graph_file();
}


TEST(GetNodeFromGraph, GivenGraphAndKmer_CorrectNodeReturned) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );
    auto kmer = "AATGT";

    Node node;
    bool found;
    std::tie(node, found) = get_node(kmer, graph);

    auto result = graph.toString(node);
    auto &expected = kmer;

    EXPECT_EQ(expected, result);
    remove_graph_file();
}


TEST(GetNodeFromGraph, GivenGraphAndMissingKmer_CorrectEmptyNodeReturned) {
    Graph graph = Graph::create(
            new BankStrings("AATGTCAGG", NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );
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
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const std::string s1{"AATGTAAGG"};
    const std::string s2{"AATGTCAGG"};
    const std::string s3{"AATGTTAGG"};
    std::vector<std::string> seqs = {s1, s2, s3};

    Graph graph = Graph::create(
            new BankStrings(seqs),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node("AATGT", graph);

    auto tree = DFS(start_node, graph);

    const auto end_kmer{"AGG"};
    auto result = get_paths_between("AATGT", end_kmer, tree, graph, g_max_length);

    Paths expected_seqs(seqs.begin(), seqs.end());
    EXPECT_EQ(result, expected_seqs);
    remove_graph_file();
}


TEST(DFSTest, SimpleGraphTwoNodes_ReturnSeqPassedIn) {
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const auto seq{"ATGCAG"};
    const auto start_kmer{"ATGCA"};
    const auto end_kmer{"TGCAG"};

    const Graph graph = Graph::create(
            new BankStrings(seq, NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = DFS(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_max_length);

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(*result.begin(), seq);
    remove_graph_file();
}


TEST(DFSTest, SimpleGraphSixNodes_ReturnSeqPassedIn) {
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const auto seq{"ATGCAGTACA"};
    const auto start_kmer{"ATGCA"};
    const auto end_kmer{"GTACA"};

    const Graph graph = Graph::create(
            new BankStrings(seq, NULL),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = DFS(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_max_length);

    bool original_seq_found = false;
    // make sure all paths begin and end with correct kmer
    for (auto &path: result) {
        EXPECT_EQ(path.substr(0, g_test_kmer_size), start_kmer);
        EXPECT_EQ(path.substr(path.length() - g_test_kmer_size, path.length()), end_kmer);

        if (path == seq)
            original_seq_found = true;
    }

    EXPECT_TRUE(original_seq_found);
    remove_graph_file();
}

TEST(DFSTest, TwoReadsSameSequence_ReturnOneSequence) {
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const auto seq1{"ATGCAG"};
    const auto seq2{"ATGCAG"};
    std::vector<std::string> seqs = {seq1, seq2};
    const auto start_kmer{"ATGCA"};
    const auto end_kmer{"TGCAG"};

    const Graph graph = Graph::create(
            new BankStrings(seqs),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = DFS(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_max_length);

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(*result.begin(), seq1);
    remove_graph_file();
}

TEST(DFSTest, TwoReadsOneVariant_ReturnOriginalTwoSequences) {
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const std::string seq1{"ATGCAGTACAA"};
    const std::string seq2{"ATGCATTACAA"};
    std::vector<std::string> seqs = {seq1, seq2};
    const auto start_kmer{"ATGCA"};
    const auto end_kmer{"TACAA"};

    const Graph graph = Graph::create(
            new BankStrings(seqs),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = DFS(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_max_length);

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


TEST(DFSTest, ThreeReadsTwoVariants_ReturnOriginalSequences) {
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const std::string seq1{"ATGCAGTACAA"};
    const std::string seq2{"ATGCATTACAA"};
    const std::string seq3{"ATGCACTACAA"};
    std::vector<std::string> seqs = {seq1, seq2, seq3};
    const auto start_kmer{"ATGCA"};
    const auto end_kmer{"TACAA"};

    const Graph graph = Graph::create(
            new BankStrings(seqs),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = DFS(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_max_length);

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


TEST(DFSTest, TwoReadsTwoVariants_ReturnOriginalTwoSequencesPlusTwoMosaics) {
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const std::string seq1{"ATGCAGTACAAGGATAC"};
    const std::string seq2{"ATGCATTACAATGATAC"};
    std::vector<std::string> seqs = {seq1, seq2};
    const auto start_kmer{"ATGCA"};
    const auto end_kmer{"GATAC"};

    const Graph graph = Graph::create(
            new BankStrings(seqs),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = DFS(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_max_length);

    // add other expected paths due to variants
    const std::vector<std::string> expected_seqs = {
            seq1, seq2,
            "ATGCAGTACAATGATAC",
            "ATGCATTACAAGGATAC"
    };

    int original_seq_found = 0;
    for (auto &path: result) {
        EXPECT_EQ(path.substr(0, g_test_kmer_size), start_kmer);
        EXPECT_EQ(path.substr(path.length() - g_test_kmer_size, path.length()), end_kmer);

        if (path.length() == seq1.length()) {
            bool path_in_expected = std::find(expected_seqs.begin(), expected_seqs.end(), path) != expected_seqs.end();
            if (path_in_expected) {
                ++original_seq_found;
            }
        }
    }
    EXPECT_EQ(original_seq_found, expected_seqs.size());
    remove_graph_file();

}


TEST(DFSTest, ThreeReadsOneReverseCompliment_ReturnPathsForStrandOfStartAndEndKmers) {
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const std::string seq1{"ATGTG"};
    const std::string seq2{"TGTGC"};
    const std::string seq3{"TGCAC"};
    std::vector<std::string> seqs = {seq1, seq2, seq3};
    const auto start_kmer{"ATGTG"};
    const auto end_kmer{"GTGCA"};

    const Graph graph = Graph::create(
            new BankStrings(seqs),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = DFS(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_max_length);

    // add other expected paths due to variants
    const std::string expected_seq = "ATGTGCA";

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(*result.begin(), expected_seq);
    remove_graph_file();

}

TEST(DFSTest, SimpleCycle_ReturnPathsOfLengthsUpToMaxPathLengthCycling) {
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const std::string seq1{"ATATATATA"};
    const std::string seq2{"TATAT"};
    std::vector<std::string> seqs = {seq1, seq2};
    const auto start_kmer{"ATATA"};
    const auto end_kmer{"TATAT"};

    const Graph graph = Graph::create(
            new BankStrings(seqs),
            "-kmer-size %d -abundance-min 1 -verbose 0", g_test_kmer_size
    );

    Node start_node;
    bool found;
    std::tie(start_node, found) = get_node(start_kmer, graph);

    auto tree = DFS(start_node, graph);
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph, g_max_length);

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


TEST(hasEndingTest, hasEnding_ReturnTrue) {
    std::string test = "binary";
    std::string ending = "nary";

    EXPECT_TRUE(has_ending(test, ending));
}


TEST(hasEndingTest, doesNotHaveEnding_ReturnFalse) {
    std::string test = "tertiary";
    std::string ending = "nary";

    EXPECT_FALSE(has_ending(test, ending));
}


TEST(hasEndingTest, endingLongerThanQuery_ReturnFalse) {
    std::string test = "ry";
    std::string ending = "nary";

    EXPECT_FALSE(has_ending(test, ending));
}


TEST(FastaWriter, ReadsShorterThanLineWidth_OneReadPerLine) {
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const auto filepath = "TEST.fa";
    const auto header = ">path";
    unsigned long line_width = 90;

    Paths reads = {"ATGATGTTTTTTTTTTCGCATGCAT", "TGCATGCATGCACACACACACACAGCA"};

    write_paths_to_fasta(filepath, reads, line_width);

    std::ifstream in_file(filepath);

    std::string line;

    for (auto &read: reads) {
        std::getline(in_file, line);
        EXPECT_EQ(line, header);

        for (unsigned long i = 0; i < read.length(); i += line_width) {
            std::getline(in_file, line);
            EXPECT_EQ(line, read.substr(i, line_width));
        }
    }

    EXPECT_TRUE(in_file.peek() == std::ifstream::traits_type::eof());
    EXPECT_TRUE(std::remove(filepath) == 0);
}


TEST(FastaWriter, ReadsLongerThanLineWidth_ReadSpreadEvenlyOnLines) {
    logging::core::get()->set_filter(logging::trivial::severity >= test_log_level);

    const auto filepath = "TEST.fa";
    const auto header = ">path";
    unsigned long line_width = 10;

    Paths reads = {"ATGATGTTTTTTTTTTCGCATGCAT", "TGCATGCATGCACACACACACACAGCA"};

    write_paths_to_fasta(filepath, reads, line_width);

    std::ifstream in_file(filepath);

    std::string line;

    for (auto &read: reads) {
        std::getline(in_file, line);
        EXPECT_EQ(line, header);

        for (unsigned long i = 0; i < read.length(); i += line_width) {
            std::getline(in_file, line);
            EXPECT_EQ(line, read.substr(i, line_width));
        }
    }

    EXPECT_TRUE(in_file.peek() == std::ifstream::traits_type::eof());
    EXPECT_TRUE(std::remove(filepath) == 0);
}


TEST(ReverseComplement, SingleBase_ReturnCompliment) {
    const auto seq = "A";
    const auto expected = "T";
    auto result = reverse_complement(seq);
    EXPECT_EQ(expected, result);
}


TEST(ReverseComplement, TwoBases_ReturnCompliment) {
    const auto seq = "AA";
    const auto expected = "TT";
    auto result = reverse_complement(seq);
    EXPECT_EQ(expected, result);
}


TEST(ReverseComplement, AllBases_ReturnCompliment) {
    const auto seq = "ACTGCA";
    const auto expected = "TGCAGT";
    auto result = reverse_complement(seq);
    EXPECT_EQ(expected, result);
}


TEST(ReverseComplement, Palindrome_ReturnCompliment) {
    const auto seq = "ACGT";
    const auto expected = "ACGT";
    auto result = reverse_complement(seq);
    EXPECT_EQ(expected, result);
}

TEST(FileExists, realFile_returnsTrue) {
    const std::string filepath{"../../CMakeLists.txt"};
    const bool exists{file_exists(filepath)};

    EXPECT_TRUE(exists);
}

TEST(FileExists, fakeFile_returnsFalse) {
    const std::string filepath{"../../fake.txt"};
    const bool exists{file_exists(filepath)};

    EXPECT_FALSE(exists);
}


TEST(GraphCleaning, simpleTip_remove) {
    const int kmer_size{21};
    const std::vector<std::string> sequences{
            //>works well for k=21; part of genome10K.fasta
            "CATCGATGCGAGACGCCTGTCGCGGGGAATTGTGGGGCGGACCACGCTCTGGCTAACGAGCTACCGTTTCCTTTAACCTGCCAGACGGTGACCAGGGCCGTTCGGCGTTGCATCGAGCGGTGTCGCTAGCGCAATGCGCAAGATTTTGACATTTACAAGGCAACATTGCAGCGTCCGATGGTCCGGTGGCCTCCAGATAGTGTCCAGTCGCTCTAACTGTATGGAGACCATAGGCATTTACCTTATTCTCATCGCCACGCCCCAAGATCTTTAGGACCCAGCATTCCTTTAACCACTAACATAACGCGTGTCATCTAGTTCAACAACC",
            "TGTCATCTAGTTCAACAACCAAAAAAA", //>that's the tip
            "TGTCATCTAGTTCAACAACCGTTATGCCGTCCGACTCTTGCGCTCGGATGTCCGCAATGGGTTATCCCTATGTTCCGGTAATCTCTCATCTACTAAGCGCCCTAAAGGTCGTATGGTTGGAGGGCGGTTACACACCCTTAAGTACCGAACGATAGAGCACCCGTCTAGGAGGGCGTGCAGGGTCTCCCGCTAGCTAATGGTCACGGCCTCTCTGGGAAAGCTGAACAACGGATGATACCCATACTGCCACTCCAGTACCTGGGCCGCGTGTTGTACGCTGTGTATCTTGAGAGCGTTTCCAGCAGATAGAACAGGATCACATGTACATG" //>remaining part
    };
    Graph graph = Graph::create(
            new BankStrings(sequences),
            "-kmer-size %d -abundance-min 1 -verbose 0", kmer_size
    );
    do_graph_clean(graph);

    unsigned int num_non_deleted_nodes{0};
    unsigned int num_nodes{0};

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


TEST(GenerateStartKmers, GenerateOneKmer_ReturnFirstKCharacters) {
    const std::string sequence{"ACGTGCGATGCAT"};
    const auto k{4};
    const auto n{1};

    const auto result{generate_start_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected{"ACGT"};

    EXPECT_EQ(result, expected);
}


TEST(GenerateStartKmers, GenerateTwoKmers_ReturnFirstTwoKmers) {
    const std::string sequence{"ACGTGCGATGCAT"};
    const auto k{4};
    const auto n{2};

    const auto result{generate_start_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected{"ACGT", "CGTG"};

    EXPECT_EQ(result, expected);
}


TEST(GenerateStartKmers, GenerateMaxPossibleNumKmers_ReturnWholeSeqAsKmers) {
    const std::string sequence{"ACGTGCGA"};
    const auto k{4};
    const auto n{5};

    const auto result{generate_start_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected{"ACGT", "CGTG", "GTGC", "TGCG", "GCGA"};

    EXPECT_EQ(result, expected);
}


TEST(GenerateStartKmers, GenerateTooManyKmers_ReturnWholeSeqAsKmers) {
    const std::string sequence{"ACGTGCGA"};
    const auto k{4};
    const auto n{20};

    const auto result{generate_start_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected{"ACGT", "CGTG", "GTGC", "TGCG", "GCGA"};

    EXPECT_EQ(result, expected);
}


TEST(GenerateStartKmers, SequenceHasRepeatKmers_ReturnOnlyUniqueKmers) {
    const std::string sequence{"ACGTACGT"};
    const auto k{4};
    const auto n{20};

    const auto result{generate_start_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected{"ACGT", "CGTA", "GTAC", "TACG"};

    EXPECT_EQ(result, expected);
}


TEST(GenerateStartKmers, GenerateNoKmers_ReturnEmptySet) {
    const std::string sequence{"ACGTACGT"};
    const auto k{4};
    const auto n{0};

    const auto result{generate_start_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected;

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmers, GenerateOneKmer_ReturnLastKCharacters) {
    const std::string sequence{"ACGTGCGATGCAT"};
    const auto k{4};
    const auto n{1};

    const auto result{generate_end_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected{"GCAT"};

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmers, GenerateTwoKmers_ReturnLastTwoKmers) {
    const std::string sequence{"ACGTGCGATGCAT"};
    const auto k{4};
    const auto n{2};

    const auto result{generate_end_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected{"GCAT", "TGCA"};

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmers, GenerateMaxPossibleNumKmers_ReturnWholeSeqAsKmers) {
    const std::string sequence{"ACGTGCGA"};
    const auto k{4};
    const auto n{5};

    const auto result{generate_end_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected{"ACGT", "CGTG", "GTGC", "TGCG", "GCGA"};

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmers, GenerateTooManyKmers_ReturnWholeSeqAsKmers) {
    const std::string sequence{"ACGTGCGA"};
    const auto k{4};
    const auto n{20};

    const auto result{generate_end_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected{"ACGT", "CGTG", "GTGC", "TGCG", "GCGA"};

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmers, SequenceHasRepeatKmers_ReturnOnlyUniqueKmers) {
    const std::string sequence{"ACGTACGT"};
    const auto k{4};
    const auto n{20};

    const auto result{generate_end_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected{"ACGT", "CGTA", "GTAC", "TACG"};

    EXPECT_EQ(result, expected);
}


TEST(GenerateEndKmers, GenerateNoKmers_ReturnEmptySet) {
    const std::string sequence{"ACGTACGT"};
    const auto k{4};
    const auto n{0};

    const auto result{generate_end_kmers(sequence, k, n)};
    const std::unordered_set<std::string> expected;

    EXPECT_EQ(result, expected);
}


TEST(LocalAssemblyTest, passFakeFastqPath_dontRaiseError) {
    const std::string filepath{"FAKE.fakeq"};
    std::string start_kmer{"ATGATGATG"};
    std::string end_kmer{"ATGATGATG"};
    const std::string out_path{"../../test/test_cases/fake.fa"};
    const int k{9};
    const int max_len{30};
    local_assembly(filepath, start_kmer, end_kmer, out_path, k, max_len);
}


TEST(LocalAssemblyTest, twoIdenticalReads_onePath) {
    const std::string filepath{"../../test/test_cases/local_assembly1.fa"};
    std::string start_kmer{"ATGCGCTGA"};
    std::string end_kmer{"AGTCGGACT"};
    const std::string out_path{"../../test/test_cases/local_assembly1_paths.fa"};
    const int k{9};
    const int max_len{30};
    const bool clean{false};
    local_assembly(filepath, start_kmer, end_kmer, out_path, k, max_len, 1, clean);

    const std::unordered_set<std::string> expected{"ATGCGCTGAGAGTCGGACT"};
    std::unordered_set<std::string> result;

    // read paths file  back in and store all paths in set
    std::ifstream fin{out_path};
    std::string line;

    while (std::getline(fin, line)) {
        if (line[0] == '>') {
            line.clear();
        } else {
            result.insert(line);
            line.clear();
        }
    }

    EXPECT_EQ(result, expected);
    remove(out_path.c_str());
}

TEST(LocalAssemblyTest, twoIdenticalOneSoloReadsMinCovgOne_twoPaths) {
    const std::string filepath{"../../test/test_cases/local_assembly2.fa"};
    std::string start_kmer{"ATGCGCTGA"};
    std::string end_kmer{"AGTCGGACT"};
    const std::string out_path{"../../test/test_cases/local_assembly2_paths.fa"};
    const int k{9};
    const int max_len{30};
    const int min_coverage{1};
    const bool clean{false};
    local_assembly(filepath, start_kmer, end_kmer, out_path, k, max_len, 1, clean, min_coverage);

    const std::unordered_set<std::string> expected{"ATGCGCTGATAGTCGGACT", "ATGCGCTGAGAGTCGGACT"};
    std::unordered_set<std::string> result;

    // read paths file  back in and store all paths in set
    std::ifstream fin{out_path};
    std::string line;

    while (std::getline(fin, line)) {
        if (line[0] == '>') {
            line.clear();
        } else {
            result.insert(line);
            line.clear();
        }
    }

    EXPECT_EQ(result, expected);
    remove(out_path.c_str());
}


TEST(LocalAssemblyTest, twoIdenticalOneSoloReadsMinCovgTwo_onePath) {
    const std::string filepath{"../../test/test_cases/local_assembly2.fa"};
    std::string start_kmer{"ATGCGCTGA"};
    std::string end_kmer{"AGTCGGACT"};
    const std::string out_path{"../../test/test_cases/local_assembly2_paths.fa"};
    const int k{9};
    const int max_len{30};
    const int min_coverage{2};
    const bool clean{false};
    local_assembly(filepath, start_kmer, end_kmer, out_path, k, max_len, 1, clean, min_coverage);

    const std::unordered_set<std::string> expected{"ATGCGCTGAGAGTCGGACT"};
    std::unordered_set<std::string> result;

    // read paths file  back in and store all paths in set
    std::ifstream fin{out_path};
    std::string line;

    while (std::getline(fin, line)) {
        if (line[0] == '>') {
            line.clear();
        } else {
            result.insert(line);
            line.clear();
        }
    }

    EXPECT_EQ(result, expected);
    remove(out_path.c_str());
}

TEST(LocalAssemblyTestUsingSequenceVector, twoIdenticalReads_onePath) {
    const std::vector<std::string> sequences{"ATGCGCTGAGAGTCGGACT", "ATGCGCTGAGAGTCGGACT"};
    std::string start_kmer{"ATGCGCTGA"};
    std::string end_kmer{"AGTCGGACT"};
    const std::string out_path{"../../test/test_cases/local_assembly1_paths.fa"};
    const int k{9};
    const int max_len{30};
    const bool clean{false};
    local_assembly(sequences, start_kmer, end_kmer, out_path, k, max_len, 1, clean);

    const std::unordered_set<std::string> expected{"ATGCGCTGAGAGTCGGACT"};
    std::unordered_set<std::string> result;

    // read paths file  back in and store all paths in set
    std::ifstream fin{out_path};
    std::string line;

    while (std::getline(fin, line)) {
        if (line[0] == '>') {
            line.clear();
        } else {
            result.insert(line);
            line.clear();
        }
    }

    EXPECT_EQ(result, expected);
    remove(out_path.c_str());
}

//TEST(LocalAssemblyTest, debug) {
//    const std::string filepath {"../../test/test_cases/GC00000008_17.2568-2577.fa"};
//    std::string start_kmer = "CGGGCGGACGC";
//    std::string end_kmer = "CAGCGCAGGAC";
//    const std::string out_path {"../../test/test_cases/GC00000008_17.2568-2577_local_assembly_test.fa"};
//    local_assembly(filepath, start_kmer, end_kmer, out_path, 11, 150, 49);
//}

//TEST(LocalAssemblyTest, buildGraphFromRealReads_ExpectRefPathInResults) {
//    const std::string ref_sequence = "TCCTCAAGCACCAGGTACGC";
//    const std::string reads_filepath = "../../test/test_cases/loman_k12_merged_pass.mm2.sorted_1196-1216.fastq";
//    const std::string start_kmer = ref_sequence.substr(0, g_kmer_size);
//    const std::string end_kmer = ref_sequence.substr(ref_sequence.length() - g_kmer_size, ref_sequence.length());
//    const std::string out_path = "../../test/test_cases/local_assembly_test.fa";
//
//    local_assembly(reads_filepath, start_kmer, end_kmer, out_path);
//}


// test if path exists in graph. take all kmers of ref and query each one
// TEST(LocalAssemblyTest, buildGraphForAllSlices_writeAllPathsToFile) {
//     const std::string meta_file = "/Users/mbhall88/Projects/Pandora_variation/slice_fastq_files/padding_10/subsample_30x/ref_seqs_for_slices_padding_10.tsv";
//
//     std::ifstream fin (meta_file);
//     std::string line;
//     std::string filepath;
//     std::string ref_sequence;
//
//     while (std::getline(fin, line)) {
//         std::stringstream ss (line);
//         ss >> filepath >> ref_sequence;
//
//         std::cout << "Processing " << filepath << "\n";
//
//         const unsigned long max_length = ref_sequence.length() + 30;
//         auto start_kmer = ref_sequence.substr(0, g_local_assembly_kmer_size);
//         auto end_kmer = ref_sequence.substr(ref_sequence.length() - g_local_assembly_kmer_size, std::string::npos);
//
//         // clear the stringstream
//         ss.str(std::string());
//
//         std::ostringstream oss;
//         oss << "/Users/mbhall88/Projects/Pandora_variation/slice_fastq_files/padding_10/subsample_30x/local_assembly_paths_covg1";
//         auto idx = filepath.rfind('/');
//         oss << filepath.substr(idx, filepath.rfind('.') - idx) << "_K" << g_local_assembly_kmer_size << ".fa";
//         std::string out_path = oss.str();
//         const bool clean {false}   ;
//         const int min_coverage {1};
//
//         local_assembly(filepath, start_kmer, end_kmer, out_path, g_local_assembly_kmer_size, max_length, clean,
//                        min_coverage);
//
//         // clear the stringstream
//         oss.str(std::string());
//     }
// }

TEST(QueryAbundance, oneKmer_ReturnOne) {
    Graph graph = Graph::create(
            new BankStrings("AATGT", NULL),
            "-kmer-size 5 -abundance-min 1 -verbose 0"
    );
    auto kmer = "AATGT";
    auto node{graph.buildNode(kmer)};
    const auto covg{graph.queryAbundance(node)};

    // We get the neighbors of this real node and make sure it has the neighbours we expect
    EXPECT_EQ(covg, 1);
    remove_graph_file();
}


TEST(QueryAbundance, twoKmers_ReturnTwo) {
    Graph graph = Graph::create(
            new BankStrings("AATGTAATGT", NULL),
            "-kmer-size 5 -abundance-min 1 -verbose 0"
    );
    auto kmer = "AATGT";
    auto node{graph.buildNode(kmer)};
    const auto covg{graph.queryAbundance(node)};

    // We get the neighbors of this real node and make sure it has the neighbours we expect
    EXPECT_EQ(covg, 2);
    remove_graph_file();
}


TEST(QueryAbundance, fakeKmer_ReturnZero) {
    Graph graph = Graph::create(
            new BankStrings("AATGT", NULL),
            "-kmer-size 5 -abundance-min 1 -verbose 0"
    );
    auto kmer = "CCCCC";
    auto node{graph.buildNode(kmer)};
    const auto covg{graph.queryAbundance(node)};

    // We get the neighbors of this real node and make sure it has the neighbours we expect
    EXPECT_EQ(covg, 0);
    remove_graph_file();
}