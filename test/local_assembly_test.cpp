#include "gtest/gtest.h"
#include <gatb/gatb_core.hpp>
#include <iostream>
#include "local_assembly.h"
#include <cstdio>




const int g_test_kmer_size = 5;


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
}


TEST(GetPathsBetweenTest, OnlyReturnPathsBetweenStartAndEndKmers) {
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
    auto result = get_paths_between("AATGT", end_kmer, tree, graph);

    Paths expected_seqs(seqs.begin(), seqs.end());
    EXPECT_EQ(result, expected_seqs);
}


TEST(DFSTest, SimpleGraphTwoNodes_ReturnSeqPassedIn) {
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
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph);

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(*result.begin(), seq);
}


TEST(DFSTest, SimpleGraphSixNodes_ReturnSeqPassedIn) {
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
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph);

    bool original_seq_found = false;
    // make sure all paths begin and end with correct kmer
    for (auto &path: result) {
        EXPECT_EQ(path.substr(0, g_test_kmer_size), start_kmer);
        EXPECT_EQ(path.substr(path.length() - g_test_kmer_size, path.length()), end_kmer);

        if (path == seq)
            original_seq_found = true;
    }

    EXPECT_TRUE(original_seq_found);
}

TEST(DFSTest, TwoReadsSameSequence_ReturnOneSequence) {
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
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph);

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(*result.begin(), seq1);
}

TEST(DFSTest, TwoReadsOneVariant_ReturnOriginalTwoSequences) {
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
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph);

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

}


TEST(DFSTest, ThreeReadsTwoVariants_ReturnOriginalSequences) {
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
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph);

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
}


TEST(DFSTest, TwoReadsTwoVariants_ReturnOriginalTwoSequencesPlusTwoMosaics) {
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
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph);

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

}


TEST(DFSTest, ThreeReadsOneReverseCompliment_ReturnPathsForStrandOfStartAndEndKmers) {
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
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph);

    // add other expected paths due to variants
    const std::string expected_seq = "ATGTGCA";

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(*result.begin(), expected_seq);

}

TEST(DFSTest, SimpleCycle_ReturnPathsOfLengthsUpToMaxPathLengthCycling) {
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
    auto result = get_paths_between(start_kmer, end_kmer, tree, graph);

    const std::string min_expected_seq = "ATATAT";
    bool is_in = false;

    for (auto &path: result) {
        if (path == min_expected_seq) {
            is_in = true;
        }
    }

    EXPECT_TRUE(is_in);
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


TEST(ReverseCompliment, SingleBase_ReturnCompliment) {
    const auto seq = "A";
    const auto expected = "T";
    auto result = reverse_compliment(seq);
    EXPECT_EQ(expected, result);
}


TEST(ReverseCompliment, TwoBases_ReturnCompliment) {
    const auto seq = "AA";
    const auto expected = "TT";
    auto result = reverse_compliment(seq);
    EXPECT_EQ(expected, result);
}


TEST(ReverseCompliment, AllBases_ReturnCompliment) {
    const auto seq = "ACTGCA";
    const auto expected = "TGCAGT";
    auto result = reverse_compliment(seq);
    EXPECT_EQ(expected, result);
}


TEST(ReverseCompliment, Palindrome_ReturnCompliment) {
    const auto seq = "ACGT";
    const auto expected = "ACGT";
    auto result = reverse_compliment(seq);
    EXPECT_EQ(expected, result);
}

//TEST(LocalAssemblyTest, buildGraphFromRealReads_ExpectRefPathInResults) {
//    const std::string ref_sequence = "TCCTCAAGCACCAGGTACGC";
//    const std::string reads_filepath = "../../test/test_cases/loman_k12_merged_pass.mm2.sorted_1196-1216.fastq";
//    const std::string start_kmer = ref_sequence.substr(0, g_kmer_size);
//    const std::string end_kmer = ref_sequence.substr(ref_sequence.length() - g_kmer_size, ref_sequence.length());
//    const std::string out_path = "../../test/test_cases/local_assembly_test.fa";
//
//    local_assembly(reads_filepath, start_kmer, end_kmer, out_path);
//}


/*
// test if path exists in graph. take all kmers of ref and query each one
TEST(LocalAssemblyTest, buildGraphForAllSlices_writeAllPathsToFile) {
    const std::string meta_file = "/Users/mbhall88/Projects/Pandora_variation/slice_fastq_files/padding_10/ref_seqs_for_slices_padding_10.tsv";

    std::ifstream fin (meta_file);
    std::string line;
    std::string filepath;
    std::string ref_sequence;

    while (std::getline(fin, line)) {
        std::stringstream ss (line);
        ss >> filepath >> ref_sequence;

        std::cout << "Processing " << filepath << "\n";

        const long max_length = ref_sequence.length() + 10;
        const std::string start_kmer = ref_sequence.substr(0, g_kmer_size);
        const std::string end_kmer = ref_sequence.substr(ref_sequence.length() - g_kmer_size, std::string::npos);

        // clear the stringstream
        ss.str(std::string());

        std::ostringstream oss;
        oss << "/Users/mbhall88/Projects/Pandora_variation/slice_fastq_files/padding_10/local_assembly_paths";
        int idx = filepath.rfind('/');
        oss << filepath.substr(idx, filepath.rfind('.') - idx) << ".fa";
        std::string out_path = oss.str();

        local_assembly(filepath, start_kmer, end_kmer, out_path);

        // clear the stringstream
        oss.str(std::string());
    }
}
*/
