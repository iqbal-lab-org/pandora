#include "gtest/gtest.h"
#include "test_macro.cpp"

#include "interval.h"
#include "prg/path.h"
#include "kmergraph.h"
#include "kmernode.h"
#include "localPRG.h"
#include <stdint.h>
#include <cmath>
#include <deque>
#include "fatal_error.h"
#include "test_helpers_containers.h"
#include "test_helpers.h"

using namespace prg;

TEST(KmerGraphTest, add_node)
{
    // add node and check it's there
    KmerGraph kg;
    uint32_t sample_id = 0;

    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    uint j = 1;
    EXPECT_EQ(j, kg.nodes.size());
    EXPECT_EQ(p, kg.nodes[0]->path);
    j = 0;
    EXPECT_EQ(j, kg.nodes[0]->id);
    EXPECT_EQ(j, kg.nodes[0]->num_AT);

    // add node another time and expect nothing to happen
    kg.add_node(p);
    j = 1;
    EXPECT_EQ(j, kg.nodes.size());
    EXPECT_EQ(p, kg.nodes[0]->path);
    j = 0;
    EXPECT_EQ(j, kg.nodes[0]->id);
    EXPECT_EQ(j, kg.nodes[0]->num_AT);

    // add a second node and check gets next id
    d = { Interval(1, 4) };
    p.initialize(d);
    kg.add_node(p);
    j = 2;
    EXPECT_EQ(j, kg.nodes.size());
    EXPECT_EQ(p, kg.nodes[1]->path);
    j = 0;
    EXPECT_EQ(j, kg.nodes[0]->id);
    EXPECT_EQ(j, kg.nodes[0]->num_AT);
    j = 1;
    EXPECT_EQ(j, kg.nodes[1]->id);
}

TEST(KmerGraphTest, add_node_with_kh)
{
    // add node and check it's there
    KmerGraph kg;
    uint32_t sample_id = 0;

    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p;
    p.initialize(d);
    uint64_t kh = 469;
    kg.add_node_with_kh(p, kh);
    uint j = 1;
    EXPECT_EQ(j, kg.nodes.size());
    EXPECT_EQ(p, kg.nodes[0]->path);
    j = 0;
    EXPECT_EQ(j, kg.nodes[0]->id);
    EXPECT_EQ(j, kg.nodes[0]->num_AT);

    EXPECT_EQ(kh, kg.nodes[0]->khash);
}

TEST(KmerGraphTest, add_edge)
{
    // add edge and check it's there
    KmerGraph kg;

    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    uint j = 2;
    EXPECT_EQ(j, kg.nodes.size());

    kg.add_edge(n1, n2);
    kg.add_edge(n1, n2);

    d = { Interval(4, 7) };
    p3.initialize(d);
    auto n3 = kg.add_node(p3);
    kg.add_edge(n1, n3);
    j = 2;
    EXPECT_EQ(j, kg.nodes[0]->out_nodes.size());
    j = 1;
    EXPECT_EQ(j, kg.nodes[1]->in_nodes.size());
    EXPECT_EQ(j, kg.nodes[2]->in_nodes.size());
    j = 0;
    EXPECT_EQ(j, kg.nodes[1]->out_nodes.size());
    EXPECT_EQ(j, kg.nodes[0]->in_nodes.size());

    // repeat and nothing should happen
    kg.add_edge(n1, n3);
    j = 2;
    EXPECT_EQ(j, kg.nodes[0]->out_nodes.size());
    j = 1;
    EXPECT_EQ(j, kg.nodes[1]->in_nodes.size());
    EXPECT_EQ(j, kg.nodes[2]->in_nodes.size());
    j = 0;
    EXPECT_EQ(j, kg.nodes[1]->out_nodes.size());
    EXPECT_EQ(j, kg.nodes[0]->in_nodes.size());
}

TEST(KmerGraphTest, clear)
{
    KmerGraph kg;
    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);
    uint j = 2;
    EXPECT_EQ(j, kg.nodes.size());

    kg.clear();
    j = 0;
    EXPECT_EQ(j, kg.nodes.size());

    n1 = kg.add_node(p1);
    n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);
    j = 2;
    EXPECT_EQ(j, kg.nodes.size());
}

TEST(KmerGraphTest, equals)
{
    KmerGraph kg1, kg2;
    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kg1.add_node(p1);
    auto m1 = kg2.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kg1.add_node(p2);
    auto m2 = kg2.add_node(p2);
    kg1.add_edge(n1, n2);
    kg2.add_edge(m1, m2);

    d = { Interval(2, 5) };
    p3.initialize(d);
    auto m3 = kg2.add_node(p3);

    // same as themselves, different if different numbers of nodes
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1 == kg2), false);
    EXPECT_EQ((kg2 == kg1), false);

    auto n3 = kg1.add_node(p3);
    kg2.add_edge(m1, m3);

    // same as themselves, different if different numbers of edges
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1 == kg2), false);
    EXPECT_EQ((kg2 == kg1), false);

    kg1.add_edge(n2, n3);

    // same as themselves, different if edges in different places
    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
    EXPECT_EQ((kg1 == kg2), false);
    EXPECT_EQ((kg2 == kg1), false);
}

TEST(KmerGraphTest, copy)
{
    KmerGraph kg1;
    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kg1.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kg1.add_node(p2);
    kg1.add_edge(n1, n2);

    KmerGraph kg2(kg1);

    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
}

TEST(KmerGraphTest, assign)
{
    KmerGraph kg1;
    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    auto n = kg1.add_node(p);
    d = { Interval(0, 3) };
    prg::Path p1, p2, p3;
    p1.initialize(d);
    auto n1 = kg1.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kg1.add_node(p2);
    kg1.add_edge(n1, n2);
    d = { Interval(11, 14) };
    p3.initialize(d);
    auto n3 = kg1.add_node(p3);
    kg1.add_edge(n1, n3);
    d = { Interval(15, 18) };
    p1.initialize(d);
    n1 = kg1.add_node(p1);
    kg1.add_edge(n2, n1);
    d = { Interval(20, 20) };
    p.initialize(d);
    n = kg1.add_node(p);
    kg1.add_edge(n1, n);
    kg1.add_edge(n3, n);

    KmerGraph kg2 = kg1;

    EXPECT_EQ(kg1, kg1);
    EXPECT_EQ(kg2, kg2);
}

TEST(KmerGraphTest, sort_topologically)
{
    std::vector<KmerNodePtr> exp_sorted_nodes;
    KmerNodePtr n;

    KmerGraph kg;
    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13) };
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = { Interval(0, 1), Interval(19, 20), Interval(23, 24) };
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = { Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);
    d = { Interval(24, 24) };
    p.initialize(d);
    n = kg.add_node(p);
    exp_sorted_nodes.push_back(n);

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[0], kg.nodes[2]);
    kg.add_edge(kg.nodes[0], kg.nodes[3]);
    kg.add_edge(kg.nodes[1], kg.nodes[4]);
    kg.add_edge(kg.nodes[2], kg.nodes[5]);
    kg.add_edge(kg.nodes[3], kg.nodes[6]);
    kg.add_edge(kg.nodes[4], kg.nodes[6]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);

    // for each node, outnodes are further along vector
    std::set<KmerNodePtr>::iterator it;
    uint i = 0;
    for (auto c = kg.sorted_nodes.begin(); c != kg.sorted_nodes.end(); ++c) {
        for (const auto& d : (*c)->out_nodes) {
            it = c;
            ++it;
            while ((*it)->path != d.lock()->path and it != kg.sorted_nodes.end()) {
                it++;
            }
            EXPECT_EQ((it != kg.sorted_nodes.end()), true);
        }
        EXPECT_EQ(exp_sorted_nodes[i], *c);
        i++;
    }
}

TEST(KmerGraphTest, check)
{
    KmerGraph kg;
    std::deque<Interval> d = { Interval(0, 0) };
    prg::Path p;
    p.initialize(d);
    kg.add_node(p);
    d = { Interval(0, 1), Interval(4, 5), Interval(8, 9) };
    p.initialize(d);
    kg.add_node(p);
    d = { Interval(0, 1), Interval(4, 5), Interval(12, 13) };
    p.initialize(d);
    kg.add_node(p);
    d = { Interval(0, 1), Interval(19, 20), Interval(23, 24) };
    p.initialize(d);
    kg.add_node(p);
    d = { Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kg.add_node(p);
    d = { Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24) };
    p.initialize(d);
    kg.add_node(p);
    d = { Interval(24, 24) };
    p.initialize(d);
    kg.add_node(p);

    kg.add_edge(kg.nodes[0], kg.nodes[1]);
    kg.add_edge(kg.nodes[0], kg.nodes[2]);
    kg.add_edge(kg.nodes[0], kg.nodes[3]);
    kg.add_edge(kg.nodes[1], kg.nodes[4]);
    kg.add_edge(kg.nodes[2], kg.nodes[5]);
    kg.add_edge(kg.nodes[3], kg.nodes[6]);
    kg.add_edge(kg.nodes[4], kg.nodes[6]);
    kg.add_edge(kg.nodes[5], kg.nodes[6]);

    kg.sorted_nodes = { kg.nodes[0], kg.nodes[1], kg.nodes[2], kg.nodes[3], kg.nodes[4],
        kg.nodes[5], kg.nodes[6] };
    kg.check();
    kg.sorted_nodes = { kg.nodes[0], kg.nodes[1], kg.nodes[4], kg.nodes[3], kg.nodes[2],
        kg.nodes[5], kg.nodes[6] };
    kg.check();
    kg.sorted_nodes = { kg.nodes[6], kg.nodes[5], kg.nodes[0], kg.nodes[3], kg.nodes[2],
        kg.nodes[1], kg.nodes[4] };
    // There is no way to expect death here since kg.sorted_nodes are always sorted now
    // EXPECT_DEATH(kg.check(), "");
    kg.check();
}

TEST(KmerGraphTest, remove_shortcut_edges)
{
    auto index = std::make_shared<Index>();
    uint32_t w = 14, k = 15;

    auto s
        = " 5 "
          "CATGCGCCAGGGCGCCAATCATGCGGGCGCTCATCAGGGCGAACATCGAATAAGACCGGGTTGCGGCGAGGCAGGA"
          "AAACGCGAGGATCAGCATCAGCCCGACCAGCAGCGCCTTGCGGGAAATACGCGCCGGCATTGCGCCGGAAAGCAGA"
          "GCCGCCAGGGCGCCTACCCAGCCATAGGCGGTGACGGCGAGGCCCACGCCGGATTCCGTCTGGTGAAAATCCGCCG"
          "CCAGGGCGTTGAGCATGCCCACCGGCGCCAGTTCGCTGGTGACGATCGAAAAGGCGCAGATCCCGAGCGCAACGAC"
          "GGCAGTCCAGACGCGCGCCGGCGCCGGGTGGAGGGGTAAAGCAATCTCTTTCAT 6  6  7  8 "
          "AAAGGCGCAGATCCCGAGTGCAACGACGGCTATCCAGACGCGCGCCGGCGCCGGGTGGAGGGGTAA 7 "
          "AGCAATCTCTTTCAT 5 ATCAGGC 9 C 10 G 9 TATCCTTAGGAAAGG 11 T 12 A 11 GCGTTCCG "
          "13  15 T 16 C 15 GCGGTGCACG 17 A 18 G 17  14  19 CA 20 CG 19 CGGTACACGG 13 "
          "ACGTTCAGGTGA 21  23 T 24 G 23 GAGAGAGCAG 25 GCGACCG 26 GCGACCA 26 ACGACCA "
          "26 GCGATCG 25  22 GGAGAGCACAGGCGATCG 22 GGAGAGAGCAAGCGACCG 22 "
          "GGGGAGAGCAGGTGACCG 21 GATGGCCTG 27 T 28 G 27 TTGTCTCCG 29  31 CGAA 32 TGAG "
          "32 CGAG 31 TGGCGTGCAGTATCATCCC 33 TT 34 TG 34 CG 33 CAAAATTGATAAAAAAGAGC 35 "
          "A 36 G 35 GAAAACGGAG 37 AGCTG 38 GGCCG 38 AGCTA 38 AGCCG 38 ATCCG 37 "
          "TTTTCCATA 39  41 AAC 42 CAT 42 AAT 41 GGAAAAGAG 40  43 T 44 A 44 C 43 "
          "ATGGAAAATAG 39  30  45 CGAA 46 CGAG 45 TGGCGTGCAGTATCATCCCTGCGAAA 47 A 48 C "
          "47 TGATAAAAAAGAGCGGAAAACGGAG 49 AGCT 50 AGCC 50 AGTC 50 GGCC 49 GTTTTCCATA "
          "51 T 52 A 52 C 51 ATGGAAAA 53 TAG 54 GAG 53  30  55  57 CA 58 CG 57 "
          "AGTGGCGTG 59 T 60 C 59  56 CAAGTGGTGTGC 55 AGTATCATCCCTG 61 T 62 C 61 "
          "GAAACTGA 63 T 64 A 63 AAAAAATAGCGGAAAACGGA 65 GAGT 66 TAGC 65 "
          "CGTTTTCCATAAATGGAAAACAG 30 "
          "CGAGTGGCGTGCAGTATCATCCCTGCGAAAATGATAAAAAAGAGTGGAAAACGGATAGCCGTTTTCCATAAATGGA"
          "AAA 67 TAG 68 CAG 67  30  69 CGAA 70 CGAG 69 TGGCGTGCAGTA 30  71 CGAATGGC "
          "72 CGAGTGGT 72 CGAGTGGC 71 GTGCAGTATCATCCCTGCGAAACTGATAAAAAAGAGC 73 A 74 G "
          "73 GAAAACGGAGAGCCGTTTTCCATAAA 75 T 76 C 75 GGAAAAGAG 29 ";
    auto l1 = std::make_shared<LocalPRG>(LocalPRG(1, "Cluster_6369", s));
    l1->minimizer_sketch(index.get(), w, k);

    s = "TTATAAAGTTCTGCAAATGGCGCCATCAAAGCGCCATTGACAGAGTTTTATTTCAATCACCTTTTTCGAGGTATCAAA"
        "AATCACGGGGTTTTAATCCCTTCCTCCAATAAGTACCAGTTTAATATTCTGAATGCCCGTCACGGGGCAACATAACCA"
        "CAGAGCCTTGCGGGGTGGGTCTATGGGGTAGGCAGTAATGCTTTCACTCTGTGGGCTGCTTTTATCCGCGTGAACTTA"
        "GGCTCACCACCGAAAGGAAAAGCA";
    auto l2 = std::make_shared<LocalPRG>(LocalPRG(1, "Cluster_15213", s));
    l2->minimizer_sketch(index.get(), w, k);
}

TEST(KmerGraphTest, save)
{
    auto l = std::make_shared<LocalPRG>(LocalPRG(1, "test localPRG", "ACGT"));

    KmerGraph kg;
    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);
    EXPECT_EQ((uint)0, kg.nodes[0]->num_AT);

    // kg.save("../test/test_cases/kmergraph_test.gfa");
    kg.save("kmergraph_test.gfa", l);
}

TEST(KmerGraphTest, save_no_prg)
{
    KmerGraph kg;
    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);
    EXPECT_EQ((uint)0, kg.nodes[0]->num_AT);

    // kg.save("../test/test_cases/kmergraph_test.gfa");
    kg.save("kmergraph_test2.gfa");
}

TEST(KmerGraphTest, load)
{
    KmerGraph kg, read_kg;
    std::deque<Interval> d = { Interval(0, 3) };
    prg::Path p1, p2;
    p1.initialize(d);
    auto n1 = kg.add_node(p1);
    d = { Interval(1, 4) };
    p2.initialize(d);
    auto n2 = kg.add_node(p2);
    kg.add_edge(n1, n2);

    // read_kg.load("../test/test_cases/kmergraph_test.gfa");
    read_kg.load("kmergraph_test2.gfa");
    EXPECT_EQ(kg, read_kg);
}

TEST(KmerGraphTest, load_prg_FatalRuntimeError)
{
    KmerGraph read_kg;
    ASSERT_EXCEPTION(
        read_kg.load("kmergraph_test.gfa"), FatalRuntimeError, "Error reading GFA");
}
