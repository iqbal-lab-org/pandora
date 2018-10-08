#include <cstdint>
#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "pangenome/ns.cpp"
#include "pangenome/pannode.h"
#include "pangenome/pansample.h"
#include "pangenome/pangraph.h"
#include "pangenome/panread.h"
#include "minihit.h"
#include "localPRG.h"


using namespace pangenome;

TEST(PangenomeNodeTest, create) {

    Node pn(4, 3, "3");
    uint32_t j = 3;
    EXPECT_EQ(j, pn.node_id);
    EXPECT_EQ((uint) 4, pn.prg_id);
    EXPECT_EQ("3", pn.name);
    EXPECT_EQ((uint) 1, pn.covg);
    EXPECT_EQ((uint) 0, pn.reads.size());
    EXPECT_EQ((uint) 0, pn.samples.size());
}

TEST(PangenomeNodeTest, get_name) {
    Node pn1(3, 3, "3");
    Node pn2(2, 2, "2");
    Node pn3(2, 4, "2");

    EXPECT_EQ(pn1.get_name(), "3");
    EXPECT_EQ(pn2.get_name(), "2");
    EXPECT_EQ(pn3.get_name(), "2.4");
}

TEST(PangenomeNodeTest, add_path) {
    Node pn1(3, 3, "3");
    vector<KmerNodePtr> kmp;
    pn1.add_path(kmp);

    KmerGraph kg;
    deque<Interval> d = {Interval(0, 0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(8, 9), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(4, 5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4, 5), Interval(12, 13), Interval(16, 16), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0, 1), Interval(19, 20), Interval(23, 24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24, 24)};
    p.initialize(d);
    kg.add_node(p);
    EXPECT_EQ((uint) 7, kg.nodes.size());

    pn1.kmer_prg = kg;
    EXPECT_EQ((uint) 7, pn1.kmer_prg.nodes.size());
    pn1.kmer_prg.sort_topologically();
    EXPECT_EQ((uint) 7, pn1.kmer_prg.sorted_nodes.size());
    kmp = {pn1.kmer_prg.sorted_nodes[0], pn1.kmer_prg.sorted_nodes[3], pn1.kmer_prg.sorted_nodes[4],
           pn1.kmer_prg.sorted_nodes[6]};
    pn1.add_path(kmp);
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[0]->covg[0]);
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[1]->covg[0]);
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[2]->covg[0]);
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[3]->covg[0]);
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[4]->covg[0]);
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[5]->covg[0]);
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[6]->covg[0]);
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[0]->covg[1]);
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[1]->covg[1]);
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[2]->covg[1]);
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[3]->covg[1]);
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[4]->covg[1]);
    EXPECT_EQ((uint) 0, pn1.kmer_prg.sorted_nodes[5]->covg[1]);
    EXPECT_EQ((uint) 1, pn1.kmer_prg.sorted_nodes[6]->covg[1]);
}

TEST(PangenomeNodeTest, get_read_overlap_coordinates) {
    Node pn(3, 3, "3");
    pangenome::ReadPtr pr;
    MinimizerHits mhits;

    Minimizer m;
    deque<Interval> d;
    Path p;
    MiniRecord *mr;

    // read1
    m = Minimizer(0, 1, 6, 0); // kmer, start, end, strand
    d = {Interval(7, 8), Interval(10, 14)};
    p.initialize(d);
    mr = new MiniRecord(0, p, 0, 0);
    mhits.add_hit(1, m, mr); // read 1

    m = Minimizer(0, 0, 5, 0);
    d = {Interval(6, 10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    mhits.add_hit(1, m, mr);

    d = {Interval(6, 10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    mhits.add_hit(1, m, mr);

    mhits.sort();
    pr = make_shared<pangenome::Read>(1);
    pr->add_hits(3, mhits.hits);
    pn.reads.insert(pr);
    mhits.clear();

    //read 2
    m = Minimizer(0, 2, 7, 1);
    d = {Interval(6, 10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    mhits.add_hit(2, m, mr);

    m = Minimizer(0, 5, 10, 1);
    d = {Interval(6, 10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0, p, 0, 0);
    mhits.add_hit(2, m, mr);

    mhits.sort();
    pr = make_shared<pangenome::Read>(2);
    pr->add_hits(3, mhits.hits);
    pn.reads.insert(pr);
    mhits.clear();

    delete mr;

    vector<vector<uint32_t>> read_overlap_coordinates;
    pn.get_read_overlap_coordinates(read_overlap_coordinates);
    vector<vector<uint32_t>> expected_read_overlap_coordinates = {{1, 0, 6,  1},
                                                                  {2, 2, 10, 0}};
    for (const auto &coord : read_overlap_coordinates) {
        if (coord[0] == 1) {
            EXPECT_ITERABLE_EQ(vector<uint32_t>, expected_read_overlap_coordinates[0], coord);
        } else {
            EXPECT_ITERABLE_EQ(vector<uint32_t>, expected_read_overlap_coordinates[1], coord);
        }
    }
}

/*
TEST(PangenomeNodeTest,output_samples)
{
    Node pn1(3,3,"three");
    vector<KmerNodePtr> kmp;

    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p; 
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);
    kg.sort_topologically();
    kg.add_edge(kg.sorted_nodes[0],kg.sorted_nodes[1]);
    kg.add_edge(kg.sorted_nodes[1],kg.sorted_nodes[4]);
    kg.add_edge(kg.sorted_nodes[0],kg.sorted_nodes[2]);
    kg.add_edge(kg.sorted_nodes[2],kg.sorted_nodes[5]);
    kg.add_edge(kg.sorted_nodes[0],kg.sorted_nodes[3]);
    kg.add_edge(kg.sorted_nodes[3],kg.sorted_nodes[6]);
    kg.add_edge(kg.sorted_nodes[4],kg.sorted_nodes[6]);
    kg.add_edge(kg.sorted_nodes[5],kg.sorted_nodes[6]);

    EXPECT_EQ((uint)7, kg.nodes.size());

    pn1.kmer_prg = kg;
    pn1.kmer_prg.sort_topologically();
    EXPECT_EQ((uint)7, pn1.kmer_prg.nodes.size());
    EXPECT_EQ((uint)7, pn1.kmer_prg.sorted_nodes.size());

    // want this path to be the ref, so add first a couple of times so it is weighted
    kmp = {pn1.kmer_prg.sorted_nodes[0], pn1.kmer_prg.sorted_nodes[2], pn1.kmer_prg.sorted_nodes[5], pn1.kmer_prg.sorted_nodes[6]};
    pn1.add_path(kmp);
    pn1.add_path(kmp);

    EXPECT_EQ((uint)2, pn1.kmer_prg.sorted_nodes[0]->covg[0]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[1]->covg[0]);
    EXPECT_EQ((uint)2, pn1.kmer_prg.sorted_nodes[2]->covg[0]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[3]->covg[0]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[4]->covg[0]);
    EXPECT_EQ((uint)2, pn1.kmer_prg.sorted_nodes[5]->covg[0]);
    EXPECT_EQ((uint)2, pn1.kmer_prg.sorted_nodes[6]->covg[0]);
    EXPECT_EQ((uint)2, pn1.kmer_prg.sorted_nodes[0]->covg[1]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[1]->covg[1]);
    EXPECT_EQ((uint)2, pn1.kmer_prg.sorted_nodes[2]->covg[1]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[3]->covg[1]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[4]->covg[1]);
    EXPECT_EQ((uint)2, pn1.kmer_prg.sorted_nodes[5]->covg[1]);
    EXPECT_EQ((uint)2, pn1.kmer_prg.sorted_nodes[6]->covg[1]);

    // define 3 samples, one with ref path and 2 with alts
    SamplePtr ps1(make_shared<Sample>("sample1"));
    ps1->paths[3] = {kmp};
    pn1.samples.insert(ps1);
    pn1.add_path(kmp);

    kmp = {pn1.kmer_prg.sorted_nodes[0], pn1.kmer_prg.sorted_nodes[1], pn1.kmer_prg.sorted_nodes[4], pn1.kmer_prg.sorted_nodes[6]};
    SamplePtr ps2(make_shared<Sample>("sample2"));
    ps2->paths[3] = {kmp};
    pn1.samples.insert(ps2);
    pn1.add_path(kmp);

    kmp = {pn1.kmer_prg.sorted_nodes[0], pn1.kmer_prg.sorted_nodes[3], pn1.kmer_prg.sorted_nodes[6]};
    SamplePtr ps3(make_shared<Sample>("sample3"));
    ps3->paths[3] = {kmp};
    pn1.samples.insert(ps3);
    pn1.add_path(kmp);

    EXPECT_EQ((uint)5, pn1.kmer_prg.sorted_nodes[0]->covg[0]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[1]->covg[0]);
    EXPECT_EQ((uint)3, pn1.kmer_prg.sorted_nodes[2]->covg[0]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[3]->covg[0]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[4]->covg[0]);
    EXPECT_EQ((uint)3, pn1.kmer_prg.sorted_nodes[5]->covg[0]);
    EXPECT_EQ((uint)5, pn1.kmer_prg.sorted_nodes[6]->covg[0]);
    EXPECT_EQ((uint)5, pn1.kmer_prg.sorted_nodes[0]->covg[1]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[1]->covg[1]);
    EXPECT_EQ((uint)3, pn1.kmer_prg.sorted_nodes[2]->covg[1]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[3]->covg[1]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[4]->covg[1]);
    EXPECT_EQ((uint)3, pn1.kmer_prg.sorted_nodes[5]->covg[1]);
    EXPECT_EQ((uint)5, pn1.kmer_prg.sorted_nodes[6]->covg[1]);

    pn1.covg = 3;

    LocalPRG* l3;
    l3 = new LocalPRG(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T");

    // clear the kmergraph and vectors to check that the shared pointers have really been kept valid
    // within pg!
    kg.clear();
    kmp.clear();

    pn1.output_samples(l3,"../test/test_cases/updatevcf_pannode",0);

    delete l3;

    // now check is as expect
    // find date
    time_t t = time(0);
    char mbstr[10];
    strftime(mbstr, sizeof(mbstr), "%d/%m/%y", localtime(&t));
    string dat(mbstr);

    string vcffile = R"(##fileformat=VCFv4.3
##fileDate==)" + dat + R"(
##ALT=<ID=SNP,Description="SNP">
##ALT=<ID=PH_SNPs,Description="Phased SNPs">
##ALT=<ID=INDEL,Description="Insertion-deletion">
##ALT=<ID=COMPLEX,Description="Complex variant, collection of SNPs and indels">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of variant">
##ALT=<ID=SIMPLE,Description="Graph bubble is simple">
##ALT=<ID=NESTED,Description="Variation site was a nested feature in the graph">
##ALT=<ID=TOO_MANY_ALTS,Description="Variation site was a multinested feature with too many alts to include all in the VCF">
##INFO=<ID=GRAPHTYPE,Number=1,Type=String,Description="Type of graph feature">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample3	sample2	sample1
nested varsite	1	.	GT	G	.	.	SVTYPE=INDEL;GRAPHTYPE=NESTED	GT	1	.	0
nested varsite	2	.	T	C	.	.	SVTYPE=SNP;GRAPHTYPE=NESTED	GT	.	1	0
)";
    string vcffile2 = R"(##fileformat=VCFv4.3
##fileDate==)" + dat + R"(
##ALT=<ID=SNP,Description="SNP">
##ALT=<ID=PH_SNPs,Description="Phased SNPs">
##ALT=<ID=INDEL,Description="Insertion-deletion">
##ALT=<ID=COMPLEX,Description="Complex variant, collection of SNPs and indels">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of variant">
##ALT=<ID=SIMPLE,Description="Graph bubble is simple">
##ALT=<ID=NESTED,Description="Variation site was a nested feature in the graph">
##ALT=<ID=TOO_MANY_ALTS,Description="Variation site was a multinested feature with too many alts to include all in the VCF">
##INFO=<ID=GRAPHTYPE,Number=1,Type=String,Description="Type of graph feature">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample2	sample3	sample1
nested varsite	1	.	GT	G	.	.	SVTYPE=INDEL;GRAPHTYPE=NESTED	GT	.	1	0
nested varsite	2	.	T	C	.	.	SVTYPE=SNP;GRAPHTYPE=NESTED	GT	1	.	0
)";

    ifstream ifs("../test/test_cases/updatevcf_pannode.three.multisample.vcf");
    string content( (std::istreambuf_iterator<char>(ifs) ),(std::istreambuf_iterator<char>()) );
    //EXPECT_EQ(vcffile,content);
    //EXPECT_EQ(vcffile2,content);
    EXPECT_EQ((vcffile==content) or (vcffile2==content), true);

    string fafile1 = R"(>sample3
AG--T
>sample1
A-GTT
>sample2
A-GCT
)";
    string fafile2 = R"(>sample3
AG--T
>sample2
A-GCT
>sample1
A-GTT
)";
    string fafile3 = R"(>sample2
A-GCT
>sample3
AG--T
>sample1
A-GTT
)";
    string fafile4 = R"(>sample2
A-GCT
>sample1
A-GTT
>sample3
AG--T
)";
    string fafile5 = R"(>sample1
A-GTT
>sample2
A-GCT
>sample3
AG--T
)";
    string fafile6 = R"(>sample1
A-GTT
>sample3
AG--T
>sample2
A-GCT
)";

    ifstream ifs2("../test/test_cases/updatevcf_pannode.three.multisample.fa");
    string content2( (std::istreambuf_iterator<char>(ifs2) ),(std::istreambuf_iterator<char>()) );
    EXPECT_EQ((fafile1==content2) or (fafile2==content2) or(fafile3==content2) or (fafile4==content2) or (fafile5==content2) or (fafile6==content2), true);
}

TEST(PangenomeNodeTest,output_samples2)
{
    Graph pg;
    LocalPRG lp(5, "five", "A 5 G 7 C 8 T 7  6 G 5 T");

    vector<KmerNodePtr> kmp;

    KmerGraph kg;
    deque<Interval> d = {Interval(0,0)};
    Path p;
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(8, 9)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(8, 9), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(4,5), Interval(12, 13)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(4,5), Interval(12, 13), Interval(16,16), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(0,1), Interval(19,20), Interval(23,24)};
    p.initialize(d);
    kg.add_node(p);
    d = {Interval(24,24)};
    p.initialize(d);
    kg.add_node(p);
    kg.sort_topologically();
    kg.add_edge(kg.sorted_nodes[0],kg.sorted_nodes[1]);
    kg.add_edge(kg.sorted_nodes[1],kg.sorted_nodes[4]);
    kg.add_edge(kg.sorted_nodes[0],kg.sorted_nodes[2]);
    kg.add_edge(kg.sorted_nodes[2],kg.sorted_nodes[5]);
    kg.add_edge(kg.sorted_nodes[0],kg.sorted_nodes[3]);
    kg.add_edge(kg.sorted_nodes[3],kg.sorted_nodes[6]);
    kg.add_edge(kg.sorted_nodes[4],kg.sorted_nodes[6]);
    kg.add_edge(kg.sorted_nodes[5],kg.sorted_nodes[6]);

    EXPECT_EQ((uint)7, kg.nodes.size());

    lp.kmer_prg = kg;
    lp.kmer_prg.sort_topologically();
    EXPECT_EQ((uint)7, lp.kmer_prg.nodes.size());
    EXPECT_EQ((uint)7, lp.kmer_prg.sorted_nodes.size());

    kmp = {lp.kmer_prg.sorted_nodes[0], lp.kmer_prg.sorted_nodes[2], lp.kmer_prg.sorted_nodes[5], lp.kmer_prg.sorted_nodes[6]};
    //pg.add_node(5, "five", "sample0", kmp, &lp);
    pg.add_node(5, "five", "sample1", kmp, &lp);

    kmp = {lp.kmer_prg.sorted_nodes[0], lp.kmer_prg.sorted_nodes[1], lp.kmer_prg.sorted_nodes[4], lp.kmer_prg.sorted_nodes[6]};
    pg.add_node(5, "five", "sample2", kmp, &lp);

    kmp = {lp.kmer_prg.sorted_nodes[0], lp.kmer_prg.sorted_nodes[3], lp.kmer_prg.sorted_nodes[6]};
    pg.add_node(5, "five", "sample3", kmp, &lp);

    // clear the kmergraph and vectors to check that the shared pointers have really been kept valid
    // within pg!
    kg.clear();
    kmp.clear();

    pg.nodes[5]->output_samples(&lp,"../test/test_cases/updatevcf_pannode2",0);

    // now check is as expect
    // find date
    time_t t = time(0);
    char mbstr[10];
    strftime(mbstr, sizeof(mbstr), "%d/%m/%y", localtime(&t));
    string dat(mbstr);

    string vcffile = R"(##fileformat=VCFv4.3
##fileDate==)" + dat + R"(
##ALT=<ID=SNP,Description="SNP">
##ALT=<ID=PH_SNPs,Description="Phased SNPs">
##ALT=<ID=INDEL,Description="Insertion-deletion">
##ALT=<ID=COMPLEX,Description="Complex variant, collection of SNPs and indels">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of variant">
##ALT=<ID=SIMPLE,Description="Graph bubble is simple">
##ALT=<ID=NESTED,Description="Variation site was a nested feature in the graph">
##ALT=<ID=TOO_MANY_ALTS,Description="Variation site was a multinested feature with too many alts to include all in the VCF">
##INFO=<ID=GRAPHTYPE,Number=1,Type=String,Description="Type of graph feature">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample3	sample2	sample1
five	1	.	GC	G	.	.	SVTYPE=INDEL;GRAPHTYPE=NESTED	GT	1	0	.
five	2	.	C	T	.	.	SVTYPE=SNP;GRAPHTYPE=NESTED	GT	.	0	1
)";

    ifstream ifs("../test/test_cases/updatevcf_pannode2.five.multisample.vcf");
    string content( (std::istreambuf_iterator<char>(ifs) ),(std::istreambuf_iterator<char>()) );
    EXPECT_EQ(vcffile,content);
    //EXPECT_EQ(vcffile2,content);
    //EXPECT_EQ((vcffile==content) or (vcffile2==content), true);

    string fafile1 = R"(>sample3
AG--T
>sample1
A-GTT
>sample2
A-GCT
)";
    string fafile2 = R"(>sample3
AG--T
>sample2
A-GCT
>sample1
A-GTT
)";
    string fafile3 = R"(>sample2
A-GCT
>sample3
AG--T
>sample1
A-GTT
)";
    string fafile4 = R"(>sample2
A-GCT
>sample1
A-GTT
>sample3
AG--T
)";
    string fafile5 = R"(>sample1
A-GTT
>sample2
A-GCT
>sample3
AG--T
)";
    string fafile6 = R"(>sample1
A-GTT
>sample3
AG--T
>sample2
A-GCT
)";

    ifstream ifs2("../test/test_cases/updatevcf_pannode2.five.multisample.fa");
    string content2( (std::istreambuf_iterator<char>(ifs2) ),(std::istreambuf_iterator<char>()) );
    EXPECT_EQ((fafile1==content2) or (fafile2==content2) or(fafile3==content2) or (fafile4==content2) or (fafile5==content2) or (fafile6==content2), true);
}
*/

TEST(PangenomeNodeTest, equals) {
    Node pn1(3, 3, "3");
    Node pn2(2, 2, "2");
    Node pn3(2, 2, "2");

    EXPECT_EQ(pn1, pn1);
    EXPECT_EQ(pn2, pn2);
    EXPECT_EQ(pn3, pn3);
    EXPECT_EQ(pn2, pn3);
    EXPECT_EQ(pn3, pn2);
    EXPECT_EQ((pn1 == pn2), false);
    EXPECT_EQ((pn1 == pn3), false);
}

TEST(PangenomeNodeTest, nequals) {
    Node pn1(3, 3, "3");
    Node pn2(2, 2, "2");
    Node pn3(2, 2, "2");

    EXPECT_EQ((pn1 != pn2), true);
    EXPECT_EQ((pn2 != pn1), true);
    EXPECT_EQ((pn1 != pn1), false);
    EXPECT_EQ((pn2 != pn2), false);
    EXPECT_EQ((pn3 != pn3), false);
    EXPECT_EQ((pn2 != pn3), false);
}

TEST(PangenomeNodeTest, less) {
    Node pn1(3, 3, "3");
    Node pn2(2, 2, "2");
    Node pn3(2, 2, "2");

    EXPECT_EQ((pn1 < pn1), false);
    EXPECT_EQ((pn2 < pn2), false);
    EXPECT_EQ((pn3 < pn3), false);
    EXPECT_EQ((pn1 < pn3), false);
    EXPECT_EQ((pn1 < pn2), false);
    EXPECT_EQ((pn2 < pn1), true);
    EXPECT_EQ((pn3 < pn1), true);

}
