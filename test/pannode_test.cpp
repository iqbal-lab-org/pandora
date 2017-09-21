#include "gtest/gtest.h"
#include "pannode.h"
#include "pansample.h"
#include "minihit.h"
#include "localPRG.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class PanNodeTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};

TEST_F(PanNodeTest,create){

    PanNode pn(4, 3, "3");
    uint32_t j=3;
    EXPECT_EQ(j, pn.node_id);
    EXPECT_EQ((uint)4, pn.prg_id);
    EXPECT_EQ("3", pn.name);
    EXPECT_EQ((uint)1, pn.covg);
    EXPECT_EQ((uint)0, pn.edges.size());
    EXPECT_EQ((uint)0, pn.reads.size());
    EXPECT_EQ((uint)0, pn.samples.size());
}

TEST_F(PanNodeTest,get_name){
    PanNode pn1(3,3,"3");
    PanNode pn2(2,2,"2");
    PanNode pn3(2,4,"2");

    EXPECT_EQ(pn1.get_name(), "3");
    EXPECT_EQ(pn2.get_name(), "2");
    EXPECT_EQ(pn3.get_name(), "2.4");
}

TEST_F(PanNodeTest,add_path){
    PanNode pn1(3,3,"3");
    vector<KmerNode*> kmp;
    pn1.add_path(kmp);

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
    EXPECT_EQ((uint)7, kg.nodes.size());
    
    pn1.kmer_prg = kg;
    EXPECT_EQ((uint)7, pn1.kmer_prg.nodes.size());
    pn1.kmer_prg.sort_topologically();
    EXPECT_EQ((uint)7, pn1.kmer_prg.sorted_nodes.size());
    kmp = {pn1.kmer_prg.sorted_nodes[0], pn1.kmer_prg.sorted_nodes[3], pn1.kmer_prg.sorted_nodes[4], pn1.kmer_prg.sorted_nodes[6]};
    pn1.add_path(kmp);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[0]->covg[0]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[1]->covg[0]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[2]->covg[0]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[3]->covg[0]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[4]->covg[0]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[5]->covg[0]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[6]->covg[0]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[0]->covg[1]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[1]->covg[1]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[2]->covg[1]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[3]->covg[1]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[4]->covg[1]);
    EXPECT_EQ((uint)0, pn1.kmer_prg.sorted_nodes[5]->covg[1]);
    EXPECT_EQ((uint)1, pn1.kmer_prg.sorted_nodes[6]->covg[1]);
}

TEST_F(PanNodeTest,output_samples_vcf)
{
    PanNode pn1(3,3,"three");
    vector<KmerNode*> kmp;

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
    PanSample *ps1;
    PanSample *ps2;
    PanSample *ps3;
    ps1 = new PanSample("sample1");
    ps1->paths[3] = {kmp};
    pn1.samples.insert(ps1);
    pn1.add_path(kmp);

    kmp = {pn1.kmer_prg.sorted_nodes[0], pn1.kmer_prg.sorted_nodes[1], pn1.kmer_prg.sorted_nodes[4], pn1.kmer_prg.sorted_nodes[6]};
    ps2 = new PanSample("sample2");
    ps2->paths[3] = {kmp};
    pn1.samples.insert(ps2);
    pn1.add_path(kmp);

    kmp = {pn1.kmer_prg.sorted_nodes[0], pn1.kmer_prg.sorted_nodes[3], pn1.kmer_prg.sorted_nodes[6]};
    ps3 = new PanSample("sample3");
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

    pn1.output_samples_vcf(l3,"../test/test_cases/updatevcf_pannode",0);

    delete ps1;
    delete ps2;
    delete ps3;
    delete l3;
}


TEST_F(PanNodeTest,equals){
    PanNode pn1(3,3,"3");
    PanNode pn2(2,2,"2");
    PanNode pn3(2,2,"2");

    EXPECT_EQ(pn1, pn1);
    EXPECT_EQ(pn2, pn2);
    EXPECT_EQ(pn3, pn3);
    EXPECT_EQ(pn2, pn3);
    EXPECT_EQ(pn3, pn2);
    EXPECT_EQ((pn1==pn2), false);
    EXPECT_EQ((pn1==pn3), false);
}

TEST_F(PanNodeTest,nequals){
    PanNode pn1(3,3,"3");
    PanNode pn2(2,2,"2");
    PanNode pn3(2,2,"2");

    EXPECT_EQ((pn1!=pn2), true);
    EXPECT_EQ((pn2!=pn1), true);
    EXPECT_EQ((pn1!=pn1), false);
    EXPECT_EQ((pn2!=pn2), false);
    EXPECT_EQ((pn3!=pn3), false);
    EXPECT_EQ((pn2!=pn3), false);
}

TEST_F(PanNodeTest,less){
    PanNode pn1(3,3,"3");
    PanNode pn2(2,2,"2");
    PanNode pn3(2,2,"2");

    EXPECT_EQ((pn1<pn1), false);
    EXPECT_EQ((pn2<pn2), false);
    EXPECT_EQ((pn3<pn3), false);
    EXPECT_EQ((pn1<pn3), false);
    EXPECT_EQ((pn1<pn2), false);
    EXPECT_EQ((pn2<pn1), true);
    EXPECT_EQ((pn3<pn1), true);

}
