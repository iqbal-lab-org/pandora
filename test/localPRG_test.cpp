#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "localPRG.h"
#include "minimizer.h"
#include "interval.h"
#include "localgraph.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class LocalPRGTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }

  LocalPRG l0 = LocalPRG(0,"empty", "", 1, 3);
  LocalPRG l1 = LocalPRG(1,"simple", "AGCT", 1, 3);
  LocalPRG l2 = LocalPRG(2,"varsite", "A 5 GC 6 G 5 T", 1, 3);
  LocalPRG l3 = LocalPRG(3,"nested varsite", "A 5 G 7 C 8 T 7  6 G 5 T", 1, 3);

};

TEST_F(LocalPRGTest, isalphaEmptyString){
    bool a0 = l0.isalpha_string("");
    bool a1 = l1.isalpha_string("");
    bool a2 = l2.isalpha_string("");
    bool a3 = l3.isalpha_string("");
    EXPECT_EQ(a0, 1) << "isalpha_string thinks the empty string is not alphabetic for l0";
    EXPECT_EQ(a1, 1) << "isalpha_string thinks the empty string is not alphabetic for l1";
    EXPECT_EQ(a2, 1) << "isalpha_string thinks the empty string is not alphabetic for l2";
    EXPECT_EQ(a3, 1) << "isalpha_string thinks the empty string is not alphabetic for l3";
}

TEST_F(LocalPRGTest, isalphaSpaceString){
    bool a0 = l0.isalpha_string("AGCT T");
    bool a1 = l1.isalpha_string("AGCT T");
    bool a2 = l2.isalpha_string("AGCT T");
    bool a3 = l3.isalpha_string("AGCT T");         
    EXPECT_EQ(a0, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l0";
    EXPECT_EQ(a1, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l1";
    EXPECT_EQ(a2, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l2";
    EXPECT_EQ(a3, 0) << "isalpha_string thinks a string  containing a space is alphabetic for l3";
}

TEST_F(LocalPRGTest, isalphaNumberString){
    bool a0 = l0.isalpha_string("AGCT 8 T");
    bool a1 = l1.isalpha_string("AGCT 8 T");
    bool a2 = l2.isalpha_string("AGCT 8 T");
    bool a3 = l3.isalpha_string("AGCT 8 T");
    EXPECT_EQ(a0, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l0";
    EXPECT_EQ(a1, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l1";
    EXPECT_EQ(a2, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l2";
    EXPECT_EQ(a3, 0) << "isalpha_string thinks a string  containing a number is alphabetic for l3";
}

TEST_F(LocalPRGTest, splitBySiteNoSites){
    vector<Interval> v0, v1;
    v0.push_back(Interval(0,0));
    EXPECT_ITERABLE_EQ( vector< Interval >, v0, l0.splitBySite(Interval(0,0)));// << "Failed to split empty string with input Interval";
    v1.push_back(Interval(0,4));
    EXPECT_ITERABLE_EQ( vector< Interval >, v1, l1.splitBySite(Interval(0,4)));// << "Failed to split string with input Interval";
    v1.clear();
    v1.push_back(Interval(0,2));
    EXPECT_ITERABLE_EQ( vector< Interval >, v1, l1.splitBySite(Interval(0,2)));// << "Failed to split string with short input Interval";
    v1.clear();
    v1.push_back(Interval(1,3));
    EXPECT_ITERABLE_EQ( vector< Interval >, v1, l1.splitBySite(Interval(1,3)));// << "Failed to split string with middle input Interval";
}

TEST_F(LocalPRGTest, splitBySiteSite){
    vector<Interval> v2;
    v2.push_back(Interval(0,1));
    l2.next_site = 5;
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,1)));// << "Failed to split string in Interval" << Interval(0,1);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,2)));// << "Failed to split string in Interval" << Interval(0,2);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,3)));// << "Failed to split string in Interval" << Interval(0,3);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,4)));// << "Failed to split string in Interval" << Interval(0,4);
    v2.push_back(Interval(4,6));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,6)));// << "Failed to split string in Interval" << Interval(0,6);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,7)));// << "Failed to split string in Interval" << Interval(0,7);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,8)));// << "Failed to split string in Interval" << Interval(0,8);
    v2.push_back(Interval(9,10));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,10)));// << "Failed to split string in Interval" << Interval(0,10);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,11)));// << "Failed to split string in Interval" << Interval(0,11);
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,12)));// << "Failed to split string in Interval" << Interval(0,12);
    v2.push_back(Interval(13,14));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(0,14)));// << "Failed to split string in Interval" << Interval(0,14);
    v2.clear();
    v2.push_back(Interval(5,6));
    EXPECT_ITERABLE_EQ( vector< Interval >,v2, l2.splitBySite(Interval(5,8)));// << "Failed to split string in mid Interval" << Interval(5,8);
}

TEST_F(LocalPRGTest, splitBySiteNestedSite){
    vector<Interval> v3;
    v3.push_back(Interval(0,1));
    l3.next_site = 5;
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,1)));// << "Failed to split string in Interval" << Interval(0,1);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,2)));// << "Failed to split string in Interval" << Interval(0,2);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,3)));// << "Failed to split string in Interval" << Interval(0,3);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,4)));// << "Failed to split string in Interval" << Interval(0,4);
    v3.push_back(Interval(4,16));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,16)));// << "Failed to split string in Interval" << Interval(0,6);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,17)));// << "Failed to split string in Interval" << Interval(0,7);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,18)));// << "Failed to split string in Interval" << Interval(0,8);
    v3.push_back(Interval(19,20));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,20)));// << "Failed to split string in Interval" << Interval(0,10);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,21)));// << "Failed to split string in Interval" << Interval(0,11);
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,22)));// << "Failed to split string in Interval" << Interval(0,12);
    v3.push_back(Interval(23,24));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(0,24)));// << "Failed to split string in Interval" << Interval(0,14);
    l3.next_site = 7;    
    v3.clear();
    v3.push_back(Interval(4,5));
    v3.push_back(Interval(8,9));
    v3.push_back(Interval(12,13));
    v3.push_back(Interval(16,16));
    EXPECT_ITERABLE_EQ( vector< Interval >,v3, l3.splitBySite(Interval(4,16)));// << "Failed to split string in mid Interval" << Interval(5,8);
}

/*TEST_F(LocalPRGTest, buildGraph)
{
    LocalGraph lg0;
    lg0.add_node(0,"",Interval(0,0));
    EXPECT_ITERABLE_EQ( map<uint32_t, LocalNode*>, lg0, l0.nodes);

    LocalGraph lg1;
    lg1.add_node(0,"AGCT", Interval(0,4));
    EXPECT_ITERABLE_EQ( map<uint32_t, LocalNode*>, lg1, l1.nodes);

    LocalGraph lg2;
    lg2.add_node(0,"A", Interval(0,1));
    lg2.add_node(1,"GC", Interval(4,6));
    lg2.add_node(2,"G", Interval(7,8));
    lg2.add_node(3,"T", Interval(13,14));
    lg2.add_edge(0,1);
    lg2.add_edge(0,2);
    lg2.add_edge(1,3);
    lg2.add_edge(2,3);
    EXPECT_ITERABLE_EQ( map<uint32_t, LocalNode*>, lg2, l2.nodes);

    LocalGraph lg3;
    lg3.add_node(0,"A", Interval(0,1));
    lg3.add_node(1,"G", Interval(4,5));
    lg3.add_node(2,"C", Interval(8,9));
    lg3.add_node(3,"T", Interval(12,13));
    lg3.add_node(4,"", Interval(16,16));
    lg3.add_node(5,"G", Interval(19,20));
    lg3.add_node(6,"T", Interval(23,24));
    lg3.add_edge(0,1);
    lg3.add_edge(0,5);
    lg3.add_edge(1,2);
    lg3.add_edge(1,3);
    lg3.add_edge(2,4);
    lg3.add_edge(3,4);
    lg3.add_edge(4,6);
    lg3.add_edge(5,6);
    EXPECT_ITERABLE_EQ( map<uint32_t, LocalNode*>, lg3, l3.nodes);
}*/
//could test that no int is included in more than one node

/*TEST_F(LocalPRGTest,unpackLinearString){
    Seq s1 = Seq(0,"0", "AGCTAATGCGTT", 11, 3);
    Seq s2 = Seq(0,"0", "AGCTAATGCGTT", 10, 3);
    Seq s3 = Seq(0,"0", "AGCTAATGCGTT", 9, 3);
    Seq s4 = Seq(0,"0", "AGCTAGTGCGTT", 9, 3);
    uint32_t j = 0;
    EXPECT_EQ(s1.sketch.size(),j) << "Have " << s1.sketch.size() << " minimizer when string is too short";
    ++j;
    EXPECT_EQ(s2.sketch.size(),j) << "Have " << s2.sketch.size() << " minimizers when should have 1";
    EXPECT_EQ(s3.sketch.size(),j) << "Have " << s3.sketch.size() << " minimizers when should have 1";
    ++j;
    EXPECT_EQ(s4.sketch.size(),j) << "Have " << s4.sketch.size() << " minimizers when should have 2";
}

TEST_F(SeqTest,sketchIncludesEveryLetter){
    Seq s1 = Seq(0,"0", "AGCTAATGTGTT", 3, 3);
    Seq s2 = Seq(0,"0", "AGCTAATGTGTT", 2, 3);
    Seq s3 = Seq(0,"0", "AGCTAATGTGTT", 1, 3);
    Seq s4 = Seq(0,"0", "AGCTAATGTGAT", 3, 3);

    set<int> pos_inc;
    for(set<Minimizer*>::iterator it=s4.sketch.begin(); it != s4.sketch.end(); ++it)
    {
        for (std::deque<Interval>::iterator it2=((*it)->path).path.begin(); it2!=((*it)->path).path.end(); ++it2)
        {
            for (uint32_t j = it2->start; j<it2->end; ++j)
            {
                pos_inc.insert(j);
            }
        }
    } 
    set<int> expected = {0,1,2,3,4,5,6,7,8,9,10,11};
    EXPECT_EQ(pos_inc, expected) << "sketch misses a letter";
 
    uint32_t j = 10;
    EXPECT_EQ(s3.sketch.size(), j) << "sketch with w=1 has incorrect size " << s3.sketch.size();
    
    pos_inc.clear();
    for(set<Minimizer*>::iterator it=s2.sketch.begin(); it != s2.sketch.end(); ++it)
    {
        for (std::deque<Interval>::iterator it2=((*it)->path).path.begin(); it2!=((*it)->path).path.end(); ++it2)
        {
            for (uint32_t j = it2->start; j<it2->end; ++j)
            {
                pos_inc.insert(j);
            }
        }
    }
    EXPECT_EQ(pos_inc, expected) << "sketch for s2 includes/misses wrong letter";

    pos_inc.clear();
    for(set<Minimizer*>::iterator it=s1.sketch.begin(); it != s1.sketch.end(); ++it)
    {
        for (std::deque<Interval>::iterator it2=((*it)->path).path.begin(); it2!=((*it)->path).path.end(); ++it2)
        {
            for (uint32_t j = it2->start; j<it2->end; ++j)
            {
                pos_inc.insert(j);
            }
        }
    }
    expected = {0,1,2,3,4,5,6,7,8,9};
    EXPECT_EQ(pos_inc, expected) << "sketch for s1 includes/misses wrong letter";

}*/

//TEST_F(SeqTest,sketchCorrect){
//}
