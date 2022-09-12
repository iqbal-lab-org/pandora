#include "gtest/gtest.h"
#include "seq.h"
#include "interval.h"
#include <cstdint>
#include <iostream>

using namespace std;

TEST(SeqTest, create)
{
    Seq s1(0, "0", "AGCTAATGCGTT", 11, 3);
    EXPECT_EQ((uint)0, s1.id);
    EXPECT_EQ("0", s1.name);
    const std::vector<std::string> expected_seq { "AGCTAATGCGTT" };
    EXPECT_EQ(expected_seq, s1.seq);
}

TEST(SeqTest, initialize)
{
    Seq s1(0, "0", "AGCTAATGCGTT", 11, 3);
    s1.initialize(1, "new", "AGCTAATGCATA", 9, 3);
    EXPECT_EQ((uint)1, s1.id);
    EXPECT_EQ("new", s1.name);
    const std::vector<std::string> expected_seq { "AGCTAATGCATA" };
    EXPECT_EQ(expected_seq, s1.seq);
}

TEST(SeqTest, sketchShortReads)
{
    Seq s1(0, "0", "AGCTAATGCGTT", 11, 3);
    Seq s2(0, "0", "AGCTAATGCGTT", 10, 3);
    Seq s3(0, "0", "AGCTAATGCGTT", 9, 3);
    Seq s4(0, "0", "AGCTAATGCATA", 9, 3);
    uint32_t j = 0;
    EXPECT_EQ(s1.sketch.size(), j)
        << "Have " << s1.sketch.size() << " minimizer when string is too short";
    ++j;
    EXPECT_EQ(s2.sketch.size(), j)
        << "Have " << s2.sketch.size() << " minimizers when should have 1";
    ++j;
    EXPECT_EQ(s3.sketch.size(), j)
        << "Have " << s3.sketch.size() << " minimizers when should have 2";
    j = 1;
    EXPECT_EQ(s4.sketch.size(), j)
        << "Have " << s4.sketch.size() << " minimizers when should have 1";
}

TEST(SeqTest, sketchIncludesEveryLetter)
{
    Seq s1(0, "0", "AGCTAATGTGTT", 3, 3);
    Seq s2(0, "0", "AGCTAATGTGTT", 2, 3);
    Seq s3(0, "0", "AGCTAATGTGTT", 1, 3);
    Seq s4(0, "0", "AGCTAATGTGAT", 3, 3);

    set<int> pos_inc;
    for (auto it = s4.sketch.begin(); it != s4.sketch.end(); ++it) {
        for (uint32_t j = (*it).pos_of_kmer_in_read.start;
             j < (*it).pos_of_kmer_in_read.get_end(); ++j) {
            pos_inc.insert(j);
        }
    }
    for (uint32_t i = 3 - 1; i != 12 - 3 + 1;
         ++i) // first or last w-1 may not be included
    {
        EXPECT_EQ((pos_inc.find(i) != pos_inc.end()), true);
    }

    uint32_t j = 10;
    EXPECT_EQ(s3.sketch.size(), j)
        << "sketch with w=1 has incorrect size " << s3.sketch.size();

    pos_inc.clear();
    for (auto it = s2.sketch.begin(); it != s2.sketch.end(); ++it) {
        for (uint32_t j = (*it).pos_of_kmer_in_read.start;
             j < (*it).pos_of_kmer_in_read.get_end(); ++j) {
            pos_inc.insert(j);
        }
    }
    for (uint32_t i = 2 - 1; i != 12 - 2 + 1;
         ++i) // first or last w-1 may not be included
    {
        EXPECT_EQ((pos_inc.find(i) != pos_inc.end()), true);
    }

    pos_inc.clear();
    for (auto it = s1.sketch.begin(); it != s1.sketch.end(); ++it) {
        for (uint32_t j = (*it).pos_of_kmer_in_read.start;
             j < (*it).pos_of_kmer_in_read.get_end(); ++j) {
            pos_inc.insert(j);
        }
    }
    for (uint32_t i = 3 - 1; i != 12 - 3 + 1;
         ++i) // first or last w-1 may not be included
    {
        EXPECT_EQ((pos_inc.find(i) != pos_inc.end()), true);
    }
}

TEST(SeqTest, sketchSkipsAmbiguousBaseAtStart)
{
    const string seq = "NGCTAATGTGTT";
    const auto w { 1 };
    const auto k { 3 };
    Seq s1(0, "0", seq, w, k);

    set<int> pos_exclude { 0 };
    set<int> pos_include {};
    for (auto it = s1.sketch.begin(); it != s1.sketch.end(); ++it) {
        for (uint32_t j = (*it).pos_of_kmer_in_read.start;
             j < (*it).pos_of_kmer_in_read.get_end(); ++j) {
            EXPECT_TRUE(pos_exclude.find(j) == pos_exclude.end()) << (*it);
            pos_include.insert(j);
        }
    }
    EXPECT_EQ(pos_include.size(), seq.length() - 1);
}

TEST(SeqTest, lengthNoAmbiguous)
{
    const std::string s { "AGCTAATGCGTT" };
    const Seq seq(0, "0", s, 3, 3);

    EXPECT_EQ(seq.length(), s.length());
}

TEST(SeqTest, lengthOneAmbiguous)
{
    const std::string s { "AGCTAATGNGTT" };
    const Seq seq(0, "0", s, 3, 3);

    EXPECT_EQ(seq.length(), s.length() - 1);
}

TEST(SeqTest, lengthTwoAmbiguous)
{
    const std::string s { "AWGCTAATGNGTT" };
    const Seq seq(0, "0", s, 3, 3);

    EXPECT_EQ(seq.length(), s.length() - 2);
}

