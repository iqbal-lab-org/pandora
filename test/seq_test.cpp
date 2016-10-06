#include "gtest/gtest.h"
#include "seq.h"
#include "minimizer.h"
#include "interval.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class SeqTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(SeqTest,sketchShortReads){
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

}

//TEST_F(SeqTest,sketchCorrect){
//}
