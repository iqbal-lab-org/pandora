#include <stdint.h> 
#include <cstring>
#include <iostream>
#include "gtest/gtest.h"
#include "pangraph.h"
#include "estimate_parameters.h"

using namespace std;

class EstimateParametersTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(EstimateParametersTest, find_mean_covg){
    //NB this finds the position in vector at which max of the second peak occurs
    std::vector<uint> v1 = {30, 24, 12, 3, 6, 2, 14, 15, 16, 18, 40, 26, 35, 14};
    EXPECT_EQ(uint(10), find_mean_covg(v1));
    // not thrown by a single weird one not near peak
    std::vector<uint> v2 = {30, 24, 12, 3, 70, 2, 14, 15, 16, 18, 40, 26, 35, 14};
    EXPECT_EQ(uint(10), find_mean_covg(v2));
    // doesn't matter if second peak much lower
    std::vector<uint> v3 = {30, 24, 12, 3, 6, 2, 14, 15, 16, 18, 14, 8, 9, 1};
    EXPECT_EQ(uint(9), find_mean_covg(v3));
    // do need an increase three times
    std::vector<uint> v4 = {30, 24, 12, 3, 6, 2, 11, 10, 9, 8, 4, 3, 2, 1};
    EXPECT_EQ(uint(0), find_mean_covg(v4)); 
}

TEST_F(EstimateParametersTest, find_prob_thresh)
{
    //NB this finds the position in vector at which min occurs between 2 peaks
    std::vector<uint> v1 = {30, 24, 18, 16, 12, 3, 6, 2, 1, 15, 16, 18, 12, 26, 35, 40};
    EXPECT_EQ(8-200, find_prob_thresh(v1));
    // not thrown by low values outside of valley
    std::vector<uint> v2 = {1, 30, 24, 12, 3, 6, 2, 0, 15, 16, 18, 12, 26, 35, 40, 0};
    EXPECT_EQ(7-200, find_prob_thresh(v2));
}
