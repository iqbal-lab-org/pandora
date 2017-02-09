#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "maxpath.h"
#include "localnode.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class MaxPathTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(MaxPathTest,create)
{
    LocalNode *ln1, *ln2, *ln3, *ln4;
    ln1 = new LocalNode("ACGTA", Interval(0,5), 0);
    ln2 = new LocalNode("AGCTA", Interval(0,5), 0);
    ln3 = new LocalNode("ACGTA", Interval(0,4), 0);
    ln4 = new LocalNode("ACGTA", Interval(0,5), 1);

    vector<LocalNode*> v = {ln1, ln2, ln3, ln4};
    vector<int> y = {0,1,1};
    
    MaxPath mp(v, y, 0);

    EXPECT_ITERABLE_EQ(vector<int>, mp.kmers_on_path, y);
    EXPECT_EQ((uint)0, mp.num_equivalent_paths);
    EXPECT_EQ(v.size(), mp.npath.size());
    for (uint32_t i = 0; i!= v.size(); ++i)
    {
	EXPECT_EQ(*(v[i]), *(mp.npath[i]));
    }
}

TEST_F(MaxPathTest,extend)
{
    LocalNode *ln1, *ln2, *ln3, *ln4;
    ln1 = new LocalNode("ACGTA", Interval(0,5), 0);
    ln2 = new LocalNode("AGCTA", Interval(0,5), 0);
    ln3 = new LocalNode("ACGTA", Interval(0,4), 0);
    ln4 = new LocalNode("ACGTA", Interval(0,5), 1);

    vector<LocalNode*> v = {ln1, ln2, ln3};
    vector<int> y = {0,1,1};
    MaxPath mp1(v, y, 0);

    v = {ln2, ln3, ln4};
    y = {1,1,0};
    MaxPath mp2(v, y, 0);

    mp1.extend(mp2);
    v = {ln1, ln2, ln3, ln4};
    y = {0,1,0};
    MaxPath mp3(v, y, 0);
    
    EXPECT_ITERABLE_EQ(vector<int>, mp1.kmers_on_path, mp3.kmers_on_path);
    EXPECT_EQ(mp1.npath.size(), mp3.npath.size());
    for (uint32_t i = 0; i!= mp1.npath.size(); ++i)
    {
        EXPECT_EQ(*(mp1.npath[i]), *(mp3.npath[i]));
    }
}

/*TEST_F(MaxPathTest,getprob)
{
}

TEST_F(MaxPathTest,getmeanprob)
{
}*/
