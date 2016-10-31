#include "gtest/gtest.h"
#include "minirecord.h"
#include "path.h"
#include "index.h"
#include "interval.h"
#include "inthash.h"
#include <vector>
#include <stdint.h>
#include <iostream>

using namespace std;

class IndexTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(IndexTest,addRecord){
    Index idx;
    deque<Interval> d = {Interval(3,5), Interval(9,12)};
    Path p;
    p.initialize(d);
    uint64_t kh = kmerhash("ACGTA",5);
    idx.add_record(kh, 1, p);
    uint32_t j=1;
    EXPECT_EQ(j, idx.minhash.size());

    // add again - should stay same size
    idx.add_record(kh, 1, p);
    EXPECT_EQ(j, idx.minhash.size());
    EXPECT_EQ(j, idx.minhash[kh].size());

    // add a new record with different key
    idx.add_record(kmerhash("ACTGA",5), 2, p);
    j=2;
    EXPECT_EQ(j, idx.minhash.size());

    // and a new record which is different but has same key
    idx.add_record(kh, 4, p);
    EXPECT_EQ(j, idx.minhash.size());
    EXPECT_EQ(j, idx.minhash[kh].size());
}

TEST_F(IndexTest, clear){
    Index idx;
    deque<Interval> d = {Interval(3,5), Interval(9,12)};
    Path p;
    p.initialize(d);
    idx.add_record(kmerhash("ACGTA",5), 1, p);
    idx.add_record(kmerhash("ACTGA",5), 2, p);
    idx.add_record(kmerhash("ACGTA",5), 4, p);
    idx.clear();
    uint32_t j = 0;
    EXPECT_EQ(j, idx.minhash.size());
}
