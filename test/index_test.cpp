#include "gtest/gtest.h"
#include "minirecord.h"
#include "path.h"
#include "index.h"
#include "interval.h"
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
    Index idx = Index();
    deque<Interval> d = {Interval(3,5), Interval(9,12)};
    Path p = Path();
    p.initialize(d);
    idx.add_record("hello", 1, p);
    EXPECT_EQ(1, idx.minhash.size());

    // add again - should stay same size
    idx.add_record("hello", 1, p);
    EXPECT_EQ(1, idx.minhash.size());

    // add a new record with different key
    idx.add_record("henno", 2, p);
    EXPECT_EQ(2, idx.minhash.size());

    // and a new record which is different but has same key
    idx.add_record("hello", 4, p);
    EXPECT_EQ(3, idx.minhash.size());
}
