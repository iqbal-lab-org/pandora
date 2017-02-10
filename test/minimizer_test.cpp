//#include <limits.h>
#include "gtest/gtest.h"
#include "minimizer.h"
#include "path.h"
#include "interval.h"
#include "inthash.h"
#include <set>
#include <vector>
#include <stdint.h> 
#include <iostream>

using std::set;
using namespace std;

struct Interval;
class Path;
struct Minimizer;

class MinimizerTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(MinimizerTest,create){
    KmerHash hash;
    pair<uint64_t,uint64_t> kh = hash.kmerhash("ACGTA", 5);
    Minimizer m1(kh.first, 0,5, 0);
    kh = hash.kmerhash("ACGTG", 5);
    Minimizer m2(kh.first, 1,6, 0);
    kh = hash.kmerhash("ACGTA", 5);
    Minimizer m3(kh.first, 5,10, 0);

    EXPECT_EQ(m1.kmer, kh.first);
    EXPECT_EQ(m3.kmer, kh.first);
    kh = hash.kmerhash("ACGTG", 5);
    EXPECT_EQ(m2.kmer, kh.first);

    uint32_t j = 0;
    EXPECT_EQ(m1.pos.start, j);
    j = 1;
    EXPECT_EQ(m2.pos.start, j);
    j = 5;
    EXPECT_EQ(m3.pos.start, j);

    EXPECT_EQ(m1.pos.end, j);
    j = 6;
    EXPECT_EQ(m2.pos.end, j);
    j = 10;
    EXPECT_EQ(m3.pos.end, j);

    EXPECT_DEATH(Minimizer(kh.first, 0,2,0),""); // interval too short to be valid
    //EXPECT_DEATH(Minimizer(kh.first, 0,8,0),""); // interval too long to be valid
    EXPECT_DEATH(Minimizer(kh.first, 2,0,0),""); // doesn't generate an interval as 2>0
}

TEST_F(MinimizerTest,comparisonCheck){
    KmerHash hash;
    pair<uint64_t,uint64_t> kh1 = hash.kmerhash("AGGTG", 5);
    Minimizer m1(kh1.first, 0,5,0);
    pair<uint64_t,uint64_t> kh2 = hash.kmerhash("ACGTA", 5);
    Minimizer m2(kh2.first, 1,6,0);
    Minimizer m3(kh1.first, 5,10,0); 
    Minimizer m4(kh2.first, 0,5,0);
    pair<uint64_t,uint64_t> kh3 = hash.kmerhash("ACGTG", 5);
    Minimizer m5(kh3.first, 0,5,0);

    //cout << kh1 << " " << kh2 << " " << kh3 << endl;
    set<Minimizer> s;
    s.insert(m1);
    s.insert(m2);
    s.insert(m3);
    s.insert(m4);
    s.insert(m5);

    uint32_t j = 5;
    EXPECT_EQ(s.size(),j) << "size of set of minimizers " << s.size() << " is not equal to 5.";

    // note this is a bad test as need to know the order of hash.kmerhash values to set this up
    vector<Minimizer> v = {m1, m3, m4, m2, m5};
    int i = 0;
    for (std::set<Minimizer>::iterator it=s.begin(); it!=s.end(); ++it)
    {
	EXPECT_EQ(it->kmer, v[i].kmer) << "kmers do not agree: " << it->kmer << ", " << v[i].kmer;
	EXPECT_EQ(it->pos.start, v[i].pos.start) << "start positions do not agree: " << it->pos.start << ", " << v[i].pos.start;
	EXPECT_EQ(it->pos.end, v[i].pos.end) << "end positions do not agree: " << it->pos.end << ", " << v[i].pos.end;
    	++i;
    }
}
