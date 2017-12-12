#include "gtest/gtest.h"
#include "minihits.h"
#include "minihit.h"
#include "minimizer.h"
#include "minirecord.h"
#include "interval.h"
#include "path.h"
#include "inthash.h"
#include <stdint.h>
#include <iostream>
#include <algorithm>
#include <utility>

using namespace std;

class MinimizerHitsTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(MinimizerHitsTest,add_hit){
    // tests both add_hit and that sort doesn't break. Doesn't test resut of sort
    MinimizerHits mhits;
    KmerHash hash;
    Minimizer* m;
    pair<uint64_t,uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = new Minimizer(min(kh.first,kh.second), 1,6,0);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MiniRecord* mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr);
    EXPECT_EQ((uint)1, mhits.uhits.size());   
    //mhits.add_hit(1, m, mr);
    //EXPECT_EQ((uint)1, mhits.uhits.size());
    mhits.add_hit(2, m, mr);
    EXPECT_EQ((uint)2, mhits.uhits.size());

    delete m;
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);
    mhits.add_hit(1, m, mr);
    EXPECT_EQ((uint)3, mhits.uhits.size());

    d = {Interval(6,10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr);
    EXPECT_EQ((uint)4, mhits.uhits.size());

    d = {Interval(6,10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr);
    EXPECT_EQ((uint)5, mhits.uhits.size());

    uint32_t j(5);
    mhits.sort();
    EXPECT_EQ(j, mhits.hits.size());

    delete m;
    delete mr;
}

TEST_F(MinimizerHitsTest, pComp) {
    MinimizerHits mhits;
    vector<MinimizerHit> expected;
    KmerHash hash;
    Minimizer* m;
    pair<uint64_t,uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = new Minimizer(min(kh.first,kh.second), 1,6,0);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MiniRecord* mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr);
    expected.push_back(MinimizerHit(1, m, mr));
    mhits.add_hit(0, m, mr);
    expected.push_back(MinimizerHit(0, m, mr));

    delete m;
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);

    d = {Interval(6,10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr);
    expected.push_back(MinimizerHit(1, m, mr));

    d = {Interval(6,10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr);
    expected.push_back(MinimizerHit(1, m, mr));

    mhits.sort();
    uint32_t j(1);
    for (set<MinimizerHitPtr, pComp>::iterator it=mhits.hits.begin(); it!=--mhits.hits.end(); ++it)
    {
        EXPECT_EQ(expected[j], **it);
        j++;
    }
    EXPECT_EQ(expected[0], **(--mhits.hits.end()));

    delete m;
    delete mr;
}

TEST_F(MinimizerHitsTest, pComp_path) {
    set<MinimizerHitPtr, pComp_path> mhitspath;
    MinimizerHits mhits;
    deque<MinimizerHit> expected;
    KmerHash hash;
    Minimizer* m;
    pair<uint64_t,uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = new Minimizer(min(kh.first,kh.second), 1,6,0);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MiniRecord* mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(0, m, mr);
    expected.push_back(MinimizerHit(0, m, mr));
    mhits.add_hit(1, m, mr);
    expected.push_back(MinimizerHit(1, m, mr));

    delete m;
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);
    mhits.add_hit(2, m, mr);
    expected.push_back(MinimizerHit(2, m, mr));

    d = {Interval(6,10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr);
    expected.push_front(MinimizerHit(1, m, mr));

    d = {Interval(6,10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mhits.add_hit(1, m, mr);
    expected.push_front(MinimizerHit(1, m, mr));

    mhits.sort();
    for (set<MinimizerHitPtr, pComp>::iterator it=mhits.hits.begin(); it!=--mhits.hits.end(); ++it)
    {
	mhitspath.insert(*it);
    }
    uint32_t j(0);
    for (set<MinimizerHitPtr, pComp_path>::iterator it=mhitspath.begin(); it!=mhitspath.end(); ++it)
    {
        EXPECT_EQ(expected[j], **it);
        j++;
    }

    delete m;
    delete mr;
}

TEST_F(MinimizerHitsTest, clusterComp){
    set<MinimizerHitPtr, pComp> hits;
    set<pair<set<MinimizerHitPtr, pComp>::iterator,
             set<MinimizerHitPtr, pComp>::iterator>,clusterComp> clusters_of_hits;

    vector<vector<MinimizerHitPtr>> expected(6);

    KmerHash hash;
    Minimizer* m;
    pair<uint64_t,uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = new Minimizer(min(kh.first,kh.second), 1,6,0);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MiniRecord* mr;
    mr = new MiniRecord(0,p,0,0);
    MinimizerHitPtr mh (make_shared<MinimizerHit>(1, m, mr));
    hits.insert(mh);
    expected[4].push_back(mh);

    mh = make_shared<MinimizerHit>(2, m, mr);
    hits.insert(mh);
    expected[5].push_back(mh);

    delete m;
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);
    mh = make_shared<MinimizerHit>(1, m, mr);
    hits.insert(mh);
    expected[0].push_back(mh);
    expected[1].push_back(mh);


    d = {Interval(6,10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mh = make_shared<MinimizerHit>(1, m, mr);
    hits.insert(mh);
    expected[1].push_back(mh);
    expected[2].push_back(mh);

    d = {Interval(6,10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mh = make_shared<MinimizerHit>(1, m, mr);
    hits.insert(mh);
    expected[1].push_back(mh);
    expected[3].push_back(mh);

    clusters_of_hits.insert(make_pair(next(hits.begin(), 3),next(hits.begin(), 3))); // hit 1
    clusters_of_hits.insert(make_pair(next(hits.begin(), 4),next(hits.begin(), 4))); // hit 2
    clusters_of_hits.insert(make_pair(hits.begin(), next(hits.begin(), 2))); // hits 3 4
    clusters_of_hits.insert(make_pair(next(hits.begin(), 2),next(hits.begin(), 2))); // hit 5
    clusters_of_hits.insert(make_pair(hits.begin(), next(hits.begin(), 2))); // hits 3 4 5
    clusters_of_hits.insert(make_pair(next(hits.begin(), 1),next(hits.begin(), 1))); // hit 4

    // have inserted 6 clusters
    uint32_t j(6);
    EXPECT_EQ(j, clusters_of_hits.size());

    // check clusters are in the order we expect
    j = 0;
    for (auto it : clusters_of_hits)
    {
        EXPECT_EQ((uint)distance(it.first, it.second), expected[j].size());
        if((uint)distance(it.first, it.second)==expected[j].size())
        {
            uint k=0;
            for (auto r=it.first; r!=it.second; ++r)
            {
                EXPECT_EQ(*r, expected[j][k]);
                k+=1;
            }
        }
        j+=1;
    }

    delete m;
    delete mr;
}
