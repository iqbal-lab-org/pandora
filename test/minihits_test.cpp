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

using namespace std;

TEST(MinimizerHitsTest,add_hit){
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

TEST(MinimizerHitsTest, pComp) {
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

TEST(MinimizerHitsTest, pComp_path) {
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

TEST(MinimizerHitsTest, clusterComp){
    set<set<MinimizerHitPtr, pComp>,clusterComp> clusters_of_hits;
    set<MinimizerHitPtr, pComp> current_cluster;

    vector<MinimizerHitPtr> expected1, expected2;

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
    current_cluster.insert(mh);
    expected1.push_back(mh);

    mh = make_shared<MinimizerHit>(2, m, mr);
    current_cluster.insert(mh);
    expected1.push_back(mh);
    clusters_of_hits.insert(current_cluster);

    current_cluster.clear();
    delete m;
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);
    mh = make_shared<MinimizerHit>(1, m, mr);
    current_cluster.insert(mh);
    expected2.push_back(mh);

    d = {Interval(6,10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mh = make_shared<MinimizerHit>(1, m, mr);
    current_cluster.insert(mh);
    expected2.push_back(mh);

    d = {Interval(6,10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0,0);
    mh = make_shared<MinimizerHit>(1, m, mr);
    current_cluster.insert(mh);
    expected2.push_back(mh);

    clusters_of_hits.insert(current_cluster);

    // have inserted 2 clusters
    uint32_t j(2);
    EXPECT_EQ(j, clusters_of_hits.size());
    // expect the cluster of 3 added second to be first
    j=3;
    EXPECT_EQ(j, clusters_of_hits.begin()->size());
    
    // check that this is indeed the cluster we think
    // note that can't sort the set of hits by pcomp without changing the defintion of cluster comparison function
    for (set<MinimizerHitPtr, pComp>::iterator it=clusters_of_hits.begin()->begin(); it!=clusters_of_hits.begin()->end(); ++it)
    {
        EXPECT_EQ(((*expected2[0] == **it) or (*expected2[1] == **it) or (*expected2[2] == **it)), true);
    }

    delete m;
    delete mr;
}
