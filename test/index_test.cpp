#include "gtest/gtest.h"
#include "minirecord.h"
#include "prg/path.h"
#include "index.h"
#include "interval.h"
#include "inthash.h"
#include "utils.h"
#include <vector>
#include <stdint.h>
#include <iostream>
#include <algorithm>

using namespace std;

const std::string TEST_CASE_DIR = "../../test/test_cases/";

TEST(IndexTest, add_record)
{
    Index idx;
    pandora::KmerHash hash;
    deque<Interval> d = { Interval(3, 5), Interval(9, 12) };
    prg::Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    idx.add_record(min(kh.first, kh.second), 1, p, 0, 0);
    uint32_t j = 1;
    EXPECT_EQ(j, idx.minhash.size());

    // add again - should stay same size
    idx.add_record(min(kh.first, kh.second), 1, p, 0, 0);
    EXPECT_EQ(j, idx.minhash.size());
    EXPECT_EQ(j, idx.minhash[min(kh.first, kh.second)]->size());

    // add a new record with different key
    pair<uint64_t, uint64_t> kh2 = hash.kmerhash("ACTGA", 5);
    idx.add_record(min(kh2.first, kh2.second), 2, p, 0, 0);
    j = 2;
    EXPECT_EQ(j, idx.minhash.size());

    // and a new record which is different but has same key
    idx.add_record(min(kh.first, kh.second), 4, p, 0, 0);
    EXPECT_EQ(j, idx.minhash.size());
    EXPECT_EQ(j, idx.minhash[min(kh.first, kh.second)]->size());
}

TEST(IndexTest, clear)
{
    Index idx;
    pandora::KmerHash hash;
    deque<Interval> d = { Interval(3, 5), Interval(9, 12) };
    prg::Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    idx.add_record(min(kh.first, kh.second), 1, p, 0, 0);
    kh = hash.kmerhash("ACTGA", 5);
    idx.add_record(min(kh.first, kh.second), 2, p, 0, 0);
    kh = hash.kmerhash("ACGTA", 5);
    idx.add_record(min(kh.first, kh.second), 4, p, 0, 0);
    idx.clear();
    uint32_t j = 0;
    EXPECT_EQ(j, idx.minhash.size());
}

TEST(IndexTest, save)
{
    Index idx;
    pandora::KmerHash hash;
    deque<Interval> d = { Interval(3, 5), Interval(9, 12) };
    prg::Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    idx.add_record(min(kh.first, kh.second), 1, p, 0, 0);
    kh = hash.kmerhash("ACTGA", 5);
    idx.add_record(min(kh.first, kh.second), 2, p, 0, 0);
    kh = hash.kmerhash("ACGTA", 5);
    idx.add_record(min(kh.first, kh.second), 4, p, 0, 0);
    idx.save("indextext", 1, 5);
    ASSERT_TRUE(fopen("indextext.k5.w1.idx", "r") != NULL);
}

TEST(IndexTest, load)
{
    Index idx1, idx2;
    pandora::KmerHash hash;
    deque<Interval> d = { Interval(3, 5), Interval(9, 12) };
    prg::Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh1 = hash.kmerhash("ACGTA", 5);
    idx1.add_record(min(kh1.first, kh1.second), 1, p, 0, 0);
    pair<uint64_t, uint64_t> kh2 = hash.kmerhash("ACTGA", 5);
    idx1.add_record(min(kh2.first, kh2.second), 2, p, 0, 0);
    idx1.add_record(min(kh1.first, kh1.second), 4, p, 0, 0);

    idx2.load("indextext", 1, 5);
    EXPECT_EQ(idx1.minhash.size(), idx2.minhash.size());
    EXPECT_EQ(idx1.minhash[min(kh1.first, kh1.second)]->size(),
        idx2.minhash[min(kh1.first, kh1.second)]->size());
    EXPECT_EQ(idx1.minhash[min(kh2.first, kh2.second)]->size(),
        idx2.minhash[min(kh2.first, kh2.second)]->size());
    EXPECT_EQ(idx1.minhash[min(kh1.first, kh1.second)]->at(0),
        idx2.minhash[min(kh1.first, kh1.second)]->at(0));
    EXPECT_EQ(idx1.minhash[min(kh1.first, kh1.second)]->at(1),
        idx2.minhash[min(kh1.first, kh1.second)]->at(1));
    EXPECT_EQ(idx1.minhash[min(kh2.first, kh2.second)]->at(0),
        idx2.minhash[min(kh2.first, kh2.second)]->at(0));
}

TEST(IndexTest, equals)
{
    Index idx1, idx2;
    pandora::KmerHash hash;
    deque<Interval> d = { Interval(3, 5), Interval(9, 12) };
    prg::Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh1 = hash.kmerhash("ACGTA", 5);
    idx1.add_record(min(kh1.first, kh1.second), 1, p, 0, 0);
    pair<uint64_t, uint64_t> kh2 = hash.kmerhash("ACTGA", 5);
    idx1.add_record(min(kh2.first, kh2.second), 2, p, 0, 0);
    idx1.add_record(min(kh1.first, kh1.second), 4, p, 0, 0);

    idx2.load("indextext", 1, 5);
    EXPECT_EQ(idx1, idx2);
    EXPECT_EQ(idx2, idx1);
}

TEST(IndexTest, equals_fails)
{
    Index idx1, idx2;
    pandora::KmerHash hash;
    deque<Interval> d = { Interval(3, 5), Interval(9, 12) };
    prg::Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh1 = hash.kmerhash("ACGTA", 5);
    // idx1.add_record(min(kh1.first, kh1.second), 1, p, 0, 0);
    pair<uint64_t, uint64_t> kh2 = hash.kmerhash("ACTGA", 5);
    idx1.add_record(min(kh2.first, kh2.second), 2, p, 0, 0);
    // idx1.add_record(min(kh1.first, kh1.second), 4, p, 0, 0);

    idx2.load("indextext", 1, 5);
    EXPECT_NE(idx1, idx2);
    EXPECT_NE(idx2, idx1);

    idx1.add_record(min(kh1.first, kh1.second), 1, p, 0, 0);
    EXPECT_NE(idx1, idx2);
    EXPECT_NE(idx2, idx1);

    idx1.add_record(min(kh1.first, kh1.second), 3, p, 0, 0);
    EXPECT_NE(idx1, idx2);
    EXPECT_NE(idx2, idx1);
}

TEST(IndexTest, merging_indexes)
{
    uint32_t w = 2, k = 3;
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    auto index = std::make_shared<Index>();
    auto outdir = TEST_CASE_DIR + "kgs/";

    read_prg_file(prgs, TEST_CASE_DIR + "prg1.fa", 1);
    index->index_prgs(prgs, w, k, outdir);
    index->save(TEST_CASE_DIR + "prg1.fa.idx");

    prgs.clear();
    index->clear();
    read_prg_file(prgs, TEST_CASE_DIR + "prg2.fa", 2);
    index->index_prgs(prgs, w, k, outdir);
    index->save(TEST_CASE_DIR + "prg2.fa.idx");

    prgs.clear();
    index->clear();
    read_prg_file(prgs, TEST_CASE_DIR + "prg3.fa", 3);
    index->index_prgs(prgs, w, k, outdir);
    index->save(TEST_CASE_DIR + "prg3.fa.idx");

    // merge
    auto index_merged = std::make_shared<Index>();
    index_merged->load(TEST_CASE_DIR + "prg1.fa.idx");
    index_merged->load(TEST_CASE_DIR + "prg2.fa.idx");
    index_merged->load(TEST_CASE_DIR + "prg3.fa.idx");

    // now an index from all 4 in
    prgs.clear();
    auto index_all = std::make_shared<Index>();
    read_prg_file(prgs, TEST_CASE_DIR + "prg0123.fa");
    index_all->index_prgs(prgs, w, k, outdir);
}
