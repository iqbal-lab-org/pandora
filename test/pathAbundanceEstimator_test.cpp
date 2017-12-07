#include "gtest/gtest.h"
#include "test_macro.cpp"
#include <cmath>
#include "pathAbundanceEstimator.h"

class PathAbundanceEstimatorTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }
  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(PathAbundanceEstimatorTest, constructor) {
  /*
    List of path IDs and their length:
    0 --> 5
    1 --> 3
    2 --> 7
  */
  std::vector<std::deque<KmerNodePtr>> paths;
  KmerNodePtr dummyNode;
  std::deque<KmerNodePtr> p1;
  p1.push_front(dummyNode);
  p1.push_front(dummyNode);
  p1.push_front(dummyNode);
  p1.push_front(dummyNode);
  p1.push_front(dummyNode);
  std::deque<KmerNodePtr> p2;
  p2.push_front(dummyNode);
  p2.push_front(dummyNode);
  p2.push_front(dummyNode);
  std::deque<KmerNodePtr> p3;
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);

  paths.push_back(p1);
  paths.push_back(p2);
  paths.push_back(p3);


  /*
    Here we have two reads.
    Each read has list of <path IDs and # of hits>
         | 0, 4 |       |2, 6|
    r1 = |      |  r2 = |0, 5|
         | 1, 2 |       |1, 1|

    Expected output:
         | 0, 4 |       |2, 6|
    r1 = |      |  r2 = |0, 5|
         | 1, 2 |       |1, 1|
   */

  double eps = 0.0001;
  std::vector<std::vector<std::pair<uint16_t, uint16_t>>> hitCntPerRead4Paths;
  std::vector<std::pair<uint16_t, uint16_t>> read1Hits;
  read1Hits.emplace_back(0,4);
  read1Hits.emplace_back(1,2);
  hitCntPerRead4Paths.push_back(read1Hits);

  std::vector<std::pair<uint16_t, uint16_t>> read2Hits;
  read2Hits.emplace_back(2,6);
  read2Hits.emplace_back(0,5);
  read2Hits.emplace_back(1,1);
  hitCntPerRead4Paths.push_back(read2Hits);


  PathAbundanceEstimator pae(hitCntPerRead4Paths, paths);

  EXPECT_EQ(pae.readProbs.size(), static_cast<size_t>(2));

  // validate probabilites for read 1
  EXPECT_EQ(pae.readProbs[0].size(), static_cast<size_t>(2));
  EXPECT_EQ(pae.readProbs[0][0].first, 0); // first path ID for the first read
  EXPECT_LE(std::abs(pae.readProbs[0][0].second - 4.0/5.0), eps); // first path prob for the first read
  EXPECT_EQ(pae.readProbs[0][1].first, 1); // second path ID for the first read
  EXPECT_LE(std::abs(pae.readProbs[0][1].second - 2.0/3.0), eps); // second path prob for the first read

  // validate probabilities for read 2
  EXPECT_EQ(pae.readProbs[1].size(), static_cast<size_t>(3));
  EXPECT_EQ(pae.readProbs[1][0].first, 2); // first path ID for the second read
  EXPECT_LE(std::abs(pae.readProbs[1][0].second - 6.0/7.0), eps); // first path prob for the second read
  EXPECT_EQ(pae.readProbs[1][1].first, 0); // second path ID for the second read
  EXPECT_LE(std::abs(pae.readProbs[1][1].second - 5.0/5.0), eps); // second path prob for the second read
  EXPECT_EQ(pae.readProbs[1][2].first, 1); // third path ID for the second read
  EXPECT_LE(std::abs(pae.readProbs[1][2].second - 1.0/3.0), eps); // third path prob for the second read

  // validate total number of pathCnts and that all are initialized to 1
  EXPECT_EQ(pae.pathCnts.size(), static_cast<size_t>(3));
  for (auto const & p :  pae.pathCnts)
    EXPECT_EQ(p, 1.0);

}

TEST_F(PathAbundanceEstimatorTest, runEM) {
    /*
    List of path IDs and their length:
    0 --> 5
    1 --> 3
    2 --> 7
  */
  std::vector<std::deque<KmerNodePtr>> paths;
  KmerNodePtr dummyNode;
  std::deque<KmerNodePtr> p1;
  p1.push_front(dummyNode);
  p1.push_front(dummyNode);
  p1.push_front(dummyNode);
  p1.push_front(dummyNode);
  p1.push_front(dummyNode);
  std::deque<KmerNodePtr> p2;
  p2.push_front(dummyNode);
  p2.push_front(dummyNode);
  p2.push_front(dummyNode);
  std::deque<KmerNodePtr> p3;
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);
  p3.push_front(dummyNode);

  paths.push_back(p1);
  paths.push_back(p2);
  paths.push_back(p3);


  /*
    Here we have two reads.
    Each read has list of <path IDs and # of hits>
         | 0, 4 |       |2, 6|
    r1 = |      |  r2 = |0, 5|
         | 1, 2 |       |1, 1|

    Expected output:
         | 0, 4 |       |2, 6|
    r1 = |      |  r2 = |0, 5|
         | 1, 2 |       |1, 1|
   */

  double eps = 0.0001;
  std::vector<std::vector<std::pair<uint16_t, uint16_t>>> hitCntPerRead4Paths;
  std::vector<std::pair<uint16_t, uint16_t>> read1Hits;
  read1Hits.emplace_back(0,4);
  read1Hits.emplace_back(1,2);
  hitCntPerRead4Paths.push_back(read1Hits);

  std::vector<std::pair<uint16_t, uint16_t>> read2Hits;
  read2Hits.emplace_back(2,6);
  read2Hits.emplace_back(0,5);
  read2Hits.emplace_back(1,1);
  hitCntPerRead4Paths.push_back(read2Hits);


  PathAbundanceEstimator pae(hitCntPerRead4Paths, paths, 1e-8, 1);
  std::vector<double> pathCnts = pae.runEM();

  EXPECT_EQ(pathCnts.size(), static_cast<size_t>(3));
  EXPECT_LE(std::abs(pathCnts[0]-((12.0/22.0)+(105.0/230.0))), eps);
  EXPECT_LE(std::abs(pathCnts[1]-((10.0/22.0)+(35.0/230.0))), eps);
  EXPECT_LE(std::abs(pathCnts[2]-(90.0/230.0)), eps);

}
