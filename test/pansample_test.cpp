#include "gtest/gtest.h"
#include "pannode.h"
#include "panedge.h"
#include "pansample.h"
#include <stdint.h>
#include <iostream>

using namespace std;

class PanSampleTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};

TEST_F(PanSampleTest,create){

    PanSample ps("sample");
    EXPECT_EQ("sample", ps.name);
    EXPECT_EQ((uint)0, ps.edges.size());
    EXPECT_EQ((uint)0, ps.paths.size());

    // do the same creating a pointer
    PanSample *ps1;
    ps1 = new PanSample("sample");
    EXPECT_EQ("sample", ps1->name);
    EXPECT_EQ((uint)0, ps1->edges.size());
    EXPECT_EQ((uint)0, ps1->paths.size());
    delete ps1;
}

TEST_F(PanSampleTest, add_path)
{
    PanSample ps("sample");
    std::vector<KmerNode*> kmp;
    ps.add_path(2,kmp);
    EXPECT_EQ((uint)1, ps.paths.size());
    EXPECT_EQ((uint)1, ps.paths[2].size());

    ps.add_path(2,kmp);
    EXPECT_EQ((uint)1, ps.paths.size());
    EXPECT_EQ((uint)2, ps.paths[2].size());

    ps.add_path(3,kmp);
    EXPECT_EQ((uint)2, ps.paths.size());
    EXPECT_EQ((uint)2, ps.paths[2].size());
    EXPECT_EQ((uint)1, ps.paths[3].size());
}

TEST_F(PanSampleTest,equals){
    PanSample ps1("1");
    PanSample ps2("2");
    EXPECT_EQ(ps1, ps1);
    EXPECT_EQ(ps2, ps2);
    EXPECT_EQ((ps1==ps2), false);
    EXPECT_EQ((ps2==ps1), false);   
}

TEST_F(PanSampleTest,nequals){
    PanSample ps1("1");
    PanSample ps2("2");
    EXPECT_EQ((ps1!=ps1), false);
    EXPECT_EQ((ps2!=ps2), false);
    EXPECT_EQ((ps1!=ps2), true);
    EXPECT_EQ((ps2!=ps1), true);
}

TEST_F(PanSampleTest,less){
    PanSample ps1("1");
    PanSample ps2("2");
    EXPECT_EQ((ps1<ps1), false);
    EXPECT_EQ((ps2<ps2), false);
    EXPECT_EQ((ps1<ps2), true);
    EXPECT_EQ((ps2<ps1), false);
}
