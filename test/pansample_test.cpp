#include "gtest/gtest.h"
#include "pangenome/pannode.h"
#include "pangenome/pansample.h"
#include <stdint.h>
#include <iostream>

using namespace pangenome;

TEST(PangenomeSampleTest, create)
{

    Sample ps("sample", 0);
    EXPECT_EQ("sample", ps.name);
    EXPECT_EQ((uint)0, ps.paths.size());

    // do the same creating a pointer
    SamplePtr ps1(std::make_shared<Sample>("sample", 0));
    EXPECT_EQ("sample", ps1->name);
    EXPECT_EQ((uint)0, ps1->paths.size());
}

TEST(PangenomeSampleTest, add_path)
{
    Sample ps("sample", 0);
    std::vector<KmerNodePtr> kmp;
    ps.add_path(2, kmp);
    EXPECT_EQ((uint)1, ps.paths.size());
    EXPECT_EQ((uint)1, ps.paths[2].size());

    ps.add_path(2, kmp);
    EXPECT_EQ((uint)1, ps.paths.size());
    EXPECT_EQ((uint)2, ps.paths[2].size());

    ps.add_path(3, kmp);
    EXPECT_EQ((uint)2, ps.paths.size());
    EXPECT_EQ((uint)2, ps.paths[2].size());
    EXPECT_EQ((uint)1, ps.paths[3].size());
}

TEST(PangenomeSampleTest, equals)
{
    Sample ps1("1", 0);
    Sample ps2("2", 0);
    EXPECT_EQ(ps1, ps1);
    EXPECT_EQ(ps2, ps2);
    EXPECT_EQ((ps1 == ps2), false);
    EXPECT_EQ((ps2 == ps1), false);
}

TEST(PangenomeSampleTest, nequals)
{
    Sample ps1("1", 0);
    Sample ps2("2", 0);
    EXPECT_EQ((ps1 != ps1), false);
    EXPECT_EQ((ps2 != ps2), false);
    EXPECT_EQ((ps1 != ps2), true);
    EXPECT_EQ((ps2 != ps1), true);
}

TEST(PangenomeSampleTest, less)
{
    Sample ps1("1", 0);
    Sample ps2("2", 0);
    EXPECT_EQ((ps1 < ps1), false);
    EXPECT_EQ((ps2 < ps2), false);
    EXPECT_EQ((ps1 < ps2), true);
    EXPECT_EQ((ps2 < ps1), false);
}
