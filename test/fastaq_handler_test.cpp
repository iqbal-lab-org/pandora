#include <iostream>
#include <fstream>
#include <cstdint>
#include "gtest/gtest.h"
#include "fastaq_handler.h"

using namespace std;

const std::string TEST_CASE_DIR = "../../test/test_cases/";

TEST(FastaqHandlerTest, create_fa)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fa");
    EXPECT_EQ((uint32_t)0, fh.num_reads_parsed);

    EXPECT_TRUE(fh.fastaq_file);
}

TEST(FastaqHandlerTest, create_fq)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fq");
    EXPECT_EQ((uint)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file);
}

TEST(FastaqHandlerTest, create_fagz)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fa.gz");
    EXPECT_EQ((uint)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file);
}

TEST(FastaqHandlerTest, create_fqgz)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fq.gz");
    EXPECT_EQ((uint)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file);
}

TEST(FastaqHandlerTest, get_next)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fa");
    fh.get_next();
    EXPECT_EQ((uint32_t)1, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read0");
    EXPECT_EQ(fh.read, "to be ignored");

    fh.get_next();
    EXPECT_EQ((uint32_t)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_next();
    EXPECT_EQ((uint32_t)3, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read2");
    EXPECT_EQ(fh.read, "this time we should get *is time *");

    fh.get_next();
    EXPECT_EQ((uint)4, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read3");
    EXPECT_EQ(fh.read, "nonsense");

    fh.get_next();
    EXPECT_EQ((uint)5, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read4");
    EXPECT_EQ(fh.read, "another junk line");

    fh.get_next();
    EXPECT_EQ((uint)5, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read4");
    EXPECT_EQ(fh.read, "another junk line");
}

TEST(FastaqHandlerTest, eof)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fa");
    EXPECT_FALSE(fh.eof());
    fh.get_next();
    EXPECT_FALSE(fh.eof());
    fh.get_next();
    EXPECT_FALSE(fh.eof());
    fh.get_next();
    EXPECT_FALSE(fh.eof());
    fh.get_next();
    EXPECT_FALSE(fh.eof());
    fh.get_next();
    EXPECT_TRUE(fh.eof());
    fh.get_next();
    EXPECT_TRUE(fh.eof());
}

TEST(FastaqHandlerTest, get_id_fa)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fa");

    fh.get_id(1);
    EXPECT_EQ((uint32_t)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_id(0);
    EXPECT_EQ((uint32_t)1, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read0");
    EXPECT_EQ(fh.read, "to be ignored");

    fh.get_id(2);
    EXPECT_EQ((uint32_t)3, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read2");
    EXPECT_EQ(fh.read, "this time we should get *is time *");

    fh.get_id(1);
    EXPECT_EQ((uint32_t)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_id(0);
    EXPECT_EQ((uint32_t)1, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read0");
    EXPECT_EQ(fh.read, "to be ignored");

    fh.get_id(1);
    EXPECT_EQ((uint32_t)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_id(2);
    EXPECT_EQ((uint32_t)3, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read2");
    EXPECT_EQ(fh.read, "this time we should get *is time *");
}

TEST(FastaqHandlerTest, get_id_fq)
{
    const std::string filepath = std::tmpnam(nullptr);
    {
        std::ofstream outstream(filepath);
        outstream << "@read1 comment\nACGT\n+\n^^^^\n";
        outstream << "@read2 comment\nGCGT\n+\n^^^^\n";
        outstream << "@read3 comment\nTCGT\n+\n^^^^\n";
    }

    FastaqHandler fh(filepath);
    fh.get_id(1);
    uint32_t expected_num_parsed = 2;
    std::string expected_name = "read2";
    std::string expected_read = "GCGT";
    EXPECT_EQ(expected_num_parsed, fh.num_reads_parsed);
    EXPECT_EQ(expected_name, fh.name);
    EXPECT_EQ(expected_read, fh.read);

    fh.get_id(0);
    expected_num_parsed = 1;
    expected_name = "read1";
    expected_read = "ACGT";
    EXPECT_EQ(expected_num_parsed, fh.num_reads_parsed);
    EXPECT_EQ(expected_name, fh.name);
    EXPECT_EQ(expected_read, fh.read);

    fh.get_id(2);
    expected_num_parsed = 3;
    expected_name = "read3";
    expected_read = "TCGT";
    EXPECT_EQ(expected_num_parsed, fh.num_reads_parsed);
    EXPECT_EQ(expected_name, fh.name);
    EXPECT_EQ(expected_read, fh.read);
}

TEST(FastaqHandlerTest, get_id_past_end_throws_error)
{
    const std::string filepath = std::tmpnam(nullptr);
    {
        std::ofstream outstream(filepath);
        outstream << "@read1 comment\nACGT\n+\n^^^^\n";
        outstream << "@read2 comment\nGCGT\n+\n^^^^\n";
    }

    FastaqHandler fh(filepath);
    EXPECT_THROW(fh.get_id(10), std::out_of_range);
}

TEST(FastaqHandlerTest, get_id_fagz)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fa.gz");

    fh.get_id(1);
    EXPECT_EQ((uint)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_id(0);
    EXPECT_EQ((uint)1, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read0");
    EXPECT_EQ(fh.read, "to be ignored");

    fh.get_id(2);
    EXPECT_EQ((uint)3, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read2");
    EXPECT_EQ(fh.read, "this time we should get *is time *");

    fh.get_id(1);
    EXPECT_EQ((uint)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_id(0);
    EXPECT_EQ((uint)1, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read0");
    EXPECT_EQ(fh.read, "to be ignored");

    fh.get_id(1);
    EXPECT_EQ((uint)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_id(2);
    EXPECT_EQ((uint)3, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read2");
    EXPECT_EQ(fh.read, "this time we should get *is time *");
}

TEST(FastaqHandlerTest, get_id_fqgz)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fq.gz");

    fh.get_id(1);
    EXPECT_EQ((uint)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_id(0);
    EXPECT_EQ((uint)1, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read0");
    EXPECT_EQ(fh.read, "to be ignored");

    fh.get_id(2);
    EXPECT_EQ((uint)3, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read2");
    EXPECT_EQ(fh.read, "this time we should get *is time *");

    fh.get_id(1);
    EXPECT_EQ((uint)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_id(0);
    EXPECT_EQ((uint)1, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read0");
    EXPECT_EQ(fh.read, "to be ignored");

    fh.get_id(1);
    EXPECT_EQ((uint)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_id(2);
    EXPECT_EQ((uint)3, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read2");
    EXPECT_EQ(fh.read, "this time we should get *is time *");
}

TEST(FastaqHandlerTest, close)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fa");
    EXPECT_EQ((uint32_t)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file);
    fh.close();
    EXPECT_TRUE(fh.is_closed());
}

TEST(FastaqHandlerTest, close_fqgz)
{
    FastaqHandler fh(TEST_CASE_DIR + "reads.fq.gz");
    EXPECT_EQ((uint)0, fh.num_reads_parsed);
    EXPECT_FALSE(fh.is_closed());
    fh.close();
    EXPECT_TRUE(fh.is_closed());
}
