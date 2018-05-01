#include <iostream>
#include <fstream>
#include <cstdint>
#include "gtest/gtest.h"
#include "fastaq_handler.h"

using namespace std;

TEST(FastaqHandlerTest, create) {
    FastaqHandler fh("../../test/test_cases/reads.fa");
    EXPECT_EQ((uint32_t)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file.is_open());
}

TEST(FastaqHandlerTest, get_next) {
    FastaqHandler fh("../../test/test_cases/reads.fa");
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
    EXPECT_EQ((uint32_t)3, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read2");
    EXPECT_EQ(fh.read, "this time we should get *is time *");
}

TEST(FastaqHandlerTest, get_id) {
    FastaqHandler fh("../../test/test_cases/reads.fa");

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

TEST(FastaqHandlerTest, get_id_fq) {
    FastaqHandler fh("../../test/test_cases/reads.fq");

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

TEST(FastaqHandlerTest, close) {
    FastaqHandler fh("../../test/test_cases/reads.fa");
    EXPECT_EQ((uint32_t)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file.is_open());
    fh.close();
    EXPECT_FALSE(fh.fastaq_file.is_open());
}
