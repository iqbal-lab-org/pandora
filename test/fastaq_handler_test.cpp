#include <iostream>
#include <fstream>
#include "gtest/gtest.h"
#include "fastaq_handler.h"

using namespace std;

TEST(FastaqHandlerTest, create_fa) {
    FastaqHandler fh("../../test/test_cases/reads.fa");
    EXPECT_EQ((uint)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file.is_open());
}

TEST(FastaqHandlerTest, create_fq) {
    FastaqHandler fh("../../test/test_cases/reads.fq");
    EXPECT_EQ((uint)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file.is_open());
}

TEST(FastaqHandlerTest, create_fagz) {
    FastaqHandler fh("../../test/test_cases/reads.fa.gz");
    EXPECT_EQ((uint)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file.is_open());
}

TEST(FastaqHandlerTest, create_fqgz) {
    FastaqHandler fh("../../test/test_cases/reads.fq.gz");
    EXPECT_EQ((uint)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file.is_open());
}

TEST(FastaqHandlerTest, getline_fa) {
    FastaqHandler fh("../../test/test_cases/reads.fa");
    bool foundline = false;
    while(std::getline(fh.instream, fh.line)) {
        std::cout << fh.line << std::endl;
        foundline = true;
    }
    EXPECT_TRUE(foundline);
}

TEST(FastaqHandlerTest, getline_fagz) {
    FastaqHandler fh("../../test/test_cases/reads.fa.gz");
    bool foundline = false;
    while(std::getline(fh.instream, fh.line)) {
        std::cout << fh.line << std::endl;
        foundline = true;
    }
    EXPECT_TRUE(foundline);
}

TEST(FastaqHandlerTest, get_next) {
    FastaqHandler fh("../../test/test_cases/reads.fa");
    fh.get_next();
    EXPECT_EQ((uint)1, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read0");
    EXPECT_EQ(fh.read, "to be ignored");

    fh.get_next();
    EXPECT_EQ((uint)2, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read1");
    EXPECT_EQ(fh.read, "should copy the phrase *should*");

    fh.get_next();
    EXPECT_EQ((uint)3, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read2");
    EXPECT_EQ(fh.read, "this time we should get *is time *");

    fh.get_next();
    EXPECT_EQ((uint)3, fh.num_reads_parsed);
    EXPECT_EQ(fh.name, "read2");
    EXPECT_EQ(fh.read, "this time we should get *is time *");
}

TEST(FastaqHandlerTest, get_id_fa) {
    FastaqHandler fh("../../test/test_cases/reads.fa");

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

TEST(FastaqHandlerTest, get_id_fq) {
    FastaqHandler fh("../../test/test_cases/reads.fq");

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

TEST(FastaqHandlerTest, get_id_fagz) {
    FastaqHandler fh("../../test/test_cases/reads.fa.gz");

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

TEST(FastaqHandlerTest, get_id_fqgz) {
    FastaqHandler fh("../../test/test_cases/reads.fq.gz");

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

TEST(FastaqHandlerTest, close) {
    FastaqHandler fh("../../test/test_cases/reads.fa");
    EXPECT_EQ((uint)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file.is_open());
    fh.close();
    EXPECT_FALSE(fh.fastaq_file.is_open());
}

TEST(FastaqHandlerTest, close_fqgz) {
    FastaqHandler fh("../../test/test_cases/reads.fq.gz");
    EXPECT_EQ((uint)0, fh.num_reads_parsed);
    EXPECT_TRUE(fh.fastaq_file.is_open());
    fh.close();
    EXPECT_FALSE(fh.fastaq_file.is_open());
}
