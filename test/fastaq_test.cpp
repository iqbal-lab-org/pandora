#include "gtest/gtest.h"
#include "fastaq.h"
#include <iostream>
#include "test_helpers_containers.h"
#include "test_helpers.h"

using namespace std;

TEST(FastaqTest, create_null)
{
    Fastaq f1;
    EXPECT_FALSE(f1.gzipped);
    EXPECT_FALSE(f1.fastq);
    EXPECT_EQ(f1.names.size(), (uint)0);
    EXPECT_EQ(f1.sequences.size(), (uint)0);
    EXPECT_EQ(f1.scores.size(), (uint)0);
}

TEST(FastaqTest, create_with_args)
{
    Fastaq f1(true, false);
    EXPECT_TRUE(f1.gzipped);
    EXPECT_FALSE(f1.fastq);
    EXPECT_EQ(f1.names.size(), (uint)0);
    EXPECT_EQ(f1.sequences.size(), (uint)0);
    EXPECT_EQ(f1.scores.size(), (uint)0);

    Fastaq f2(false, true);
    EXPECT_FALSE(f2.gzipped);
    EXPECT_TRUE(f2.fastq);
    EXPECT_EQ(f2.names.size(), (uint)0);
    EXPECT_EQ(f2.sequences.size(), (uint)0);
    EXPECT_EQ(f2.scores.size(), (uint)0);

    Fastaq f3(true, true);
    EXPECT_TRUE(f3.gzipped);
    EXPECT_TRUE(f3.fastq);
    EXPECT_EQ(f3.names.size(), (uint)0);
    EXPECT_EQ(f3.sequences.size(), (uint)0);
    EXPECT_EQ(f3.scores.size(), (uint)0);
}

TEST(FastaqTest, covg_to_score)
{
    Fastaq f;
    char ascii_range[] = "!\"#$%&'()*+,-./"
                         "0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`"
                         "abcdefghijklmnopqrstuvwxyz{|}~";
    for (uint i = 0; i < 40; ++i) {
        EXPECT_EQ(f.covg_to_score(i, 40), ascii_range[i]);
    }
}

TEST(FastaqTest, covg_to_score_with_rounding)
{
    Fastaq f;
    char ascii_range[] = "!\"#$%&'()*+,-./"
                         "0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`"
                         "abcdefghijklmnopqrstuvwxyz{|}~";
    for (uint i = 0; i < 40; ++i) {
        EXPECT_EQ(f.covg_to_score(3 * i, 119), ascii_range[i]);
    }
}

TEST(AltCovgToScore, CovgToScoreWithAltTrue_ReturnAltCovgToScoreResult)
{
    Fastaq f;
    const uint_least16_t covg { 0 };

    const auto result { f.covg_to_score(covg, 0, true) };
    const char expected { '!' };

    EXPECT_EQ(result, expected);
}

TEST(AltCovgToScore, CovgZero_ReturnFirstPrintableAscii)
{
    Fastaq f;
    const uint_least16_t covg { 0 };

    const auto result { f.alt_covg_to_score(covg) };
    const char expected { '!' };

    EXPECT_EQ(result, expected);
}

TEST(AltCovgToScore, CovgFive_ReturnSixthPrintableAscii)
{
    Fastaq f;
    const uint_least16_t covg { 5 };

    const auto result { f.alt_covg_to_score(covg) };
    const char expected { '&' };

    EXPECT_EQ(result, expected);
}

TEST(AltCovgToScore, CovgNinetyThree_ReturnLastPrintableAscii)
{
    Fastaq f;
    const uint_least16_t covg { 93 };

    const auto result { f.alt_covg_to_score(covg) };
    const char expected { '~' };

    EXPECT_EQ(result, expected);
}

TEST(AltCovgToScore, CovgNinetyFour_ReturnLastPrintableAscii)
{
    Fastaq f;
    const uint_least16_t covg { 94 };

    const auto result { f.alt_covg_to_score(covg) };
    const char expected { '~' };

    EXPECT_EQ(result, expected);
}

TEST(AltCovgToScore, CovgNinetyTwo_ReturnSecondLastPrintableAscii)
{
    Fastaq f;
    const uint_least16_t covg { 92 };

    const auto result { f.alt_covg_to_score(covg) };
    const char expected { '}' };

    EXPECT_EQ(result, expected);
}

TEST(AltCovgToScore, CrazyHighCovg_ReturnLastPrintableAscii)
{
    Fastaq f;
    const uint_least16_t covg { 920 };

    const auto result { f.alt_covg_to_score(covg) };
    const char expected { '~' };

    EXPECT_EQ(result, expected);
}

TEST(FastaqTest, add_entry_FatalRuntimeError)
{
    Fastaq f;
    ASSERT_EXCEPTION(f.add_entry("", "ACGT", { 0, 1, 2, 3 }, 40), FatalRuntimeError,
        "Error adding entry to Fasta/q file");
    ASSERT_EXCEPTION(f.add_entry("dummy", "ACGT", { 0, 1, 2 }, 40), FatalRuntimeError,
        "Error adding entry to Fasta/q file");
    ASSERT_EXCEPTION(f.add_entry("dummy", "ACG", { 0, 1, 2, 3 }, 40), FatalRuntimeError,
        "Error adding entry to Fasta/q file");
}

TEST(FastaqTest, add_entry_works)
{
    Fastaq f;
    f.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    bool found_name = find(f.names.begin(), f.names.end(), "dummy") != f.names.end();
    EXPECT_TRUE(found_name);
    bool added_seq = f.sequences.find("dummy") != f.sequences.end();
    EXPECT_TRUE(added_seq);
    EXPECT_EQ(f.sequences["dummy"], "ACGTA");
    bool added_score = f.scores.find("dummy") != f.scores.end();
    EXPECT_TRUE(added_score);
    EXPECT_EQ(f.scores["dummy"], "#$%&'");
}

TEST(FastaqTest, different_fastaq_equals_false)
{
    Fastaq f1(false, true);
    f1.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    Fastaq f2(false, false);
    f2.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    EXPECT_FALSE(f1 == f2);
    EXPECT_FALSE(f2 == f1);
}

TEST(FastaqTest, different_gzipped_equals_true)
{
    Fastaq f1(true, true);
    f1.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    Fastaq f2(false, true);
    f2.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    EXPECT_TRUE(f1 == f2);
    EXPECT_TRUE(f2 == f1);
}

TEST(FastaqTest, different_names_equals_false)
{
    Fastaq f1(false, true);
    f1.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    Fastaq f2(false, true);
    f2.add_entry("dummer", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    EXPECT_FALSE(f1 == f2);
    EXPECT_FALSE(f2 == f1);
}

TEST(FastaqTest, different_num_seqs_equals_false)
{
    Fastaq f1(false, true);
    f1.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    Fastaq f2(false, true);
    f2.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    f2.add_entry("dummer", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    EXPECT_FALSE(f1 == f2);
    EXPECT_FALSE(f2 == f1);
}

TEST(FastaqTest, different_seqs_equals_false)
{
    Fastaq f1(false, true);
    f1.add_entry("dummy", "ACGTT", { 2, 3, 4, 5, 6 }, 40);
    Fastaq f2(false, true);
    f2.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    EXPECT_FALSE(f1 == f2);
    EXPECT_FALSE(f2 == f1);
}

TEST(FastaqTest, different_scores_equals_false)
{
    Fastaq f1(false, true);
    f1.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 7 }, 40);
    Fastaq f2(false, true);
    f2.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    EXPECT_FALSE(f1 == f2);
    EXPECT_FALSE(f2 == f1);
}

TEST(FastaqTest, same_equals_false)
{
    Fastaq f1(false, true);
    f1.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    Fastaq f2(false, true);
    f2.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    EXPECT_TRUE(f1 == f2);
    EXPECT_TRUE(f2 == f1);
}

TEST(FastaqTest, ostream)
{
    Fastaq f_out(false, true);
    f_out.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
}

TEST(FastaqTest, istream_fq)
{
    Fastaq f_in;
    istringstream is("@dummy\nACGTA\n+\n#$%&'");
    is >> f_in;

    EXPECT_TRUE(f_in.fastq);
    EXPECT_FALSE(f_in.gzipped);
    bool found_name
        = find(f_in.names.begin(), f_in.names.end(), "dummy") != f_in.names.end();
    EXPECT_TRUE(found_name);
    bool added_seq = f_in.sequences.find("dummy") != f_in.sequences.end();
    EXPECT_TRUE(added_seq);
    EXPECT_EQ(f_in.sequences["dummy"], "ACGTA");
    bool added_score = f_in.scores.find("dummy") != f_in.scores.end();
    EXPECT_TRUE(added_score);
    EXPECT_EQ(f_in.scores["dummy"], "#$%&'");
}

TEST(FastaqTest, istream_fa)
{
    Fastaq f_in;
    istringstream is(">dummy\nACGTA\n>dummer\nGTGGC");
    is >> f_in;

    EXPECT_FALSE(f_in.fastq);
    EXPECT_FALSE(f_in.gzipped);
    bool found_name
        = find(f_in.names.begin(), f_in.names.end(), "dummy") != f_in.names.end();
    EXPECT_TRUE(found_name);
    bool added_seq = f_in.sequences.find("dummy") != f_in.sequences.end();
    EXPECT_TRUE(added_seq);
    EXPECT_EQ(f_in.sequences["dummy"], "ACGTA");
    bool added_score = f_in.scores.find("dummy") != f_in.scores.end();
    EXPECT_FALSE(added_score);
    found_name
        = find(f_in.names.begin(), f_in.names.end(), "dummer") != f_in.names.end();
    EXPECT_TRUE(found_name);
    added_seq = f_in.sequences.find("dummer") != f_in.sequences.end();
    EXPECT_TRUE(added_seq);
    EXPECT_EQ(f_in.sequences["dummer"], "GTGGC");
}

TEST(FastaqTest, istream_fa_with_extra_header)
{
    Fastaq f_in;
    istringstream is(">dummy with header\nACGTA\n>dummer also with header\nGTGGC");
    is >> f_in;

    EXPECT_FALSE(f_in.fastq);
    EXPECT_FALSE(f_in.gzipped);
    bool found_name
        = find(f_in.names.begin(), f_in.names.end(), "dummy") != f_in.names.end();
    EXPECT_TRUE(found_name);
    bool added_seq = f_in.sequences.find("dummy") != f_in.sequences.end();
    EXPECT_TRUE(added_seq);
    EXPECT_EQ(f_in.sequences["dummy"], "ACGTA");
    bool added_score = f_in.scores.find("dummy") != f_in.scores.end();
    EXPECT_FALSE(added_score);
    bool added_header = f_in.headers.find("dummy") != f_in.headers.end();
    EXPECT_TRUE(added_header);
    EXPECT_EQ(f_in.headers["dummy"], " with header");
    found_name
        = find(f_in.names.begin(), f_in.names.end(), "dummer") != f_in.names.end();
    EXPECT_TRUE(found_name);
    added_seq = f_in.sequences.find("dummer") != f_in.sequences.end();
    EXPECT_TRUE(added_seq);
    EXPECT_EQ(f_in.sequences["dummer"], "GTGGC");
    added_header = f_in.headers.find("dummer") != f_in.headers.end();
    EXPECT_TRUE(added_header);
    EXPECT_EQ(f_in.headers["dummer"], " also with header");
}

TEST(FastaqTest, iostream)
{
    Fastaq f_out(false, true);
    f_out.add_entry("dummy", "ACGTA", { 2, 3, 4, 5, 6 }, 40);
    Fastaq f_in;
    stringstream out;
    out << f_out;
    out >> f_in;

    EXPECT_TRUE(f_in.fastq);
    EXPECT_FALSE(f_in.gzipped);
    bool found_name
        = find(f_in.names.begin(), f_in.names.end(), "dummy") != f_in.names.end();
    EXPECT_TRUE(found_name);
    bool added_seq = f_in.sequences.find("dummy") != f_in.sequences.end();
    EXPECT_TRUE(added_seq);
    EXPECT_EQ(f_in.sequences["dummy"], "ACGTA");
    bool added_score = f_in.scores.find("dummy") != f_in.scores.end();
    EXPECT_TRUE(added_score);
    EXPECT_EQ(f_in.scores["dummy"], "#$%&'");
}
