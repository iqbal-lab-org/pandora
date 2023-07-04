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


// This mocking class just makes some private members in Index public to easen testing
class IndexMock : public Index {
public:
    IndexMock(uint32_t w, uint32_t k, std::vector<std::string> prg_names, std::vector<uint32_t> prg_lengths) :
        Index(w, k, std::move(prg_names), std::move(prg_lengths)) {}
    IndexMock(uint32_t w, uint32_t k, std::vector<std::string> prg_names,
        std::vector<uint32_t> prg_lengths, std::vector<uint32_t> prg_min_path_lengths,
        const std::shared_ptr<ZipFileReader> &zip_file) :
        Index(w, k, std::move(prg_names), std::move(prg_lengths),
            std::move(prg_min_path_lengths), zip_file) {}
    void save_minhash(ZipFileWriter &index_archive) const {
        Index::save_minhash(index_archive);
    }
    void load_minhash() {
        Index::load_minhash();
    }
};


class IndexTest___Fixture : public ::testing::Test {
protected:
    IndexMock* index;
    IndexTest___Fixture() : index(nullptr){}
    void SetUp() override
    {
        std::vector<std::string> prg_names{"prg_1", "prg_2", "prg_3"};
        std::vector<uint32_t> prg_lengths{10, 20, 30};
        index = new IndexMock(1, 5, prg_names, prg_lengths);
    }
    void TearDown() override {
        delete index;
    }
};


TEST_F(IndexTest___Fixture, add_record)
{
    pandora::KmerHash hash;
    deque<Interval> d = { Interval(3, 5), Interval(9, 12) };
    prg::Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    index->add_record(min(kh.first, kh.second), 1, p, 0, 0);
    uint32_t j = 1;
    EXPECT_EQ(j, index->minhash.size());

    // add again - should stay same size
    index->add_record(min(kh.first, kh.second), 1, p, 0, 0);
    EXPECT_EQ(j, index->minhash.size());
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());

    // add a new record with different key
    pair<uint64_t, uint64_t> kh2 = hash.kmerhash("ACTGA", 5);
    index->add_record(min(kh2.first, kh2.second), 2, p, 0, 0);
    j = 2;
    EXPECT_EQ(j, index->minhash.size());

    // and a new record which is different but has same key
    index->add_record(min(kh.first, kh.second), 4, p, 0, 0);
    EXPECT_EQ(j, index->minhash.size());
    EXPECT_EQ(j, index->minhash[min(kh.first, kh.second)]->size());
}

TEST_F(IndexTest___Fixture, clear)
{
    pandora::KmerHash hash;
    deque<Interval> d = { Interval(3, 5), Interval(9, 12) };
    prg::Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    index->add_record(min(kh.first, kh.second), 1, p, 0, 0);
    kh = hash.kmerhash("ACTGA", 5);
    index->add_record(min(kh.first, kh.second), 2, p, 0, 0);
    kh = hash.kmerhash("ACGTA", 5);
    index->add_record(min(kh.first, kh.second), 4, p, 0, 0);
    index->clear();
    uint32_t j = 0;
    EXPECT_EQ(j, index->minhash.size());
}


void build_load_and_getters___check_prg(Index &loaded_index, uint32_t id, const std::string& name, const std::string& seq, const std::string& gfa)
{
    LocalPRG expected_local_prg(id, name, seq);
    expected_local_prg.kmer_prg.load(gfa);
    auto actual_local_prg = loaded_index.get_prg_given_name(name);
    EXPECT_EQ(*actual_local_prg, expected_local_prg);
    actual_local_prg = loaded_index.get_prg_given_id(id);
    EXPECT_EQ(*actual_local_prg, expected_local_prg);
}
TEST_F(IndexTest___Fixture, build_load_and_getters)
{
    const auto index_filepath = TEST_CASE_DIR + "/sample_example/pangenome.prg.fa.panidx.zip";
    const auto prg_filepath = TEST_CASE_DIR + "/sample_example/pangenome.prg.fa";

    Index::build_index_on_disk(14, 15, prg_filepath, index_filepath);
    Index loaded_index = Index::load(index_filepath);
    loaded_index.load_all_prgs();

    EXPECT_EQ(loaded_index.get_kmer_size(), 15);
    EXPECT_EQ(loaded_index.get_window_size(), 14);
    EXPECT_EQ(loaded_index.get_prg_names(), std::vector<std::string>({"GC00006032", "GC00010897", "empty"}));
    EXPECT_EQ(loaded_index.get_prg_name_given_id(0), "GC00006032");
    EXPECT_EQ(loaded_index.get_prg_name_given_id(1), "GC00010897");
    EXPECT_EQ(loaded_index.get_prg_name_given_id(2), "empty");
    EXPECT_EQ(loaded_index.get_prg_min_path_lengths(), std::vector<::uint32_t>({33, 63, 6}));
    EXPECT_EQ(loaded_index.get_prg_min_path_lengths_given_id(0), 33);
    EXPECT_EQ(loaded_index.get_prg_min_path_lengths_given_id(1), 63);
    EXPECT_EQ(loaded_index.get_prg_min_path_lengths_given_id(2), 6);
    EXPECT_EQ(loaded_index.get_prg_lengths(), std::vector<::uint32_t>({312, 491, 40}));
    EXPECT_EQ(loaded_index.get_prg_length_given_id(0), 312);
    EXPECT_EQ(loaded_index.get_prg_length_given_id(1), 491);
    EXPECT_EQ(loaded_index.get_prg_length_given_id(2), 40);
    EXPECT_EQ(loaded_index.minhash.size(), 114);
    EXPECT_EQ(loaded_index.get_number_of_prgs(), 3);
    EXPECT_EQ(loaded_index.get_loaded_prgs().size(), 3);
    build_load_and_getters___check_prg(loaded_index, 0, "GC00006032",
        "TTGAGTAAAACAATCCCCCGCGCTTATATAAGCGCGTTGATATTTTTAATTATTAACAAGCAACATCATGCTAATACAGACATACAAGGAGATCATCTCTCTTTGCCTGTTTTTTATTATTTCAGGAGTGTAAACACATTTTCCG 5 T 6 C 5 CTCCCTGGCTAAT 7 C 8 A 7 ACCACATTGGCATTTATGGAGCACATCACAATATTTCAATACCATTAAAGCACTGCA 9 C 10 T 9 CAAAATGAAACACTGCGA 11 C 12 T 11 ATTAAAATT 13 A 14 C 13 TTTCAATT",
        "H\tVN:Z:1.0\tbn:Z:--linear --singlearr\n"
        "S\t0\t1{[0, 0)}\tFC:i:0\t\tRC:i:0\n"
        "L\t0\t+\t1\t+\t0M\n"
        "S\t1\t1{[9, 24)}\tFC:i:0\t\tRC:i:0\n"
        "L\t1\t+\t2\t+\t0M\n"
        "S\t2\t1{[10, 25)}\tFC:i:0\t\tRC:i:0\n"
        "L\t2\t+\t3\t+\t0M\n"
        "S\t3\t1{[15, 30)}\tFC:i:0\t\tRC:i:0\n"
        "L\t3\t+\t4\t+\t0M\n"
        "S\t4\t1{[25, 40)}\tFC:i:0\t\tRC:i:0\n"
        "L\t4\t+\t5\t+\t0M\n"
        "S\t5\t1{[29, 44)}\tFC:i:0\t\tRC:i:0\n"
        "L\t5\t+\t6\t+\t0M\n"
        "S\t6\t1{[43, 58)}\tFC:i:0\t\tRC:i:0\n"
        "L\t6\t+\t7\t+\t0M\n"
        "S\t7\t1{[52, 67)}\tFC:i:0\t\tRC:i:0\n"
        "L\t7\t+\t8\t+\t0M\n"
        "S\t8\t1{[58, 73)}\tFC:i:0\t\tRC:i:0\n"
        "L\t8\t+\t9\t+\t0M\n"
        "S\t9\t1{[61, 76)}\tFC:i:0\t\tRC:i:0\n"
        "L\t9\t+\t10\t+\t0M\n"
        "S\t10\t1{[66, 81)}\tFC:i:0\t\tRC:i:0\n"
        "L\t10\t+\t11\t+\t0M\n"
        "S\t11\t1{[67, 82)}\tFC:i:0\t\tRC:i:0\n"
        "L\t11\t+\t12\t+\t0M\n"
        "S\t12\t1{[79, 94)}\tFC:i:0\t\tRC:i:0\n"
        "L\t12\t+\t13\t+\t0M\n"
        "S\t13\t1{[82, 97)}\tFC:i:0\t\tRC:i:0\n"
        "L\t13\t+\t14\t+\t0M\n"
        "S\t14\t1{[86, 101)}\tFC:i:0\t\tRC:i:0\n"
        "L\t14\t+\t15\t+\t0M\n"
        "S\t15\t1{[93, 108)}\tFC:i:0\t\tRC:i:0\n"
        "L\t15\t+\t16\t+\t0M\n"
        "S\t16\t1{[105, 120)}\tFC:i:0\t\tRC:i:0\n"
        "L\t16\t+\t17\t+\t0M\n"
        "S\t17\t1{[115, 130)}\tFC:i:0\t\tRC:i:0\n"
        "L\t17\t+\t18\t+\t0M\n"
        "S\t18\t1{[125, 140)}\tFC:i:0\t\tRC:i:0\n"
        "L\t18\t+\t19\t+\t0M\n"
        "S\t19\t1{[129, 144)}\tFC:i:0\t\tRC:i:0\n"
        "L\t19\t+\t20\t+\t0M\n"
        "L\t19\t+\t21\t+\t0M\n"
        "S\t20\t3{[133, 145)[148, 149)[156, 158)}\tFC:i:0\t\tRC:i:0\n"
        "L\t20\t+\t22\t+\t0M\n"
        "L\t20\t+\t23\t+\t0M\n"
        "S\t21\t3{[137, 145)[152, 153)[156, 162)}\tFC:i:0\t\tRC:i:0\n"
        "L\t21\t+\t24\t+\t0M\n"
        "S\t22\t3{[156, 169)[172, 173)[180, 181)}\tFC:i:0\t\tRC:i:0\n"
        "L\t22\t+\t25\t+\t0M\n"
        "S\t23\t3{[138, 145)[148, 149)[156, 163)}\tFC:i:0\t\tRC:i:0\n"
        "L\t23\t+\t26\t+\t0M\n"
        "S\t24\t3{[144, 145)[152, 153)[156, 169)}\tFC:i:0\t\tRC:i:0\n"
        "L\t24\t+\t25\t+\t0M\n"
        "L\t24\t+\t26\t+\t0M\n"
        "S\t25\t3{[158, 169)[172, 173)[180, 183)}\tFC:i:0\t\tRC:i:0\n"
        "L\t25\t+\t27\t+\t0M\n"
        "S\t26\t3{[161, 169)[176, 177)[180, 186)}\tFC:i:0\t\tRC:i:0\n"
        "L\t26\t+\t28\t+\t0M\n"
        "S\t27\t3{[161, 169)[172, 173)[180, 186)}\tFC:i:0\t\tRC:i:0\n"
        "L\t27\t+\t29\t+\t0M\n"
        "S\t28\t3{[166, 169)[176, 177)[180, 191)}\tFC:i:0\t\tRC:i:0\n"
        "L\t28\t+\t30\t+\t0M\n"
        "S\t29\t2{[172, 173)[180, 194)}\tFC:i:0\t\tRC:i:0\n"
        "L\t29\t+\t30\t+\t0M\n"
        "S\t30\t1{[186, 201)}\tFC:i:0\t\tRC:i:0\n"
        "L\t30\t+\t31\t+\t0M\n"
        "S\t31\t1{[200, 215)}\tFC:i:0\t\tRC:i:0\n"
        "L\t31\t+\t32\t+\t0M\n"
        "S\t32\t1{[214, 229)}\tFC:i:0\t\tRC:i:0\n"
        "L\t32\t+\t33\t+\t0M\n"
        "S\t33\t1{[218, 233)}\tFC:i:0\t\tRC:i:0\n"
        "L\t33\t+\t34\t+\t0M\n"
        "L\t33\t+\t35\t+\t0M\n"
        "S\t34\t3{[229, 237)[245, 246)[249, 255)}\tFC:i:0\t\tRC:i:0\n"
        "L\t34\t+\t36\t+\t0M\n"
        "S\t35\t1{[222, 237)}\tFC:i:0\t\tRC:i:0\n"
        "L\t35\t+\t37\t+\t0M\n"
        "S\t36\t3{[230, 237)[245, 246)[249, 256)}\tFC:i:0\t\tRC:i:0\n"
        "L\t36\t+\t38\t+\t0M\n"
        "S\t37\t3{[229, 237)[240, 241)[249, 255)}\tFC:i:0\t\tRC:i:0\n"
        "L\t37\t+\t39\t+\t0M\n"
        "S\t38\t3{[236, 237)[245, 246)[249, 262)}\tFC:i:0\t\tRC:i:0\n"
        "L\t38\t+\t40\t+\t0M\n"
        "L\t38\t+\t41\t+\t0M\n"
        "S\t39\t1{[250, 265)}\tFC:i:0\t\tRC:i:0\n"
        "L\t39\t+\t42\t+\t0M\n"
        "L\t39\t+\t41\t+\t0M\n"
        "S\t40\t3{[259, 267)[271, 272)[281, 287)}\tFC:i:0\t\tRC:i:0\n"
        "L\t40\t+\t44\t+\t0M\n"
        "S\t41\t3{[259, 267)[276, 277)[281, 287)}\tFC:i:0\t\tRC:i:0\n"
        "L\t41\t+\t43\t+\t0M\n"
        "L\t41\t+\t44\t+\t0M\n"
        "S\t42\t2{[253, 267)[271, 272)}\tFC:i:0\t\tRC:i:0\n"
        "L\t42\t+\t40\t+\t0M\n"
        "S\t43\t5{[264, 267)[276, 277)[281, 290)[294, 295)[304, 305)}\tFC:i:0\t\tRC:i:0\n"
        "L\t43\t+\t44\t+\t0M\n"
        "S\t44\t1{[312, 312)}\tFC:i:0\t\tRC:i:0\n");
    build_load_and_getters___check_prg(loaded_index, 1, "GC00010897",
        "ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCACCGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC 5 C 6 G 5 GAAATGCCCACGACCGGGCTGGATGAGCTGACAGAGGCCGAACGCGAGAC 7 A 8 G 7 CTGGCCGAATGGCTGGCCAAGCGCCGGGAAGCCTCGCAGAAGTCGCAGGAGGCCTACACGGCCATGTCTGCGGATCGGTGGCTGGTCACGCTGGCCAAGGCCATCAGGGAAGGGCAGGA 9 GCTA 10 ACTG 9 CGCCCCGAACAGGCGGCCGCGATCTGGCACGGCATGGGGGA 11 A 12 G 11 GTCGGCAAGGCCTTGCGCAAGGCTGGTCACGCGAAGCCCAAGGCGGTCAGAAAGGGCAAGCCGGTCGATCCGGCTGATCCCAAGGATCAAGGGGAGGGGGCACCAAAGGGGAAATGA",
        "H\tVN:Z:1.0\tbn:Z:--linear --singlearr\n"
        "S\t0\t1{[0, 0)}\tFC:i:0\t\tRC:i:0\n"
        "L\t0\t+\t1\t+\t0M\n"
        "S\t1\t1{[6, 21)}\tFC:i:0\t\tRC:i:0\n"
        "L\t1\t+\t2\t+\t0M\n"
        "S\t2\t1{[13, 28)}\tFC:i:0\t\tRC:i:0\n"
        "L\t2\t+\t3\t+\t0M\n"
        "S\t3\t1{[24, 39)}\tFC:i:0\t\tRC:i:0\n"
        "L\t3\t+\t4\t+\t0M\n"
        "S\t4\t1{[25, 40)}\tFC:i:0\t\tRC:i:0\n"
        "L\t4\t+\t5\t+\t0M\n"
        "S\t5\t1{[29, 44)}\tFC:i:0\t\tRC:i:0\n"
        "L\t5\t+\t6\t+\t0M\n"
        "S\t6\t1{[42, 57)}\tFC:i:0\t\tRC:i:0\n"
        "L\t6\t+\t7\t+\t0M\n"
        "S\t7\t1{[49, 64)}\tFC:i:0\t\tRC:i:0\n"
        "L\t7\t+\t8\t+\t0M\n"
        "S\t8\t1{[51, 66)}\tFC:i:0\t\tRC:i:0\n"
        "L\t8\t+\t9\t+\t0M\n"
        "S\t9\t1{[52, 67)}\tFC:i:0\t\tRC:i:0\n"
        "L\t9\t+\t10\t+\t0M\n"
        "S\t10\t1{[57, 72)}\tFC:i:0\t\tRC:i:0\n"
        "L\t10\t+\t11\t+\t0M\n"
        "S\t11\t1{[68, 83)}\tFC:i:0\t\tRC:i:0\n"
        "L\t11\t+\t12\t+\t0M\n"
        "S\t12\t1{[75, 90)}\tFC:i:0\t\tRC:i:0\n"
        "L\t12\t+\t13\t+\t0M\n"
        "S\t13\t1{[76, 91)}\tFC:i:0\t\tRC:i:0\n"
        "L\t13\t+\t14\t+\t0M\n"
        "S\t14\t1{[78, 93)}\tFC:i:0\t\tRC:i:0\n"
        "L\t14\t+\t15\t+\t0M\n"
        "S\t15\t1{[88, 103)}\tFC:i:0\t\tRC:i:0\n"
        "L\t15\t+\t16\t+\t0M\n"
        "S\t16\t1{[89, 104)}\tFC:i:0\t\tRC:i:0\n"
        "L\t16\t+\t17\t+\t0M\n"
        "L\t16\t+\t18\t+\t0M\n"
        "S\t17\t3{[103, 110)[113, 114)[121, 128)}\tFC:i:0\t\tRC:i:0\n"
        "L\t17\t+\t19\t+\t0M\n"
        "S\t18\t3{[99, 110)[117, 118)[121, 124)}\tFC:i:0\t\tRC:i:0\n"
        "L\t18\t+\t20\t+\t0M\n"
        "S\t19\t1{[127, 142)}\tFC:i:0\t\tRC:i:0\n"
        "L\t19\t+\t21\t+\t0M\n"
        "S\t20\t3{[107, 110)[117, 118)[121, 132)}\tFC:i:0\t\tRC:i:0\n"
        "L\t20\t+\t21\t+\t0M\n"
        "S\t21\t1{[130, 145)}\tFC:i:0\t\tRC:i:0\n"
        "L\t21\t+\t22\t+\t0M\n"
        "S\t22\t1{[132, 147)}\tFC:i:0\t\tRC:i:0\n"
        "L\t22\t+\t23\t+\t0M\n"
        "S\t23\t1{[143, 158)}\tFC:i:0\t\tRC:i:0\n"
        "L\t23\t+\t24\t+\t0M\n"
        "S\t24\t1{[147, 162)}\tFC:i:0\t\tRC:i:0\n"
        "L\t24\t+\t25\t+\t0M\n"
        "S\t25\t1{[150, 165)}\tFC:i:0\t\tRC:i:0\n"
        "L\t25\t+\t26\t+\t0M\n"
        "L\t25\t+\t27\t+\t0M\n"
        "S\t26\t3{[159, 171)[174, 175)[182, 184)}\tFC:i:0\t\tRC:i:0\n"
        "L\t26\t+\t28\t+\t0M\n"
        "S\t27\t3{[163, 171)[178, 179)[182, 188)}\tFC:i:0\t\tRC:i:0\n"
        "L\t27\t+\t29\t+\t0M\n"
        "S\t28\t3{[162, 171)[174, 175)[182, 187)}\tFC:i:0\t\tRC:i:0\n"
        "L\t28\t+\t29\t+\t0M\n"
        "S\t29\t1{[184, 199)}\tFC:i:0\t\tRC:i:0\n"
        "L\t29\t+\t30\t+\t0M\n"
        "S\t30\t1{[195, 210)}\tFC:i:0\t\tRC:i:0\n"
        "L\t30\t+\t31\t+\t0M\n"
        "S\t31\t1{[209, 224)}\tFC:i:0\t\tRC:i:0\n"
        "L\t31\t+\t32\t+\t0M\n"
        "S\t32\t1{[220, 235)}\tFC:i:0\t\tRC:i:0\n"
        "L\t32\t+\t33\t+\t0M\n"
        "S\t33\t1{[224, 239)}\tFC:i:0\t\tRC:i:0\n"
        "L\t33\t+\t34\t+\t0M\n"
        "S\t34\t1{[228, 243)}\tFC:i:0\t\tRC:i:0\n"
        "L\t34\t+\t35\t+\t0M\n"
        "S\t35\t1{[239, 254)}\tFC:i:0\t\tRC:i:0\n"
        "L\t35\t+\t36\t+\t0M\n"
        "S\t36\t1{[251, 266)}\tFC:i:0\t\tRC:i:0\n"
        "L\t36\t+\t37\t+\t0M\n"
        "S\t37\t1{[261, 276)}\tFC:i:0\t\tRC:i:0\n"
        "L\t37\t+\t38\t+\t0M\n"
        "S\t38\t1{[262, 277)}\tFC:i:0\t\tRC:i:0\n"
        "L\t38\t+\t39\t+\t0M\n"
        "S\t39\t1{[263, 278)}\tFC:i:0\t\tRC:i:0\n"
        "L\t39\t+\t40\t+\t0M\n"
        "S\t40\t1{[271, 286)}\tFC:i:0\t\tRC:i:0\n"
        "L\t40\t+\t41\t+\t0M\n"
        "S\t41\t1{[273, 288)}\tFC:i:0\t\tRC:i:0\n"
        "L\t41\t+\t42\t+\t0M\n"
        "S\t42\t1{[275, 290)}\tFC:i:0\t\tRC:i:0\n"
        "L\t42\t+\t43\t+\t0M\n"
        "S\t43\t1{[281, 296)}\tFC:i:0\t\tRC:i:0\n"
        "L\t43\t+\t44\t+\t0M\n"
        "L\t43\t+\t45\t+\t0M\n"
        "S\t44\t3{[292, 301)[312, 316)[319, 321)}\tFC:i:0\t\tRC:i:0\n"
        "L\t44\t+\t46\t+\t0M\n"
        "S\t45\t3{[293, 301)[304, 308)[319, 322)}\tFC:i:0\t\tRC:i:0\n"
        "L\t45\t+\t47\t+\t0M\n"
        "S\t46\t3{[293, 301)[312, 316)[319, 322)}\tFC:i:0\t\tRC:i:0\n"
        "L\t46\t+\t48\t+\t0M\n"
        "S\t47\t2{[307, 308)[319, 333)}\tFC:i:0\t\tRC:i:0\n"
        "L\t47\t+\t49\t+\t0M\n"
        "S\t48\t2{[313, 316)[319, 331)}\tFC:i:0\t\tRC:i:0\n"
        "L\t48\t+\t49\t+\t0M\n"
        "S\t49\t1{[322, 337)}\tFC:i:0\t\tRC:i:0\n"
        "L\t49\t+\t50\t+\t0M\n"
        "S\t50\t1{[336, 351)}\tFC:i:0\t\tRC:i:0\n"
        "L\t50\t+\t51\t+\t0M\n"
        "S\t51\t1{[343, 358)}\tFC:i:0\t\tRC:i:0\n"
        "L\t51\t+\t52\t+\t0M\n"
        "L\t51\t+\t53\t+\t0M\n"
        "S\t52\t3{[357, 360)[364, 365)[374, 385)}\tFC:i:0\t\tRC:i:0\n"
        "L\t52\t+\t54\t+\t0M\n"
        "S\t53\t3{[347, 360)[369, 370)[374, 375)}\tFC:i:0\t\tRC:i:0\n"
        "L\t53\t+\t55\t+\t0M\n"
        "S\t54\t1{[376, 391)}\tFC:i:0\t\tRC:i:0\n"
        "L\t54\t+\t56\t+\t0M\n"
        "S\t55\t3{[353, 360)[369, 370)[374, 381)}\tFC:i:0\t\tRC:i:0\n"
        "L\t55\t+\t54\t+\t0M\n"
        "S\t56\t1{[379, 394)}\tFC:i:0\t\tRC:i:0\n"
        "L\t56\t+\t57\t+\t0M\n"
        "S\t57\t1{[391, 406)}\tFC:i:0\t\tRC:i:0\n"
        "L\t57\t+\t58\t+\t0M\n"
        "S\t58\t1{[399, 414)}\tFC:i:0\t\tRC:i:0\n"
        "L\t58\t+\t59\t+\t0M\n"
        "S\t59\t1{[407, 422)}\tFC:i:0\t\tRC:i:0\n"
        "L\t59\t+\t60\t+\t0M\n"
        "S\t60\t1{[410, 425)}\tFC:i:0\t\tRC:i:0\n"
        "L\t60\t+\t61\t+\t0M\n"
        "S\t61\t1{[415, 430)}\tFC:i:0\t\tRC:i:0\n"
        "L\t61\t+\t62\t+\t0M\n"
        "S\t62\t1{[426, 441)}\tFC:i:0\t\tRC:i:0\n"
        "L\t62\t+\t63\t+\t0M\n"
        "S\t63\t1{[433, 448)}\tFC:i:0\t\tRC:i:0\n"
        "L\t63\t+\t64\t+\t0M\n"
        "S\t64\t1{[434, 449)}\tFC:i:0\t\tRC:i:0\n"
        "L\t64\t+\t65\t+\t0M\n"
        "S\t65\t1{[440, 455)}\tFC:i:0\t\tRC:i:0\n"
        "L\t65\t+\t66\t+\t0M\n"
        "S\t66\t1{[452, 467)}\tFC:i:0\t\tRC:i:0\n"
        "L\t66\t+\t67\t+\t0M\n"
        "S\t67\t1{[465, 480)}\tFC:i:0\t\tRC:i:0\n"
        "L\t67\t+\t68\t+\t0M\n"
        "S\t68\t1{[472, 487)}\tFC:i:0\t\tRC:i:0\n"
        "L\t68\t+\t69\t+\t0M\n"
        "S\t69\t1{[476, 491)}\tFC:i:0\t\tRC:i:0\n"
        "L\t69\t+\t70\t+\t0M\n"
        "S\t70\t1{[491, 491)}\tFC:i:0\t\tRC:i:0\n");
    build_load_and_getters___check_prg(loaded_index, 2, "empty",
        "CACACACACAGAGAGAGAGAGAGAGAGAGATATATATATA",
        "H\tVN:Z:1.0\tbn:Z:--linear --singlearr\n"
        "S\t0\t1{[0, 0)}\tFC:i:0\t\tRC:i:0\n"
        "L\t0\t+\t1\t+\t0M\n"
        "S\t1\t1{[7, 22)}\tFC:i:0\t\tRC:i:0\n"
        "L\t1\t+\t2\t+\t0M\n"
        "S\t2\t1{[9, 24)}\tFC:i:0\t\tRC:i:0\n"
        "L\t2\t+\t3\t+\t0M\n"
        "S\t3\t1{[11, 26)}\tFC:i:0\t\tRC:i:0\n"
        "L\t3\t+\t4\t+\t0M\n"
        "S\t4\t1{[13, 28)}\tFC:i:0\t\tRC:i:0\n"
        "L\t4\t+\t5\t+\t0M\n"
        "S\t5\t1{[15, 30)}\tFC:i:0\t\tRC:i:0\n"
        "L\t5\t+\t6\t+\t0M\n"
        "S\t6\t1{[40, 40)}\tFC:i:0\t\tRC:i:0\n");
}

TEST_F(IndexTest___Fixture, lazy_load)
{
    const auto index_filepath
        = TEST_CASE_DIR + "/sample_example/pangenome.prg.fa.panidx.zip";
    const auto prg_filepath = TEST_CASE_DIR + "/sample_example/pangenome.prg.fa";

    Index::build_index_on_disk(14, 15, prg_filepath, index_filepath);
    Index loaded_index = Index::load(index_filepath);

    // load one PRG and check
    EXPECT_EQ(loaded_index.get_loaded_prgs().size(), 0);
    loaded_index.get_prg_given_id(1);
    EXPECT_EQ(loaded_index.get_loaded_prgs().size(), 1);

    // load another PRG
    loaded_index.get_prg_given_id(0);
    EXPECT_EQ(loaded_index.get_loaded_prgs().size(), 2);

    // load repeated PRGs and check we did not load anything further
    loaded_index.get_prg_given_id(0);
    loaded_index.get_prg_given_id(1);
    EXPECT_EQ(loaded_index.get_loaded_prgs().size(), 2);

    // load last PRG
    loaded_index.get_prg_given_id(2);
    EXPECT_EQ(loaded_index.get_loaded_prgs().size(), 3);
}

TEST_F(IndexTest___Fixture, save_load_minhash)
{
    pandora::KmerHash hash;
    deque<Interval> d = { Interval(3, 5), Interval(9, 12) };
    prg::Path p;
    p.initialize(d);
    pair<uint64_t, uint64_t> kh1 = hash.kmerhash("ACGTA", 5);
    index->add_record(min(kh1.first, kh1.second), 1, p, 0, 0);
    pair<uint64_t, uint64_t> kh2 = hash.kmerhash("ACTGA", 5);
    index->add_record(min(kh2.first, kh2.second), 2, p, 0, 0);
    index->add_record(min(kh1.first, kh1.second), 4, p, 0, 0);

    {
        ZipFileWriter index_archive("save_load_minhash.panidx.zip");
        index->save_minhash(index_archive);
    }

    auto zip_file = std::make_shared<ZipFileReader>("save_load_minhash.panidx.zip");
    auto index_2 = std::make_shared<IndexMock>(1, 5, std::vector<std::string>(),
        std::vector<uint32_t>(), std::vector<uint32_t>(), zip_file);
    index_2->load_minhash();

    EXPECT_EQ(index->minhash.size(), index_2->minhash.size());
    EXPECT_EQ(index->minhash[min(kh1.first, kh1.second)]->size(),
        index_2->minhash[min(kh1.first, kh1.second)]->size());
    EXPECT_EQ(index->minhash[min(kh2.first, kh2.second)]->size(),
        index_2->minhash[min(kh2.first, kh2.second)]->size());
    EXPECT_EQ(index->minhash[min(kh1.first, kh1.second)]->at(0),
        index_2->minhash[min(kh1.first, kh1.second)]->at(0));
    EXPECT_EQ(index->minhash[min(kh1.first, kh1.second)]->at(1),
        index_2->minhash[min(kh1.first, kh1.second)]->at(1));
    EXPECT_EQ(index->minhash[min(kh2.first, kh2.second)]->at(0),
        index_2->minhash[min(kh2.first, kh2.second)]->at(0));
}
