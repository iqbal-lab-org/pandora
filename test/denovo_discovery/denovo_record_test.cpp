#include "denovo_discovery/denovo_record.h"
#include "gtest/gtest.h"


TEST(DenovoVariantRecord, to_string) {
    DenovoVariantRecord denovo_variant_record(23, "AAAAA", "CCC");
    std::string actual = denovo_variant_record.to_string();
    std::string expected {"23\tAAAAA\tCCC"};
    EXPECT_EQ(actual, expected);
}


template <typename T,
    template <typename ELEM_TYPE, typename = std::allocator<ELEM_TYPE>> class CONT_TYPE>
bool equal_containers(const CONT_TYPE<T>& lhs, const CONT_TYPE<T>& rhs)
{
    return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

class DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture : public ::testing::Test {
public:
    DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture()
        : ref_as_str("AAGATATGTTCCGTTATGCGCAGCCCACA") { }
    std::string ref_as_str;
};

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
// Trivial tests
TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, no_variant)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, ref_as_str);
    std::vector<DenovoVariantRecord> expected;
    EXPECT_TRUE(equal_containers(actual, expected));
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
// SNPs tests
TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, SNP_in_the_start)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "CAGATATGTTCCGTTATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(1, "A", "C")
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, SNP_in_the_middle)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTGATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(15, "T", "G")
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, SNP_in_the_end)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTATGCGCAGCCCACT");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(29, "A", "T")
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, MSNP_in_the_start)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "CGTGTATGTTCCGTTATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(1, "AAGA", "CGTG")
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, MSNP_in_the_middle)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTGGCGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(15, "TAT", "GGC")
    };
    EXPECT_EQ(actual, expected);
}


TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, MSNP_in_the_end)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTATGCGCAGCCCAGT");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(28, "CA", "GT")
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, two_distant_SNPs)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATCTTCCGTTATGCGCAGGCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(8, "G", "C"),
        DenovoVariantRecord(24, "C", "G")
    };
    EXPECT_EQ(actual, expected);
}


TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, three_distant_SNPs)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATCTTCCGTTATGCGCAGGCCACT");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(8, "G", "C"),
        DenovoVariantRecord(24, "C", "G"),
        DenovoVariantRecord(29, "A", "T")
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, two_SNPs_and_2_MSNPs)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AACATATGAACCGTTTTGCGTTTTTCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(3, "G", "C"),
        DenovoVariantRecord(9, "TT", "AA"),
        DenovoVariantRecord(16, "A", "T"),
        DenovoVariantRecord(21, "CAGCC", "TTTTT")
    };
    EXPECT_EQ(actual, expected);
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
// deletion tests
TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, deletion_in_the_start)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AGATATGTTCCGTTATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(1, "A", ""),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, deletion_in_the_middle)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTTGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(16, "A", ""),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, deletion_in_the_end)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTATGCGCAGCCCAC");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(29, "A", ""),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, multiple_deletion_in_the_start)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "TATGTTCCGTTATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(1, "AAGA", ""),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, multiple_deletion_in_the_middle)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(14, "TTAT", ""),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, multiple_deletion_in_the_end)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTATGCGCAGCCCA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(28, "CA", ""),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, two_distant_deletions)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATTTCCGTTATGCGCACCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(8, "G", ""),
        DenovoVariantRecord(23, "G", "")
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, three_distant_deletions)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATTTCCGTTATGCGCAGCCAC");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(8, "G", ""),
        DenovoVariantRecord(24, "C", ""),
        DenovoVariantRecord(29, "A", "")
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, two_deletions_and_2_multiple_deletions)
{
    std::vector<DenovoVariantRecord> actual = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAATATGCCGTTTGCGCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(3, "G", ""),
        DenovoVariantRecord(9, "TT", ""),
        DenovoVariantRecord(16, "A", ""),
        DenovoVariantRecord(21, "CAGCC", "")
    };
    EXPECT_EQ(actual, expected);
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//// insertion tests
TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, insertion_in_the_start)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "TAAGATATGTTCCGTTATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(1, "", "T"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, insertion_in_the_middle)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTCATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(16, "", "C"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, insertion_in_the_end)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTATGCGCAGCCCACAG");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(30, "", "G"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, multiple_insertion_in_the_start)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "CGTGAAGATATGTTCCGTTATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(1, "", "CGTG"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, multiple_insertion_in_the_middle)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTGGCATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(16, "", "GGC"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, multiple_insertion_in_the_end)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTATGCGCAGCCCACAGT");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(30, "", "GT"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, two_distant_insertions)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATCGTTCCGTTATGCGCAGACCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(8, "", "C"),
        DenovoVariantRecord(24, "", "A"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, three_distant_insertions)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATGATGTTCCGTTACTGCGCAGCCCTACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(6, "", "G"),
        DenovoVariantRecord(17, "", "C"),
        DenovoVariantRecord(27, "", "T"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(
    DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, two_insertions_and_2_multiple_insertions)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AATGATATGAATTCCGTCTATGCGCAGCCTTTTTCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(3, "", "T"),
        DenovoVariantRecord(9, "", "AA"),
        DenovoVariantRecord(15, "", "C"),
        DenovoVariantRecord(26, "", "TTTTT"),
    };
    EXPECT_EQ(actual, expected);
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
// Mixed variants tests
TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, SNP_and_insertion_in_the_start)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "CTTAGATATGTTCCGTTATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(1, "A", "CTT"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, SNP_and_deletion_in_the_start)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "CATATGTTCCGTTATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(1, "AAG", "C"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, SNP_and_insertion_in_the_middle)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTGCCATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(15, "T", "GCC"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, SNP_and_deletion_in_the_middle)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTCGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(15, "TAT", "C"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, SNP_and_insertion_in_the_end)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTATGCGCAGCCCACGTGC");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(29, "A", "GTGC"),
    };
    EXPECT_EQ(actual, expected);
}


TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, SNP_and_deletion_in_the_end)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTATGCGCAGCCG");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(26, "CACA", "G"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, MSNP_and_insertion_in_the_start)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "CCTTTATATGTTCCGTTATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(1, "AAG", "CCTTT"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, MSNP_and_deletion_in_the_start)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "CCATGTTCCGTTATGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(1, "AAGAT", "CC"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, MSNP_and_insertion_in_the_middle)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTGCCCCGGCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(15, "TAT", "GCCCCG"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, MSNP_and_deletion_in_the_middle)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTCCCGCAGCCCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(15, "TATG", "CC"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, MSNP_and_insertion_in_the_end)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTATGCGCAGCCCTTTGTGT");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(27, "ACA", "TTTGTGT"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, MSNP_and_deletion_in_the_end)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGATATGTTCCGTTATGCGCATTT");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(23, "GCCCACA", "TTT"),
    };
    EXPECT_EQ(actual, expected);
}

TEST_F(DenovoVariantRecord___get_variants_from_pair_of_sequences___Fixture, two_distant_MultiVariants)
{
    std::vector<DenovoVariantRecord> actual
        = DenovoVariantRecord::get_variants_from_pair_of_sequences(ref_as_str, "AAGACCCCCGTTATGCGCATTTTTTTTTTCACA");
    std::vector<DenovoVariantRecord> expected {
        DenovoVariantRecord(5, "TATGTT", "CCC"),
        DenovoVariantRecord(23, "GCC", "TTTTTTTTTT"),
    };
    EXPECT_EQ(actual, expected);
}
