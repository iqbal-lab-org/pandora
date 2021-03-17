#include "gtest/gtest.h"

#include "Maths.h"

class MathsTest___all_methods___Fixture : public ::testing::Test {
protected:
    void SetUp() override
    {
        four_sized_vector_int.push_back(20);
        four_sized_vector_int.push_back(8);
        four_sized_vector_int.push_back(7);
        four_sized_vector_int.push_back(10);
        one_sized_vector_int.insert(one_sized_vector_int.end(),
            four_sized_vector_int.begin(), four_sized_vector_int.begin() + 1);
        two_sized_vector_int.insert(two_sized_vector_int.end(),
            four_sized_vector_int.begin(), four_sized_vector_int.begin() + 2);
        three_sized_vector_int.insert(three_sized_vector_int.end(),
            four_sized_vector_int.begin(), four_sized_vector_int.begin() + 3);

        four_sized_vector_double.push_back(20.2);
        four_sized_vector_double.push_back(8.8);
        four_sized_vector_double.push_back(7.7);
        four_sized_vector_double.push_back(10.10);
        one_sized_vector_double.insert(one_sized_vector_double.end(),
            four_sized_vector_double.begin(), four_sized_vector_double.begin() + 1);
        two_sized_vector_double.insert(two_sized_vector_double.end(),
            four_sized_vector_double.begin(), four_sized_vector_double.begin() + 2);
        three_sized_vector_double.insert(three_sized_vector_double.end(),
            four_sized_vector_double.begin(), four_sized_vector_double.begin() + 3);
    }

    void TearDown() override { }
    std::vector<uint32_t> empty_vector_int;
    std::vector<uint32_t> one_sized_vector_int;
    std::vector<uint32_t> two_sized_vector_int;
    std::vector<uint32_t> three_sized_vector_int;
    std::vector<uint32_t> four_sized_vector_int;

    std::vector<double> empty_vector_double;
    std::vector<double> one_sized_vector_double;
    std::vector<double> two_sized_vector_double;
    std::vector<double> three_sized_vector_double;
    std::vector<double> four_sized_vector_double;
};

TEST_F(MathsTest___all_methods___Fixture, sum___empty_vector_int___returns_zero)
{
    uint32_t actual = Maths::sum(empty_vector_int.begin(), empty_vector_int.end());

    uint32_t expected = 0;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, sum___one_sized_vector_int)
{
    uint32_t actual
        = Maths::sum(one_sized_vector_int.begin(), one_sized_vector_int.end());

    uint32_t expected = 20;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, sum___two_sized_vector_int)
{
    uint32_t actual
        = Maths::sum(two_sized_vector_int.begin(), two_sized_vector_int.end());

    uint32_t expected = 28;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, sum___empty_vector_double___returns_zero)
{
    double actual = Maths::sum(empty_vector_double.begin(), empty_vector_double.end());

    double expected = 0.0;
    EXPECT_TRUE(Maths::equals(actual, expected));
}

TEST_F(MathsTest___all_methods___Fixture, sum___one_sized_vector_double)
{
    double actual
        = Maths::sum(one_sized_vector_double.begin(), one_sized_vector_double.end());

    double expected = 20.2;
    EXPECT_TRUE(Maths::equals(actual, expected));
}

TEST_F(MathsTest___all_methods___Fixture, sum___two_sized_vector_double)
{
    double actual
        = Maths::sum(two_sized_vector_double.begin(), two_sized_vector_double.end());

    double expected = 29.0;
    EXPECT_TRUE(Maths::equals(actual, expected));
}

TEST_F(MathsTest___all_methods___Fixture, mean___empty_vector_int___returns_zero)
{
    uint32_t actual = Maths::mean(empty_vector_int.begin(), empty_vector_int.end());

    uint32_t expected = 0;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, mean___one_sized_vector_int)
{
    uint32_t actual
        = Maths::mean(one_sized_vector_int.begin(), one_sized_vector_int.end());

    uint32_t expected = 20;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, mean___two_sized_vector_int)
{
    uint32_t actual
        = Maths::mean(two_sized_vector_int.begin(), two_sized_vector_int.end());

    uint32_t expected = 14;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, mean___four_sized_vector_int)
{
    uint32_t actual
        = Maths::mean(four_sized_vector_int.begin(), four_sized_vector_int.end());

    uint32_t expected = 11;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, mean___empty_vector_double___returns_zero)
{
    double actual = Maths::mean(empty_vector_double.begin(), empty_vector_double.end());

    double expected = 0.0;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, mean___one_sized_vector_double)
{
    double actual
        = Maths::mean(one_sized_vector_double.begin(), one_sized_vector_double.end());

    double expected = 20.2;
    EXPECT_TRUE(Maths::equals(actual, expected));
}

TEST_F(MathsTest___all_methods___Fixture, mean___two_sized_vector_double)
{
    double actual
        = Maths::mean(two_sized_vector_double.begin(), two_sized_vector_double.end());

    double expected = 14.5;
    EXPECT_TRUE(Maths::equals(actual, expected));
}

TEST_F(MathsTest___all_methods___Fixture, mean___four_sized_vector_double)
{
    double actual
        = Maths::mean(four_sized_vector_double.begin(), four_sized_vector_double.end());

    double expected = 11.7;
    EXPECT_TRUE(Maths::equals(actual, expected));
}

TEST_F(MathsTest___all_methods___Fixture, median___empty_vector_int___returns_zero)
{
    uint32_t actual = Maths::median(empty_vector_int.begin(), empty_vector_int.end());

    uint32_t expected = 0;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, median___one_sized_vector_int)
{
    uint32_t actual
        = Maths::median(one_sized_vector_int.begin(), one_sized_vector_int.end());

    uint32_t expected = 20;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, median___two_sized_vector_int)
{
    uint32_t actual
        = Maths::median(two_sized_vector_int.begin(), two_sized_vector_int.end());

    uint32_t expected = 14;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, median___three_sized_vector_int)
{
    uint32_t actual
        = Maths::median(three_sized_vector_int.begin(), three_sized_vector_int.end());

    uint32_t expected = 8;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, median___four_sized_vector_int)
{
    uint32_t actual
        = Maths::median(four_sized_vector_int.begin(), four_sized_vector_int.end());

    uint32_t expected = 9;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, median___empty_vector_double___returns_zero)
{
    double actual
        = Maths::median(empty_vector_double.begin(), empty_vector_double.end());

    double expected = 0.0;
    EXPECT_EQ(actual, expected);
}

TEST_F(MathsTest___all_methods___Fixture, median___one_sized_vector_double)
{
    double actual
        = Maths::median(one_sized_vector_double.begin(), one_sized_vector_double.end());

    double expected = 20.2;
    EXPECT_TRUE(Maths::equals(actual, expected));
}

TEST_F(MathsTest___all_methods___Fixture, median___two_sized_vector_double)
{
    double actual
        = Maths::median(two_sized_vector_double.begin(), two_sized_vector_double.end());

    double expected = 14.5;
    EXPECT_TRUE(Maths::equals(actual, expected));
}

TEST_F(MathsTest___all_methods___Fixture, median___three_sized_vector_double)
{
    double actual = Maths::median(
        three_sized_vector_double.begin(), three_sized_vector_double.end());

    double expected = 8.8;
    EXPECT_TRUE(Maths::equals(actual, expected));
}

TEST_F(MathsTest___all_methods___Fixture, median___four_sized_vector_double)
{
    double actual = Maths::median(
        four_sized_vector_double.begin(), four_sized_vector_double.end());

    double expected = 9.45;
    EXPECT_TRUE(Maths::equals(actual, expected));
}

TEST_F(MathsTest___all_methods___Fixture, mode___empty_vector_int___returns_zero)
{
    uint32_t actual = Maths::mode(empty_vector_int.begin(), empty_vector_int.end());

    uint32_t expected = 0;
    EXPECT_EQ(actual, expected);
}

TEST(MathsTest, mode___vector_with_only_one_value_int)
{
    std::vector<uint32_t> v { 1, 1, 1, 1 };

    uint32_t actual = Maths::mode(v.begin(), v.end());

    uint32_t expected = 1;
    EXPECT_EQ(actual, expected);
}

TEST(MathsTest, mode___vector_with_two_values_random_order_int)
{
    std::vector<uint32_t> v { 5, 1, 5, 1, 5, 5, 5 };

    uint32_t actual = Maths::mode(v.begin(), v.end());

    uint32_t expected = 5;
    EXPECT_EQ(actual, expected);
}

TEST(MathsTest, mode___vector_with_two_values_another_random_order_int)
{
    std::vector<uint32_t> v { 1, 5, 5, 5, 5, 5, 1, 1 };

    uint32_t actual = Maths::mode(v.begin(), v.end());

    uint32_t expected = 5;
    EXPECT_EQ(actual, expected);
}

TEST(MathsTest,
    mode___vector_with_two_values_random_order_same_amount_return_smallest_int)
{
    std::vector<uint32_t> v { 1, 5, 1, 1, 5, 5, 5, 5, 1, 1 };

    uint32_t actual = Maths::mode(v.begin(), v.end());

    uint32_t expected = 1;
    EXPECT_EQ(actual, expected);
}
