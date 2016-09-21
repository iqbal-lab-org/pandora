#include "gtest/gtest.h"
#include "functions.h"

TEST(IntAddition, Negative) {
	EXPECT_EQ(-5, int_addition(-2, -3)) << "This will be shown in case it fails";
	EXPECT_EQ(-3, int_addition(5, -8));
}

TEST(IntAddition, Positive) {
	EXPECT_EQ(4, int_addition(1, 3));
	EXPECT_EQ(9, int_addition(4, 5));
}

int main(int argc, char **argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
