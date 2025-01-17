#include "Definitions.hpp"
#include "Rational.hpp"

#include <gtest/gtest.h>
#include <iostream>

class RationalTest: public ::testing::Test {};

TEST_F(RationalTest, normalize) {
    const zx::PiRational r(-33, 16);
    EXPECT_EQ(r, zx::PiRational(-1, 16));
}

TEST_F(RationalTest, from_double) {
    const zx::PiRational r(-zx::PI / 8);
    EXPECT_EQ(r, zx::PiRational(-1, 8));
}

TEST_F(RationalTest, from_double_2) {
    const zx::PiRational r(-3 * zx::PI / 4);
    EXPECT_EQ(r, zx::PiRational(-3, 4));
}

TEST_F(RationalTest, from_double_3) {
    const zx::PiRational r(-7 * zx::PI / 8);
    EXPECT_EQ(r, zx::PiRational(-7, 8));
}

TEST_F(RationalTest, from_double_4) {
    const zx::PiRational r(-1 * zx::PI / 32);
    EXPECT_EQ(r, zx::PiRational(-1, 32));
}

TEST_F(RationalTest, from_double_5) {
    const zx::PiRational r(5000 * zx::PI + zx::PI / 4);
    EXPECT_EQ(r, zx::PiRational(1, 4));
}

TEST_F(RationalTest, from_double_6) {
    const zx::PiRational r(-5000 * zx::PI + 5 * zx::PI / 4);
    EXPECT_EQ(r, zx::PiRational(-3, 4));
}

// TEST_F(RationalTest, from_double_7) {
//     zx::PiRational r(0.1);
//     EXPECT_EQ(r, zx::PiRational(-3, 4));
// }

TEST_F(RationalTest, add) {
    const zx::PiRational r0(1, 8);
    const zx::PiRational r1(7, 8);
    const auto           r = r0 + r1;

    EXPECT_EQ(r, 1);
}

TEST_F(RationalTest, add_2) {
    const zx::PiRational r0(9, 8);
    const zx::PiRational r1(7, 8);
    const auto           r = r0 + r1;

    EXPECT_EQ(r, 0);
}

TEST_F(RationalTest, sub) {
    const zx::PiRational r0(9, 8);
    const zx::PiRational r1(-7, 8);
    const auto           r = r0 - r1;

    EXPECT_EQ(r, 0);
}

TEST_F(RationalTest, sub_2) {
    const zx::PiRational r0(-1, 2);
    const zx::PiRational r1(1, 2);
    const auto           r = r0 - r1;

    EXPECT_EQ(r, 1);
}

TEST_F(RationalTest, mul) {
    const zx::PiRational r0(1, 8);
    const zx::PiRational r1(1, 2);
    const auto           r = r0 * r1;

    EXPECT_EQ(r, zx::PiRational(1, 16));
}

TEST_F(RationalTest, mul_2) {
    const zx::PiRational r0(1, 8);
    const zx::PiRational r1(0, 1);
    const auto           r = r0 * r1;

    EXPECT_EQ(r, 0);
}

TEST_F(RationalTest, div) {
    const zx::PiRational r0(1, 2);
    const zx::PiRational r1(1, 2);
    const auto           r = r0 / r1;

    EXPECT_EQ(r, 1);
}

TEST_F(RationalTest, approximateDivPi) {
    const zx::PiRational r(1'000'000'000'000'000 - 1, 1'000'000'000'000'000);
    EXPECT_TRUE(r.isCloseDivPi(1, 1e-9));
}

TEST_F(RationalTest, approximate) {
    const zx::PiRational r(1, 1);
    EXPECT_TRUE(r.isClose(3.14159, 1e-5));
}
