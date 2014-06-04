#include "log_tricks.h"
#include "gtest/gtest.h"

using namespace sts::online;

const double TOLERANCE = 1e-6;

TEST(LogSum, BothLarge) {
    const double x = 0.1;
    const double y = 0.1;
    ASSERT_NEAR(std::log(0.2), logSum(std::log(x), std::log(y)), TOLERANCE);
    ASSERT_NEAR(std::log(0.2), logSum(std::log(y), std::log(x)), TOLERANCE);
}

TEST(LogSum, OneLarge1) {
    const double x = std::log(0.1);
    const double y = -500000;
    ASSERT_NEAR(x, logSum(x, y), TOLERANCE);
    ASSERT_NEAR(x, logSum(y, x), TOLERANCE);
}

TEST(LogSum, OneLarge2) {
    const double x = -10000;
    const double y = -50000;
    ASSERT_NEAR(x, logSum(x, y), TOLERANCE);
    ASSERT_NEAR(x, logSum(y, x), TOLERANCE);
}

TEST(LogSum, OneLargeOneMedium) {
    const double x = std::log(0.1);
    const double y = std::log(TOLERANCE);
    ASSERT_NEAR(std::log(0.1 + TOLERANCE), logSum(x, y), 1e-8);
    ASSERT_NEAR(std::log(0.1 + TOLERANCE), logSum(y, x), 1e-8);
}
