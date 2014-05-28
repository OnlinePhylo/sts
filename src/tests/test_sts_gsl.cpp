#include <cmath>
#include "gtest/gtest.h"
#include "gsl.h"


struct SquareFn {
    double operator()(double d) { return d * d; }
};

TEST(sts_gsl, minimize)
{
    SquareFn instance;
    ASSERT_NEAR(sts::gsl::minimize(instance, -100.0, -1000, 1000.0), 0, 1e-5);
}
