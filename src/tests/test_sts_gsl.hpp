#ifndef TEST_STS_GSL_HPP
#define TEST_STS_GSL_HPP

#include <cmath>
#include "catch.hpp"
#include "gsl.h"

namespace sts
{
namespace test
{
namespace gsl
{

struct SquareFn {
    double operator()(double d) { return d * d; }
};

TEST_CASE("sts/gsl/minimize", "sts::gsl::minimize is a wrapper for GSL minimization functions")
{
    SquareFn instance;
    REQUIRE(std::abs(sts::gsl::minimize(instance, -100.0, -1000, 1000.0)) < 1e-5);
}

} //gsl
} //test
} //sts

#endif //TEST_STS_GSL_HPP
