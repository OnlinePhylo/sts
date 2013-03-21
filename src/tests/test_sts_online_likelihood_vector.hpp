#ifndef STS_TEST_ONLINE_LIKELIHOOD_VECTOR_H
#define STS_TEST_ONLINE_LIKELIHOOD_VECTOR_H
#include "catch.hpp"

#include "likelihood_vector.h"
#include <algorithm>
#include <cmath>

namespace sts { namespace test { namespace online {


TEST_CASE("sts/online/likelihood_vector/self_product", "Product with self")
{
    using sts::online::Likelihood_vector;
    Likelihood_vector v(2, 1, 4);
    std::fill(v.get().begin(), v.get().end(), 1.0);
    REQUIRE(std::exp(v.log_dot(v)) == Approx(4.0));
}


}}} // namespaces

#endif
