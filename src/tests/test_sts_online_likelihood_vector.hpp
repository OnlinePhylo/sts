#ifndef STS_TEST_ONLINE_LIKELIHOOD_VECTOR_H
#define STS_TEST_ONLINE_LIKELIHOOD_VECTOR_H
#include "catch.hpp"

#include "likelihood_vector.h"
#include <algorithm>
#include <cmath>

namespace sts { namespace test { namespace online {


TEST_CASE("sts/online/likelihood_vector/self_product", "Product with self")
{
    using sts::online::LikelihoodVector;
    LikelihoodVector v(2, 4, 6);
    std::fill(v.get().begin(), v.get().end(), 1.0);
    CHECK(std::exp(v.logDot(v)) == Approx(1.0));
}

TEST_CASE("sts/online/likelihood_vector/log_like", "Log-likelihood from vector")
{
    using sts::online::LikelihoodVector;
    LikelihoodVector v(2, 4, 6);
    std::fill(v.get().begin(), v.get().end(), 1.0);
    CHECK(std::exp(v.logLikelihood()) == Approx(1.0));
}

}}} // namespaces

#endif
