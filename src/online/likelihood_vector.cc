#include "likelihood_vector.h"
#include <algorithm>
#include <cassert>
#include <cmath>

namespace sts { namespace online {

LikelihoodVector::LikelihoodVector(const size_t nRates,
                                   const size_t nSites,
                                   const size_t nStates) :
    nRates_(nRates),
    nSites_(nSites),
    nStates_(nStates),
    v(nRates*nSites*nStates)
{}

LikelihoodVector::LikelihoodVector(LikelihoodVector&& other) :
    nRates_(other.nRates_),
    nSites_(other.nSites_),
    nStates_(other.nStates_),
    v(std::move(other.v))
{}

LikelihoodVector& LikelihoodVector::operator=(LikelihoodVector&& other)
{
    nRates_ = other.nRates_;
    nSites_ = other.nSites_;
    nStates_ = other.nStates_;
    v = std::move(other.v);
    return *this;
}

double LikelihoodVector::logDot(const LikelihoodVector& other) const
{
    std::vector<double> site_likes(nSites(), 0.0);
    assert(other.nRates() == nRates());
    assert(other.nSites() == nSites());
    assert(other.nStates() == nStates());

    for(size_t rate = 0; rate < nRates(); rate++) {
        for(size_t site = 0; site < nSites(); site++) {
            size_t idx = index(rate, site, 0);
            site_likes[site] += std::inner_product(v.begin() + idx,
                                                   v.begin() + idx + nStates(),
                                                   other.v.begin() + idx,
                                                   0.0);
        }
    }

    double result = 0.0;
    for(const double d : site_likes) {
        result += std::log(d);
    }

    // exp(result) / n_rates^n_sites
    // **Assumes equiprobable rates.**
    return result - (std::log(nRates()) * nSites());
}

double LikelihoodVector::logDot(const LikelihoodVector& other, const std::vector<double>& weights) const
{
    assert(weights.size() == nRates());
    assert(other.nRates() == nRates());
    assert(other.nSites() == nSites());
    assert(other.nStates() == nStates());
    std::vector<double> siteLikes(nSites(), 0.0);

    for(size_t rate = 0; rate < nRates(); rate++) {
        for(size_t site = 0; site < nSites(); site++) {
            size_t idx = index(rate, site, 0);
            const double p = std::inner_product(v.begin() + idx,
                                                v.begin() + idx + nStates(),
                                                other.v.begin() + idx,
                                                0.0);
            siteLikes[site] += p * weights[rate];
        }
    }

    double result = 0.0;
    for(const double d : siteLikes) {
        result += std::log(d);
    }

    return result;
}

inline size_t LikelihoodVector::index(const size_t rate, const size_t site, const size_t state) const
{
    assert(rate < nRates());
    assert(site < nSites());
    assert(state < nStates());
    return (rate * nSites() * nStates()) + (site * nStates()) + state;
}

}}
