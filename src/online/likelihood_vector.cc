#include "likelihood_vector.h"
#include <algorithm>
#include <cassert>
#include <cmath>

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>

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

inline size_t LikelihoodVector::index(const size_t rate, const size_t site, const size_t state) const
{
    assert(rate < nRates());
    assert(site < nSites());
    assert(state < nStates());
    return (rate * nSites() * nStates()) + (site * nStates()) + state;
}

}}
