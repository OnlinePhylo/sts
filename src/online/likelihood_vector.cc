#include "likelihood_vector.h"
#include <algorithm>
#include <cassert>
#include <cmath>

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Phyl/Model/JCnuc.h>
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

double LikelihoodVector::logLikelihood() const
{
    const std::vector<double> freqs(nStates(), 1.0 / nStates());
    const std::vector<double> weights(nRates(), 1.0 / nRates());
    return logLikelihood(freqs, weights);
}

double LikelihoodVector::logLikelihood(const bpp::SubstitutionModel& model, const bpp::DiscreteDistribution& rate_dist) const
{
    assert(rate_dist.getNumberOfCategories() == nRates());
    assert(model.getNumberOfStates() == nStates());

    const std::vector<double> weights = rate_dist.getProbabilities();
    const std::vector<double>& freqs = model.getFrequencies();
    return logLikelihood(freqs, weights);
}

double LikelihoodVector::logLikelihood(const std::vector<double>& freqs, const std::vector<double>& rateWeights) const
{
    assert(rateWeights.size() == nRates());
    assert(freqs.size() == nStates());

    std::vector<double> siteStateLikes(nSites() * nStates(), 0.0); // Integral over rates

    // Integrate over rates for each site/state
    for(size_t rate = 0; rate < nRates(); rate++) {
        for(size_t site = 0; site < nSites(); site++) {
            for(size_t state = 0; state < nStates(); state++) {
                const size_t vIdx = index(rate, site, state);
                const size_t ssIdx = nStates() * site + state;
                siteStateLikes[ssIdx] += v[vIdx] * rateWeights[rate];
            }
        }
    }

    // Integrate over states for each site
    std::vector<double> siteLikes(nSites(), 0.0); // Likelihood per site
    for(size_t site = 0; site < nSites(); site++) {
        double result = 0.0;
        for(size_t state = 0; state < nStates(); state++) {
            const size_t ssIdx = nStates() * site + state;
            result += freqs[state] * siteStateLikes[ssIdx];
        }
        siteLikes[site] = std::log(result);
    }

    // Final LL = sum of log likelihoods for each site
    return std::accumulate(siteLikes.begin(), siteLikes.end(), 0.0);
}

double LikelihoodVector::logDot(const LikelihoodVector& other) const
{
    const std::vector<double> freqs(nStates(), 1.0 / nStates());
    const std::vector<double> weights(nRates(), 1.0 / nRates());
    return logDot(other, freqs, weights);
}

double LikelihoodVector::logDot(const LikelihoodVector& other, const bpp::SubstitutionModel& model,
                                const bpp::DiscreteDistribution& rate_dist) const
{
    const std::vector<double>& freqs = model.getFrequencies();
    const std::vector<double> weights = rate_dist.getProbabilities();
    return logDot(other, freqs, weights);
}

double LikelihoodVector::logDot(const LikelihoodVector& other, const std::vector<double>& freqs, const std::vector<double>& rateWeights) const
{
    assert(freqs.size() == nStates());
    assert(rateWeights.size() == nRates());
    assert(other.nRates() == nRates());
    assert(other.nSites() == nSites());
    assert(other.nStates() == nStates());

    std::vector<double> tmp(nSites() * nStates(), 0.0);
    std::vector<double> siteLikes(nSites(), 0.0);

    for(size_t rate = 0; rate < nRates(); rate++) {
        for(size_t site = 0; site < nSites(); site++) {
            for(size_t state = 0; state < nStates(); state++) {
                const size_t v_idx = index(rate, site, state);
                const size_t tmp_idx = nStates() * site + state;
                tmp[tmp_idx] += other.v[v_idx] * v[v_idx] * rateWeights[rate];
            }
        }
    }

    for(size_t site = 0; site < nSites(); site++) {
        double result = 0.0;
        for(size_t state = 0; state < nStates(); state++) {
            const size_t tmp_idx = nStates() * site + state;
            result += freqs[state] * tmp[tmp_idx];
        }
        siteLikes[site] = std::log(result);
    }

    return std::accumulate(siteLikes.begin(), siteLikes.end(), 0.0);
}

inline size_t LikelihoodVector::index(const size_t rate, const size_t site, const size_t state) const
{
    assert(rate < nRates());
    assert(site < nSites());
    assert(state < nStates());
    return (rate * nSites() * nStates()) + (site * nStates()) + state;
}

}}
