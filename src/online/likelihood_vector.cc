#include "likelihood_vector.h"
#include <cassert>
#include <cmath>

namespace sts { namespace online {

Likelihood_vector::Likelihood_vector(const size_t n_rates,
                                     const size_t n_sites,
                                     const size_t n_states) :
    n_rates_(n_rates),
    n_sites_(n_sites),
    n_states_(n_states),
    v(n_rates*n_sites*n_states)
{}

Likelihood_vector::Likelihood_vector(Likelihood_vector&& other) :
    n_rates_(other.n_rates_),
    n_sites_(other.n_sites_),
    n_states_(other.n_states_),
    v(std::move(other.v))
{}

Likelihood_vector& Likelihood_vector::operator=(Likelihood_vector&& other)
{
    n_rates_ = other.n_rates_;
    n_sites_ = other.n_sites_;
    n_states_ = other.n_states_;
    v = std::move(other.v);
    return *this;
}

double Likelihood_vector::log_dot(const Likelihood_vector& other) const
{
    std::vector<double> site_likes(n_sites(), 0.0);

    for(size_t rate = 0; rate < n_rates(); rate++) {
        for(size_t site = 0; site < n_sites(); site++) {
            for(size_t state = 0; state < n_states(); state++) {
                site_likes[site] += operator()(rate, site, state) * other(rate, site, state);
            }
        }
    }

    double result = 0.0;
    for(const double d : site_likes) {
        result += std::log(d);
    }

    return result - (n_sites() * std::log(n_rates()));
}

inline size_t Likelihood_vector::index(const size_t rate, const size_t site, const size_t state) const
{
    assert(rate < n_rates());
    assert(site < n_sites());
    assert(state < n_states());
    return (rate * n_sites() * n_states()) + (site * n_states()) + state;
}

}}
