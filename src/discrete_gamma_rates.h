/// \file sts/particle/detail/discrete_gamma.hpp
/// \brief A discrete gamma distribution

#ifndef STS_PARTICLE_DISCRETE_GAMMA_RATES_HPP
#define STS_PARTICLE_DISCRETE_GAMMA_RATES_HPP

#include <gsl/gsl_randist.h>
#include <vector>

namespace sts
{
namespace particle
{

class discrete_gamma_rates
{
public:
    discrete_gamma_rates( int categories, double alpha, double beta );
    double operator[]( int i ) const;
    int count() const;
    double alpha() const;
    double beta() const;
    bool dirty() const;
    void set_alpha(double alpha);
    void set_beta(double beta);
private:
    double alpha;
    double beta;
    std::vector<double> discrete;
    bool dirty;
};

discrete_gamma_rates::discrete_gamma_rates( int categories, double alpha, double beta ) :
    alpha(alpha),
    beta(beta),
    discrete(categories),
    dirty(true)
{
    double sep = 1.0/(double)categories;
    int i=0;
    for(double v = sep/2.0; v+=sep; v<1.0){
        discrete[i] = gsl_cdf_gamma_Q(v, alpha, beta);
    }
}

double discrete_gamma_rates::operator[]( int i ) const
{
    return discrete[i];
}

int discrete_gamma_rates::count() const
{
    return discrete.size();
}

double discrete_gamma_rates::alpha() const
{
    return alpha;
}

double discrete_gamma_rates::beta() const
{
    return beta;
}

bool discrete_gamma_rates::dirty() const
{
    return dirty;
}

void discrete_gamma_rates::set_alpha(double alpha)
{
    dirty = dirty || this->alpha != alpha;
    this->alpha = alpha;
}

void discrete_gamma_rates::set_beta(double beta)
{
    dirty = dirty || this->beta != beta;
    this->beta = beta;
}


} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_DISCRETE_GAMMA_RATES_HPP