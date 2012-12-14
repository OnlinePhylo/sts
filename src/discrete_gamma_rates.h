/// \file sts/particle/detail/discrete_gamma.hpp
/// \brief A discrete gamma distribution

#ifndef STS_PARTICLE_DISCRETE_GAMMA_RATES_HPP
#define STS_PARTICLE_DISCRETE_GAMMA_RATES_HPP

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <vector>
#include <iostream>
#include "rates.h"

namespace sts
{
namespace particle
{

class Discrete_gamma_rates : public Rates
{
public:
    Discrete_gamma_rates( int categories, double alpha, double beta );
    const double* rates() const;
    double operator[]( int i ) const;
    int count() const;
    double alpha() const;
    double beta() const;
    bool dirty() const;
    void set_alpha(double alpha);
    void set_beta(double beta);
private:
    void initialize();
    double a;
    double b;
    std::vector<double> discrete;
    bool is_dirty;
};

Discrete_gamma_rates::Discrete_gamma_rates( int categories, double alpha, double beta ) :
    a(alpha),
    b(beta),
    discrete(categories),
    is_dirty(false)
{
  initialize();
}

void Discrete_gamma_rates::initialize()
{
    double sep = 1.0/(double)discrete.size();
    unsigned i=0;
    double cumulant = 0;
    for(double v = sep/2.0; v<1.0; v+=sep){
        discrete[i] = gsl_cdf_gamma_Q(v, a, b);
	cumulant += discrete[i];
	i++;
    }
    // normalize so mean of rates is 1
    for(i=0; i<discrete.size(); i++){ 
      discrete[i] /= cumulant;
      discrete[i] *= discrete.size();
    }
    is_dirty=false;
}

const double* Discrete_gamma_rates::rates() const
{
    return discrete.data();
}

double Discrete_gamma_rates::operator[]( int i ) const
{
    return discrete[i];
}

int Discrete_gamma_rates::count() const
{
    return discrete.size();
}

double Discrete_gamma_rates::alpha() const
{
    return a;
}

double Discrete_gamma_rates::beta() const
{
    return b;
}

bool Discrete_gamma_rates::dirty() const
{
    return is_dirty;
}

void Discrete_gamma_rates::set_alpha(double alpha)
{
    is_dirty = is_dirty || this->a != alpha;
    this->a = alpha;
}

void Discrete_gamma_rates::set_beta(double beta)
{
    is_dirty = is_dirty || this->b != beta;
    this->b = beta;
}


} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_DISCRETE_GAMMA_RATES_HPP