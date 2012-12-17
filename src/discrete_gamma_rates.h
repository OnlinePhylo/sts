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



} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_DISCRETE_GAMMA_RATES_HPP