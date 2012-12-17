#include "discrete_gamma_rates.h"

namespace sts
{
namespace particle
{

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
        discrete[i] = gsl_cdf_gamma_Qinv(v, a, b);
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
    initialize();
}

void Discrete_gamma_rates::set_beta(double beta)
{
    is_dirty = is_dirty || this->b != beta;
    this->b = beta;
    initialize();
}

} // namespace particle
} // namespace sts

