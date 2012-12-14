/// \file sts/particle/detail/constant_rate.hpp
/// \brief A single constant rate parameter

#ifndef STS_PARTICLE_CONSTANT_RATE_HPP
#define STS_PARTICLE_CONSTANT_RATE_HPP
#include "rates.h"

namespace sts
{
namespace particle
{

class Constant_rate : public Rates
{
public:
    Constant_rate( double rate );
    const double* rates() const;
    double operator[]( int i ) const;
    int count() const;
    bool dirty() const;
private:
    double rate;
};

Constant_rate::Constant_rate( double rate ) :
    rate(rate)
{}

const double* Constant_rate::rates() const
{
  return &rate;
}

double Constant_rate::operator[]( int i ) const
{
    return rate;
}

int Constant_rate::count() const
{
    return 1;
}

bool Constant_rate::dirty() const
{
    return false;
}


} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_CONSTANT_RATE_HPP