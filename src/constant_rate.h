/// \file sts/particle/detail/constant_rate.hpp
/// \brief A single constant rate parameter

#ifndef STS_PARTICLE_CONSTANT_RATE_HPP
#define STS_PARTICLE_CONSTANT_RATE_HPP

namespace sts
{
namespace particle
{

class constant_rate
{
public:
    constant_rate( double rate );
    double operator[]( int i ) const;
    int count() const;
    bool dirty() const;
private:
    double rate;
};

constant_rate::constant_rate( double rate ) :
    rate(rate),
{}

double constant_rate::operator[]( int i ) const
{
    return rate;
}

int constant_rate::count() const
{
    return 1;
}

bool constant_rate::dirty() const
{
    return false;
}


} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_CONSTANT_RATE_HPP