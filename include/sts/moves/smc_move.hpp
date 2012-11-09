/// \file smc_move.hpp
/// \brief Smc_move class

#ifndef STS_MOVES_SMC_MOVE_HPP
#define STS_MOVES_SMC_MOVE_HPP

#include "smctc.hh"
#include "sts/likelihood/forest_likelihood.hpp"
#include "sts/particle/state.hpp"

namespace sts
{
namespace moves
{

/// \class Smc_move
/// A sequential monte carlo move
class Smc_move
{
public:
    explicit Smc_move(sts::likelihood::Forest_likelihood& log_likelihood) : log_likelihood(log_likelihood) {};

    /// Function call for use with smctc - calls user-defined do_move, tracks result.
    int operator()(long, smc::particle<particle::Particle>&, smc::rng*);

    /// Override in subclass with SMC move
    virtual int do_move(long, smc::particle<particle::Particle>&, smc::rng*) const = 0;

    /// Number of times this move has been performed
    int call_count;

    virtual ~Smc_move() {};
protected:
    /// Likelihood calculator
    sts::likelihood::Forest_likelihood log_likelihood;
};

int Smc_move::operator()(long t, smc::particle<particle::Particle>& p, smc::rng* r)
{
    call_count++;
    return do_move(t, p, r);
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_SMC_MOVE_HPP
