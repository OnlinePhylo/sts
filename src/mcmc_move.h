#ifndef STS_MOVES_MCMC_MOVE_HPP
#define STS_MOVES_MCMC_MOVE_HPP

#include "forest_likelihood.h"
#include "state.h"

namespace sts
{
namespace moves
{

/// Abstract class for MCMC moves in smctc.
class Mcmc_move
{
public:
    /// Create an Mcmc_move
    ///  \param log_likelihood Forest_likelihood to use for likelihood calculations
    explicit Mcmc_move(likelihood::Forest_likelihood& log_likelihood) : attempted(0), accepted(0), log_likelihood(log_likelihood) {};
    virtual ~Mcmc_move() {};
    /// Number of attempted moves
    unsigned int attempted;
    /// Number of accepted moves
    unsigned int accepted;

    int operator()(long, smc::particle<particle::Particle>&, smc::rng *);

    /// Override in subclass with MCMC move
    virtual int do_move(long, smc::particle<particle::Particle>&, smc::rng*) const = 0;

protected:
    likelihood::Forest_likelihood log_likelihood;
};
} // namespace sts::moves
} // namespace sts
#endif // STS_MOVES_MCMC_MOVE_HPP
