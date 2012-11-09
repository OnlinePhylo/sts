#ifndef STS_MOVES_MCMC_MOVE_HPP
#define STS_MOVES_MCMC_MOVE_HPP

#include "sts/likelihood/forest_likelihood.hpp"
#include "sts/particle/state.hpp"

namespace sts
{
namespace moves
{

/// Abstract class for MCMC moves in smctc.
class mcmc_move
{
public:
    /// Create an mcmc_move
    ///  \param log_likelihood forest_likelihood to use for likelihood calculations
    explicit mcmc_move(likelihood::forest_likelihood& log_likelihood) : attempted(0), accepted(0), log_likelihood(log_likelihood) {};
    virtual ~mcmc_move() {};
    /// Number of attempted moves
    unsigned int attempted;
    /// Number of accepted moves
    unsigned int accepted;

    int operator()(long time, smc::particle<particle::particle>& from, smc::rng *rng);

    /// Override in subclass with MCMC move
    virtual int do_move(long, smc::particle<particle::particle>&, smc::rng*) const = 0;

protected:
    likelihood::forest_likelihood log_likelihood;
};

/// Function call for use with smctc - calls user-defined do_move, tracks result.

/// \param time Generation number
/// \param from Source particle
/// \param rng Random number source
int mcmc_move::operator()(long time, smc::particle<particle::particle>& from, smc::rng *rng)
{
    attempted++;
    int result = do_move(time, from, rng);
    if(result) accepted++;
    return result;
}

} // namespace sts::moves
} // namespace sts
#endif // STS_MOVES_MCMC_MOVE_HPP
