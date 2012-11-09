#include "mcmc_move.h"

namespace sts
{
namespace moves
{
/// Function call for use with smctc - calls user-defined do_move, tracks result.

/// \param time Generation number
/// \param from Source particle
/// \param rng Random number source
int Mcmc_move::operator()(long time, smc::particle<particle::Particle>& from, smc::rng *rng)
{
    attempted++;
    int result = do_move(time, from, rng);
    if(result) accepted++;
    return result;
}

} // namespace sts::moves
} // namespace sts
