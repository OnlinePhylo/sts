#include "smc_move.h"

namespace sts
{
namespace moves
{
int Smc_move::operator()(long t, smc::particle<particle::Particle>& p, smc::rng* r)
{
    call_count++;
    return do_move(t, p, r);
}

} // namespace moves
} // namespace sts

