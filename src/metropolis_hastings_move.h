#ifndef STS_MOVES_METROPOLIS_HASTINGS_MOVE_H
#define STS_MOVES_METROPOLIS_HASTINGS_MOVE_H

#include "forest_likelihood.h"
#include "node.h"
#include "state.h"
#include "smctc.hh"

#include <functional>
#include <string>
#include <vector>

namespace sts
{

// Forward
namespace log { class Json_logger; }

namespace moves
{

struct Mcmc_event;

/// A callback function
typedef std::function<void(const Mcmc_event&)> Mcmc_callback;

/// Abstract class for Metropolis Hastings moves in smctc.
class Metropolis_hastings_move
{
public:
    /// Create a Metropolis_hastings_move
    ///  \param log_likelihood Forest_likelihood to use for likelihood calculations
    explicit Metropolis_hastings_move(likelihood::Forest_likelihood* log_likelihood) : attempted(0), accepted(0),
        log_likelihood(log_likelihood) {};
    virtual ~Metropolis_hastings_move() {};
    /// Number of attempted moves
    unsigned int attempted;
    /// Number of accepted moves
    unsigned int accepted;

    int operator()(long, smc::particle<particle::Particle>&, smc::rng *);
    double acceptance_probability() const;

    /// Override in subclass with MH move

    /// Function proposes a move to \c part, which has a new LL calculated for the new state.
    /// \param time Rank
    /// \param part Particle, of which <b>the current node</b> may be manipulated safely.
    /// \param rng  Random number generator
    /// \returns The size of the move attempted.
    virtual double propose_move(long time, particle::Particle& part, smc::rng* rng) const = 0;

    /// Register a callback function to be called with every proposed move.
    void add_callback(Mcmc_callback callback);

    /// Gets the name of this move
    virtual std::string get_name() const = 0;
protected:
    likelihood::Forest_likelihood* log_likelihood;
    std::vector<Mcmc_callback> callbacks;
};
} // namespace sts::moves
} // namespace sts
#endif // STS_MOVES_METROPOLIS_HASTINGS_MOVE_H
