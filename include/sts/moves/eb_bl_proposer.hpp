#ifndef STS_MOVES_ML_BL_PROPOSER_HPP
#define STS_MOVES_ML_BL_PROPOSER_HPP

#include <cassert>
#include <cmath>
#include <limits>
#include <utility>
#include "gsl/gsl_randist.h"
#include "smctc.hh"

#include "sts/likelihood/forest_likelihood.hpp"
#include "sts/particle/phylo_particle.hpp"
#include "sts/moves/branch_length_proposer.hpp"

namespace sts
{
namespace moves
{

/// \class eb_bl_proposer
/// \brief Propose branch lengths using an empirical bayes procedure
/// Wraps a branch length proposer, setting the mean value of the proposer to the value obtained from a fixed-length
/// binary search.
template <class T>
class eb_bl_proposer : public branch_length_proposer
{
public:
    eb_bl_proposer(likelihood::forest_likelihood& fl, T wrapped, int n_iters) : fl(fl), n_iters(n_iters), delta(0.1), wrapped(wrapped), initial_bl(0.11) {};
    double log_proposal_density(double);
    double operator()(particle::particle, smc::rng*);
    branch_lengths propose(particle::particle part, smc::rng *rng) { return wrapped.propose(part, rng); };

protected:
    likelihood::forest_likelihood fl;
    int n_iters;
    double delta;
    double initial_bl;
    T wrapped;
    double estimate_proposal_dist_mean(particle::particle *);
};


// Implementation
template <class T>
double eb_bl_proposer<T>::operator()(particle::particle part, smc::rng* rng)
{
    // Estimate the mean of the proposal distribution from the data
    // Initialize child lengths with initial_bl
    part->node->child1->length = initial_bl;
    part->node->child2->length = initial_bl;
    double new_mean = estimate_proposal_dist_mean(&part);
    wrapped.mean = new_mean;

    // Proceed as usual
    return branch_length_proposer::operator()(part, rng);
}

/// Generate a mean of the proposal distribution to that obtained by maximum likelihood

/// \param part Input particle
template <class T>
double eb_bl_proposer<T>::estimate_proposal_dist_mean(particle::particle *part)
{
    const auto calc = fl.get_calculator();
    std::shared_ptr<particle::phylo_node> node = (*part)->node;
    double cur_ll = fl(*part);
    double new_ll;
    double *v1 = &node->child1->length, *v2 = &node->child2->length;
    const double orig1 = *v1, orig2 = *v2;
    // Current BLs
    double best1 = orig1, best2 = orig2;

    for(int i = 0; i < n_iters; ++i) {
        // First try moving delta to the left.
        *v1 = orig1 - delta;
        *v2 = orig2 - delta;

        // Clear cache
        calc->invalidate(node);
        if((new_ll = fl(*part)) > cur_ll) {
            cur_ll = new_ll;
            best1 = *v1;
            best2 = *v2;
        } else {
            // Next try moving delta to the right.
            *v1 = orig1 + delta;
            *v2 = orig2 + delta;
            calc->invalidate(node);
            if((new_ll = fl(*part)) > cur_ll) {
                cur_ll = new_ll;
                best1 = *v1;
                best2 = *v2;
            } else {
                // Neither try worked; return to original value.
                *v1 = best1;
                *v2 = best2;
            }
        }

        delta /= 2.0;
    }

    // Return mean of two values.
    double result = (*v1 + *v2) / 2.;
    // Restore original state
    *v1 = orig1;
    *v2 = orig2;
    return result;
}

template<class T>
double eb_bl_proposer<T>::log_proposal_density(double d)
{
    return wrapped.log_proposal_density(d);
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_ML_BL_PROPOSER_HPP
