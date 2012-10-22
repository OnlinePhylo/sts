#ifndef STS_MOVES_ML_BL_PROPOSAL_HPP
#define STS_MOVES_ML_BL_PROPOSAL_HPP

#include <cassert>
#include <cmath>
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

/// \class ml_bl_proposal
/// \brief Perform some limited branch length optimizations
template <class T>
class ml_bl_proposer : public branch_length_proposer
{
public:
    ml_bl_proposer(likelihood::forest_likelihood& fl, T wrapped, int n_iters) : fl(fl), n_iters(n_iters), wrapped(wrapped), delta(0.1) {};
    double log_proposal_density(double);
    double operator()(particle::particle part, smc::rng* rng);
    branch_lengths propose(particle::particle part, smc::rng *rng) { return wrapped.propose(part, rng); };

protected:
    likelihood::forest_likelihood fl;
    int n_iters;
    double delta;
    T wrapped;
    double limited_optimization(particle::particle *part, double *to_opt, double delta);
};


// Implementation
template <class T>
double ml_bl_proposer<T>::operator()(particle::particle part, smc::rng* rng)
{
    branch_length_proposer::operator()(part, rng);

    // Try edge 1
    limited_optimization(&part, &part->node->child1->length, 0.1);

    // And edge 2
    limited_optimization(&part, &part->node->child2->length, 0.1);

    return fl(part);
}

/// Perform a limited number of attempts to optimize \c to_opt
/// TODO: include proposal density

/// \param part Input particle
/// \param to_opt Parameter to optimize
/// \param delta Amount to perturb \c to_opt
/// \param n_moves Maximum number of optimization attempts
template <class T>
double ml_bl_proposer<T>::limited_optimization(particle::particle *part, double *to_opt, double delta)
{
    auto calc = fl.get_calculator();
    int node_id = (*part)->node->id;
    double cur_ll = fl(*part);
    double new_ll;
    double orig = *to_opt;

    assert(n_iters > 0);

    for(int i = 0; i < n_iters; ++i) {
        *to_opt = orig - delta;
        if((new_ll = fl(*part)) > cur_ll) {
            cur_ll = new_ll;
        } else {
            calc->invalidate(node_id);
            *to_opt = orig + delta;
            if((new_ll = fl(*part)) > cur_ll) {
                cur_ll = new_ll;
            } else {
                *to_opt = orig;
            }
        }

        delta /= 2.0;
    }

    return cur_ll;
}

template<class T>
double ml_bl_proposer<T>::log_proposal_density(double d)
{
    return wrapped.log_proposal_density(d);
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_ML_BL_PROPOSAL_HPP
