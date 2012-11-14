/// \file eb_bl_proposer.h
/// \brief Empirical bayes branch length proposals

#ifndef STS_MOVES_ML_BL_PROPOSER_H
#define STS_MOVES_ML_BL_PROPOSER_H

#include <cassert>
#include <cmath>
#include <limits>
#include <utility>
#include "gsl/gsl_randist.h"
#include "smctc.hh"

#include "edge.h"
#include "particle.h"
#include "forest_likelihood.h"
#include "state.h"
#include "branch_length_proposer.h"
#include "lcfit.h"

namespace sts
{
namespace moves
{
namespace impl
{

class Binary_search_bl
{
public:
    Binary_search_bl(const likelihood::Forest_likelihood& fl, const particle::Particle part) :
        fl(fl),
        part(part),
        calc(fl.get_calculator()) {};
    double operator()(const double d) const;
private:
    const likelihood::Forest_likelihood fl;
    const particle::Particle part;
    const std::shared_ptr<likelihood::Online_calculator> calc;
};

/// Evaluate the log-likelihood with both child branch lengths set to \c d.
/// \returns log-likelihood.
double Binary_search_bl::operator()(const double d) const
{
    calc->invalidate(part->node);
    part->node->child1->length = d;
    part->node->child2->length = d;
    double ll = fl(part);
    return ll;
}

}

/// \class Eb_bl_proposer
/// \brief Propose branch lengths using an empirical Bayes procedure.
/// Wraps a branch length proposer, setting the mean value of the proposer to the value obtained from a deterministic
/// binary search.
template <class T>
class Eb_bl_proposer : public Branch_length_proposer
{
public:
    Eb_bl_proposer(likelihood::Forest_likelihood& fl, T wrapped, int n_iters) :
        fl(fl),
        n_iters(n_iters),
        delta(0.25),
        initial_bl(0.5),
        wrapped(wrapped) {};
    double log_proposal_density(double);
    double operator()(particle::Particle, smc::rng*);

protected:
    likelihood::Forest_likelihood fl;
    int n_iters;
    double delta;
    double initial_bl;
    T wrapped;
    double estimate_proposal_dist_mean(particle::Particle *);
    double estimate_proposal_dist_mean_lcfit(particle::Particle *part);
};


// Implementation
template <class T>
double Eb_bl_proposer<T>::operator()(particle::Particle part, smc::rng* rng)
{
    double orig_mean = wrapped.mean;
    // Estimate the mean of the proposal distribution from the data
    // Initialize child lengths with initial_bl
    part->node->child1->length = initial_bl;
    part->node->child2->length = initial_bl;
    double new_mean = estimate_proposal_dist_mean_lcfit(&part);

    if(1e-16 > new_mean || new_mean > 5. || isnan(new_mean)) {
        // Bad new_mean; re-estimate using bifurcation.
        //std::cout << new_mean << "\t";
        new_mean = estimate_proposal_dist_mean(&part);
        //std::cout << new_mean << "\n";
    }
    //else
    //   std::cout << ".";

    wrapped.mean = new_mean;

    // Proceed as usual
    double result = wrapped(part, rng);
    wrapped.mean = orig_mean;
    return result;
}

/// Generate a mean of the proposal distribution to that obtained by maximum likelihood
/// \param part Input particle
template <class T>
double Eb_bl_proposer<T>::estimate_proposal_dist_mean(particle::Particle *part)
{
    impl::Binary_search_bl f(fl, *part);
    double step = delta;
    double cur_ll = fl(*part);
    double bl = initial_bl;
    for(int i = 0; i < n_iters; ++i) {
        // First try moving right
        double new_ll = f(bl + step);
        if(new_ll > cur_ll) {
            cur_ll = new_ll;
            bl = bl + step;
            continue;
        }
        // Then left, constrained to positive values
        if(bl > step) {
            new_ll = f(bl - step);
            if(new_ll > cur_ll) {
                cur_ll = new_ll;
                bl = bl - step;
                continue;
            }
        }
        // Both attempts failed: halve step size
        step /= 2.;
    }

    //std::cout << cur_ll << "\n";

    fl.get_calculator()->invalidate((*part)->node);

    return bl;
}

/// Generate a mean of the proposal distribution to that obtained by maximum likelihood
/// \param part Input particle
template <class T>
double Eb_bl_proposer<T>::estimate_proposal_dist_mean_lcfit(particle::Particle *part)
{
    impl::binary_search_bl f(fl, *part);
    double x[3] = {1e10,1e2,1e-10};
    const size_t n_points = 5;
    double l[n_points],t[n_points] = {0.1, 0.2, 0.3, 1., 5.};
    for (int i=0; i < n_points; i++) {
        l[i] = f(t[i]);
    }
    int status = fit_ll(3, t, l, x);

    //for(int i=0; i<3; i++) std::cout << x[i] << "\t";
    //std::cout << "\n";

    fl.get_calculator()->invalidate((*part)->node);

    return ml_t(x);
}

template<class T>
double Eb_bl_proposer<T>::log_proposal_density(double d)
{
    return wrapped.log_proposal_density(d);
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_ML_BL_PROPOSER_H
