/**
 * \file lcfit_rejection_sampler.h
 * \brief lcfit C++ rejection sampler
 *
 * This file provides a rejection sampler for sampling under the likelihood
 * curve estimated by lcfit.
 */

#ifndef LCFIT_REJECTION_SAMPLER_H
#define LCFIT_REJECTION_SAMPLER_H

#include <tuple>
#include <lcfit_cpp.h>

namespace smc { class rng; }	// forward declaration

namespace lcfit {

class rejection_sampler {
private:
    smc::rng* rng_;
    lcfit::LCFitResult fit_result_;

    double ml_t_;
    double ml_ll_;

    double t_min_;
    double t_max_;

    double auc_;

    mutable size_t n_trials_;
    mutable size_t n_accepts_;

public:
    rejection_sampler(smc::rng* rng, const lcfit::LCFitResult& fitResult);
    virtual ~rejection_sampler();

    const std::pair<double, double> sample() const;

private:
    const std::pair<double, double> find_bounds_easy() const;
    const std::pair<double, double> find_bounds_dumb() const;

    double integrate() const;

    void print_model() const;
    void print_acceptance_rate() const;
};

} // namespace lcfit

#endif // LCFIT_REJECTION_SAMPLER_H
