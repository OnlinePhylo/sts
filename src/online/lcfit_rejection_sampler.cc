#include "lcfit_rejection_sampler.h"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/quadrature/adaptive.hpp>

#include <smctc.hh>
#include <lcfit_cpp.h>

namespace lcfit {

rejection_sampler::rejection_sampler(smc::rng* rng, const lcfit::LCFitResult& fit_result) :
    rng_(rng), fit_result_(fit_result), n_trials_(0), n_accepts_(0)
{
    const bsm_t* model = &(fit_result_.model_fit);
    ml_t_ = lcfit_bsm_ml_t(model);

    if (!std::isfinite(ml_t_)) {
        throw std::runtime_error("lcfit failure: non-finite ML branch length");
    }

    ml_ll_ = lcfit_bsm_log_like(ml_t_, model);

    if (!std::isfinite(ml_ll_)) {
        throw std::runtime_error("lcfit failure: non-finite ML log-likelihood");
    }

    std::tie(t_min_, t_max_) = find_bounds();

    if (t_max_ <= 0.0) {
        throw std::runtime_error("lcfit failure: invalid proposal range");
    }

    log_auc_ = std::log(integrate());

    if (!std::isnormal(log_auc_)) {
        throw std::runtime_error("lcfit failure: integration result");
    }
}

rejection_sampler::~rejection_sampler()
{ }

double rejection_sampler::sample() const
{
    double t = 0.0;
    double y = 0.0;

    do {
        t = rng_->Uniform(t_min_, t_max_);
        y = rng_->Uniform(0.0, 1.0);
        ++n_trials_;

        if (n_accepts_ == 0 && n_trials_ >= 1000) {
            throw std::runtime_error("lcfit failure: inefficient sampling");
        }

    } while (y > relative_likelihood(t));

    ++n_accepts_;
    return t;
}

double rejection_sampler::relative_log_likelihood(double t) const
{
    const bsm_t* model = &(fit_result_.model_fit);
    return lcfit_bsm_log_like(t, model) - ml_ll_;
}

double rejection_sampler::relative_likelihood(double t) const
{
    return std::exp(relative_log_likelihood(t));
}

double rejection_sampler::log_density(double t) const
{
    return relative_log_likelihood(t) - log_auc_;
}

double rejection_sampler::density(double t) const
{
    return std::exp(log_density(t));
}

const std::pair<double, double> rejection_sampler::find_bounds_easy() const
{
    // lcfit asserts that the list of selected points is sorted in
    // increasing order by x, so we don't even have to search for the min
    // and the max.
    const auto& points = fit_result_.evaluated_points;
    return std::make_pair(points.front().x, points.back().x);
}

const std::pair<double, double> rejection_sampler::find_bounds() const
{
    // preconditions: ml_t_ >= 0.0, ml_ll_ exists at ml_t_
    const double ll_threshold = -100.0;
    const double MIN_STEP = 1e-3;

    double t_min;
    if (ml_t_ == 0.0) {
        t_min = 0.0;
    } else {
        double step = MIN_STEP;
        do {
            if (!std::isfinite(step)) {
                throw std::runtime_error("lcfit failure: t_min infinite step");
            }

            t_min = ml_t_ - step;
            step *= 2.0;
        } while (t_min >= 0.0 && relative_log_likelihood(t_min) > ll_threshold);

        if (t_min < 0.0) {
            t_min = 0.0;
        }
    }

    double t_max;
    {
        double step = MIN_STEP;
        do {
            if (!std::isfinite(step)) {
                throw std::runtime_error("lcfit failure: t_max infinite step");
            }

            t_max = ml_t_ + step;
            step *= 2.0;
        } while (relative_log_likelihood(t_max) > ll_threshold);
    }

    return std::make_pair(t_min, t_max);
}

double rejection_sampler::integrate() const
{
    auto f = [=](double t) { return relative_likelihood(t); };

    double result = 0.0;
    double error = 0.0;
    boost::numeric::quadrature::adaptive()(f, t_min_, t_max_, result, error);

    if (error > 1e-3) {
        throw std::runtime_error("lcfit failure: integration error");
    }

    return result;
}

void rejection_sampler::print_acceptance_rate() const
{
    std::clog << "[rejection_sampler] acceptance rate = "
              << static_cast<double>(n_accepts_) / n_trials_
              << " [" << n_accepts_ << "/" << n_trials_ << "]\n";
}

} // namespace lcfit
