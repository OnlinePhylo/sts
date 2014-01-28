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
        std::clog << "\n** lcfit failure **\n";
        print_model();
        throw std::runtime_error("lcfit failure: non-finite ML branch length");
    }

    ml_t_ = std::max(0.0, ml_t_);
    ml_ll_ = lcfit_bsm_log_like(ml_t_, model);

    if (!std::isfinite(ml_ll_)) {
        std::clog << "\n** lcfit failure **\n";
        print_model();
        throw std::runtime_error("lcfit failure: non-finite ML log-likelihood");
    }

    std::tie(t_min_, t_max_) = find_bounds_dumb();

    if (t_max_ <= 0.0) {
        std::clog << "\n** lcfit failure **\n";
        print_model();
        throw std::runtime_error("lcfit failure: invalid proposal range");
    }

    auc_ = integrate();
}

rejection_sampler::~rejection_sampler()
{ }

const std::pair<double, double> rejection_sampler::sample() const
{
    const bsm_t* model = &(fit_result_.model_fit);

    auto f = [=](double t) -> double { return std::exp(lcfit_bsm_log_like(t, model) - ml_ll_); };

    double t = 0.0;
    double y = 0.0;

    do {
        t = rng_->Uniform(t_min_, t_max_);
        y = rng_->Uniform(0.0, 1.0);
        ++n_trials_;

        if (n_accepts_ == 0 && n_trials_ >= 1000) {
            std::clog << "\n** lcfit failure **\n";
            print_model();
            throw std::runtime_error("lcfit failure: inefficient sampling");
        }

    } while (y > f(t));

    ++n_accepts_;

    return std::make_pair(t, std::log(f(t) / auc_));
}

const std::pair<double, double> rejection_sampler::find_bounds_easy() const
{
    // lcfit asserts that the list of selected points is sorted in
    // increasing order by x, so we don't even have to search for the min
    // and the max.
    const auto& points = fit_result_.evaluated_points;
    return std::make_pair(points.front().x, points.back().x);
}

const std::pair<double, double> rejection_sampler::find_bounds_dumb() const
{
    const bsm_t* model = &(fit_result_.model_fit);

    // preconditions: ml_t_ >= 0.0, ml_ll_ exists at ml_t_
    const double ll_threshold = ml_ll_ - 100.0;
    const double MIN_STEP = 1e-3;

    double t_min;
    if (ml_t_ == 0.0) {
        t_min = 0.0;
    } else {
        double step = MIN_STEP;
        do {
            if (!std::isfinite(step)) {
                std::clog << "\n** lcfit failure **\n";
                print_model();
                throw std::runtime_error("lcfit failure: t_min infinite step");
            }

            t_min = ml_t_ - step;
            step *= 2.0;
        } while (t_min >= 0.0 && lcfit_bsm_log_like(t_min, model) > ll_threshold);

        if (t_min < 0.0) {
            t_min = 0.0;
        }
    }

    double t_max;
    {
        double step = MIN_STEP;
        do {
            if (!std::isfinite(step)) {
                std::clog << "\n** lcfit failure **\n";
                print_model();
                throw std::runtime_error("lcfit failure: t_max infinite step");
            }

            t_max = ml_t_ + step;
            step *= 2.0;
        } while (lcfit_bsm_log_like(t_max, model) > ll_threshold);
    }

    return std::make_pair(t_min, t_max);
}

double rejection_sampler::integrate() const
{
    const bsm_t* model = &(fit_result_.model_fit);

    auto f = [=](double t) -> double { return std::exp(lcfit_bsm_log_like(t, model) - ml_ll_); };

    double result = 0.0;
    double error = 0.0;
    boost::numeric::quadrature::adaptive()(f, t_min_, t_max_, result, error);

    assert(error <= 1e-3);
    return result;
}

void rejection_sampler::print_model() const
{
    const bsm_t* model = &(fit_result_.model_fit);

    std::clog << "c = " << model->c
              << ", m = " << model->m
              << ", r = " << model->r
              << ", b = " << model->b
              << '\n';
    std::clog << "ml_t = " << ml_t_
              << ", ml_ll = " << ml_ll_
              << '\n';
    std::clog << "t_min = " << t_min_
              << ", t_max = " << t_max_
              << '\n';
}

void rejection_sampler::print_acceptance_rate() const
{
    std::clog << "[rejection_sampler] acceptance rate = "
              << static_cast<double>(n_accepts_) / n_trials_
              << " [" << n_accepts_ << "/" << n_trials_ << "]\n";
}

} // namespace lcfit
