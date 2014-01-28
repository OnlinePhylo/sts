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
    rng_(rng), fit_result_(fit_result)
{
    const bsm_t* model = &(fit_result_.model_fit);
    ml_t_ = lcfit_bsm_ml_t(model);
    ml_ll_ = lcfit_bsm_log_like(ml_t_, model);

    if (!std::isfinite(ml_t_) || !std::isfinite(ml_ll_)) {
        std::clog << "** lcfit failure **\n";
        std::clog << "c = " << model->c << ", m = " << model->m << ", r = " << model->r << ", b = " << model->b << '\n';
        std::clog << "ml_t = " << ml_t_ << ", ml_ll = " << ml_ll_ << '\n';
        throw std::runtime_error("lcfit failure");
    }

    std::tie(t_min_, t_max_) = find_bounds();
    auc_ = integrate();
}

const std::pair<double, double> rejection_sampler::sample() const {
    double t = 0.0;
    double y = 0.0;

    const bsm_t* model = &(fit_result_.model_fit);
    auto f = [=](double t) -> double { return std::exp(lcfit_bsm_log_like(t, model) - ml_ll_); };

    do {
        t = rng_->Uniform(t_min_, t_max_);
        y = rng_->Uniform(0.0, 1.0);
    } while (t < 0.0 || y > f(t));

    return std::make_pair(t, std::log(f(t) / auc_));
}

const std::pair<double, double> rejection_sampler::find_bounds_easy() const {
    // lcfit asserts that the list of selected points is sorted in
    // increasing order by x, so we don't even have to search for the min
    // and the max.
    const auto& points = fit_result_.evaluated_points;
    return std::make_pair(points.front().x, points.back().x);
}

const std::pair<double, double> rejection_sampler::find_bounds(double ll_threshold) const {
    const bsm_t* model = &(fit_result_.model_fit);
    auto f = [=](double t) -> double { return lcfit_bsm_log_like(t, model) - ml_ll_ - ll_threshold; };

    std::pair<double, double> bounds;
    boost::math::tools::eps_tolerance<double> tolerance(30);
    boost::uintmax_t max_iters;

    double lower_limit = -1.0;
    while (f(lower_limit) > 0.0) lower_limit *= 2.0;

    double upper_limit = 1.0;
    while (f(upper_limit) > 0.0) upper_limit *= 2.0;

    assert(lower_limit < ml_t_ && ml_t_ < upper_limit);

    max_iters = 100;
    bounds = boost::math::tools::toms748_solve(f, lower_limit, ml_t_, tolerance, max_iters);
    const double t_min = (bounds.first + bounds.second) / 2.0;
    assert(max_iters <= 100);
    assert(std::isfinite(t_min));

    max_iters = 100;
    bounds = boost::math::tools::toms748_solve(f, ml_t_, upper_limit, tolerance, max_iters);
    const double t_max = (bounds.first + bounds.second) / 2.0;
    assert(max_iters <= 100);
    assert(std::isfinite(t_max));

    assert(t_min < t_max);
    return std::make_pair(t_min, t_max);
}

double rejection_sampler::integrate() const {
    const bsm_t* model = &(fit_result_.model_fit);
    auto f = [=](double t) -> double { return std::exp(lcfit_bsm_log_like(t, model) - ml_ll_); };

    double result = 0.0;
    double error = 0.0;
    boost::numeric::quadrature::adaptive()(f, t_min_, t_max_, result, error);

    assert(error <= 1e-3);
    return result;
}

} // namespace lcfit
