#include "lcfit_online_add_sequence_move.h"

#include <memory>
#include <string>
#include <vector>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/quadrature/adaptive.hpp>
#include <gsl/gsl_cdf.h>

#include "composite_tree_likelihood.h"
#include "guided_online_add_sequence_move.h"
#include "lcfit.h"
#include "lcfit_cpp.h"
#include "tree_particle.h"
#include "tripod_optimizer.h"

namespace sts { namespace online {

class LcfitRejectionSampler {
private:
    smc::rng* rng_;
    lcfit::LCFitResult fit_result_;

    double ml_t_;
    double ml_ll_;

    double t_min_;
    double t_max_;

    double auc_;

public:
    LcfitRejectionSampler(smc::rng* rng, const lcfit::LCFitResult& fitResult) :
        rng_(rng), fit_result_(fitResult)
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

    const std::pair<double, double> sample() const {
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

private:
    const std::pair<double, double> find_bounds_easy() const {
        // lcfit asserts that the list of selected points is sorted in
        // increasing order by x, so we don't even have to search for the min
        // and the max.
        const auto& points = fit_result_.evaluated_points;
        return std::make_pair(points.front().x, points.back().x);
    }

    const std::pair<double, double> find_bounds(double ll_threshold=-10.0) const {
        const bsm_t* model = &(fit_result_.model_fit);
        auto f = [=](double t) -> double { return lcfit_bsm_log_like(t, model) - ml_ll_ - ll_threshold; };

        std::pair<double, double> bounds;
        boost::math::tools::eps_tolerance<double> tolerance(30);
        boost::uintmax_t max_iters;

        double lower_limit = -1.0;
        while (f(lower_limit) > 0.0) lower_limit *= 2.0;

        double upper_limit = 1.0;
        while (f(upper_limit) > 0.0) upper_limit *= 2.0);

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

    double integrate() const {
        const bsm_t* model = &(fit_result_.model_fit);
        auto f = [=](double t) -> double { return std::exp(lcfit_bsm_log_like(t, model) - ml_ll_); };

        double result = 0.0;
        double error = 0.0;
        boost::numeric::quadrature::adaptive()(f, t_min_, t_max_, result, error);

        assert(error <= 1e-3);
        return result;
    }
};

LcfitOnlineAddSequenceMove::LcfitOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                                       const std::vector<std::string>& taxaToAdd,
                                                       const std::vector<double>& proposePendantBranchLengths) :
    GuidedOnlineAddSequenceMove(calculator, taxaToAdd, proposePendantBranchLengths),
    lcfit_failures_(0),
    lcfit_attempts_(0)
{ }

LcfitOnlineAddSequenceMove::~LcfitOnlineAddSequenceMove()
{
    const double lcfit_failure_rate = static_cast<double>(lcfit_failures_) / lcfit_attempts_;
    std::clog << "[LcfitOnlineAddSequenceMove] lcfit failure rate = " << lcfit_failure_rate << '\n';
}

AttachmentProposal LcfitOnlineAddSequenceMove::propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    TreeParticle* value = particle.GetValuePointer();
    std::unique_ptr<bpp::TreeTemplate<bpp::Node>>& tree = value->tree;

    // Replace node `n` in the tree with a new node containing as children `n` and `new_node`
    // Attach a new leaf, in the following configuration
    //
    //              father
    //   /          o
    //   |          | d - distal
    //   |          |
    // d | new_node o-------o new_leaf
    //   |          |
    //   |          | distal
    //   \          o
    //              n

    bpp::Node* n = nullptr;
    double edgeLogDensity;
    std::tie(n, edgeLogDensity) = chooseEdge(*tree, leafName, rng);

    assert(n);
    assert(std::isfinite(edgeLogDensity));

    double mlDistal, mlPendant;
    TripodOptimizer optim = optimizeBranchLengths(n, leafName, mlDistal, mlPendant);

    assert(std::isfinite(mlDistal));
    assert(std::isfinite(mlPendant));

    const double d = n->getDistanceToFather();
    double distalBranchLength;

    do {
        distalBranchLength = gsl_ran_rayleigh(rng->GetRaw(), mlDistal);
    } while (distalBranchLength > d);

    const double distalLogDensity = std::log(gsl_ran_rayleigh_pdf(distalBranchLength, mlDistal) /
                                             gsl_cdf_rayleigh_P(d, mlDistal));

    assert(std::isfinite(distalBranchLength));
    assert(std::isfinite(distalLogDensity));

    // Update the log-like calculator's cached values for the newly-chosen
    // distal branch length (we can throw away the result, we don't need it for
    // anything).
    optim.log_like(distalBranchLength, mlPendant, true);

    //
    // lcfit magic...
    //

    using namespace std::placeholders;
    auto pendant_ll = std::bind(&TripodOptimizer::log_like, &optim, distalBranchLength, _1, false);

    lcfit::LCFitResult pendantFit = lcfit::fit_bsm_log_likelihood(pendant_ll, DEFAULT_INIT, {0.1, 0.15, 0.5});

    double pendantBranchLength, pendantLogDensity;
    bool lcfitFailure = false;

    try {
        ++lcfit_attempts_;
        std::tie(pendantBranchLength, pendantLogDensity) = LcfitRejectionSampler(rng, pendantFit).sample();
    } catch (const std::exception& e) {
        std::clog << e.what() << '\n';

        lcfitFailure = true;
        ++lcfit_failures_;
        pendantBranchLength = gsl_ran_rayleigh(rng->GetRaw(), mlPendant);
        pendantLogDensity = std::log(gsl_ran_rayleigh_pdf(pendantBranchLength, mlPendant));
    }

    assert(std::isfinite(pendantBranchLength));
    assert(std::isfinite(pendantLogDensity));

    return AttachmentProposal { n, edgeLogDensity, distalBranchLength, distalLogDensity, pendantBranchLength, pendantLogDensity, mlDistal, mlPendant, lcfitFailure, pendantFit.evaluated_points };
}

}} // namespace sts::online
