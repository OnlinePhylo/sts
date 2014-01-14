#include "lcfit_online_add_sequence_move.h"

#include <memory>
#include <string>
#include <vector>
#include <gsl/gsl_integration.h>

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
    bsm_t model_;

    double ml_t_;
    double ml_ll_;

    double t_min_;
    double t_max_;

public:
    LcfitRejectionSampler(smc::rng* rng, const bsm_t& model) :
        rng_(rng), model_(model)
    {
        ml_t_ = lcfit_bsm_ml_t(&model);
        ml_ll_ = lcfit_bsm_log_like(ml_t_, &model);

        const double STEP = 1e-3;
        for (t_min_ = ml_t_; lcfit_bsm_log_like(t_min_, &model) - ml_ll_ > -10.0; t_min_ -= STEP);
        for (t_max_ = ml_t_; lcfit_bsm_log_like(t_max_, &model) - ml_ll_ > -10.0; t_max_ += STEP);

        assert(std::isfinite(t_min_));
        assert(std::isfinite(t_max_));
    }

    const std::pair<double, double> sample() const {
        double t = 0.0;
        double y = 0.0;

        do {
            // TODO: fix ranges of t and y
            t = rng_->Uniform(t_min_, t_max_);
            y = rng_->Uniform(0.0, 1.0);
        } while (y > p(t));

        return std::make_pair(t, std::log(y));
    }

private:

    double p(double t) const {
        return std::exp(lcfit_bsm_log_like(t, &model_) - ml_ll_);
    }
};

LcfitOnlineAddSequenceMove::LcfitOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                                       const std::vector<std::string>& taxaToAdd,
                                                       const std::vector<double>& proposePendantBranchLengths) :
    GuidedOnlineAddSequenceMove(calculator, taxaToAdd, proposePendantBranchLengths)
{ }

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
    // branch lengths
    std::tie(n, edgeLogDensity) = chooseEdge(*tree, leafName, rng);
    double mlDistal, mlPendant;
    TripodOptimizer optim = optimizeBranchLengths(n, leafName, mlDistal, mlPendant);

    const double d = n->getDistanceToFather();
    double distal = -1;

    // Handle very small branch lengths - attach with distal BL of 0
    if(d < 1e-8)
        distal = 0;
    else {
        do {
            distal = rng->NormalTruncated(mlDistal, d / 4, 0.0);
        } while(distal < 0 || distal > d);
    }
    assert(!std::isnan(distal));

    const double distalLogDensity = std::log(gsl_ran_gaussian_pdf(distal - mlDistal, d / 4));
    assert(!std::isnan(distalLogDensity));

    //
    // lcfit magic...
    //

    using namespace std::placeholders;
    auto logLike_f = std::bind(&TripodOptimizer::log_like, &optim, distal, _1, false);

    bsm_t model = DEFAULT_INIT;
    lcfit::LCFitResult result = lcfit::fit_bsm_log_likelihood(logLike_f, model, {0.1, 0.15, 0.5});

    LcfitRejectionSampler s(rng, result.model_fit);

    double pendantBranchLength, pendantLogDensity;
    std::tie(pendantBranchLength, pendantLogDensity) = s.sample();

    //const double pendantBranchLength = rng->Exponential(mlPendant);
    //const double pendantLogDensity = std::log(gsl_ran_exponential_pdf(pendantBranchLength, mlPendant));
    assert(!std::isnan(pendantLogDensity));

    return AttachmentProposal { n, edgeLogDensity, distal, distalLogDensity, pendantBranchLength, pendantLogDensity, mlDistal, mlPendant };
}

}} // namespace sts::online
