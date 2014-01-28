#include "lcfit_online_add_sequence_move.h"

#include <memory>
#include <string>
#include <vector>

#include <gsl/gsl_cdf.h>
#include <lcfit_cpp.h>

#include "composite_tree_likelihood.h"
#include "guided_online_add_sequence_move.h"
#include "lcfit_rejection_sampler.h"
#include "tree_particle.h"
#include "tripod_optimizer.h"

namespace sts { namespace online {

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
        std::tie(pendantBranchLength, pendantLogDensity) = lcfit::rejection_sampler(rng, pendantFit).sample();
    } catch (const std::exception& e) {
        std::clog << e.what() << '\n';

        lcfitFailure = true;
        ++lcfit_failures_;
        pendantBranchLength = gsl_ran_rayleigh(rng->GetRaw(), mlPendant);
        pendantLogDensity = std::log(gsl_ran_rayleigh_pdf(pendantBranchLength, mlPendant));
    }

    assert(std::isfinite(pendantBranchLength));
    assert(std::isfinite(pendantLogDensity));

    return AttachmentProposal { n, edgeLogDensity, distalBranchLength, distalLogDensity, pendantBranchLength, pendantLogDensity, mlDistal, mlPendant, lcfitFailure, pendantFit };
}

}} // namespace sts::online
