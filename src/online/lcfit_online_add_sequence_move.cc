#include "lcfit_online_add_sequence_move.h"

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gsl/gsl_cdf.h>
#include <lcfit_cpp.h>

#include "attachment_likelihood.h"
#include "composite_tree_likelihood.h"
#include "guided_online_add_sequence_move.h"
#include "lcfit_rejection_sampler.h"
#include "tree_particle.h"

namespace sts { namespace online {

LcfitOnlineAddSequenceMove::LcfitOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                                       const std::vector<std::string>& taxaToAdd,
                                                       const std::vector<double>& proposePendantBranchLengths,
                                                       const double maxLength,
                                                       const size_t subdivideTop) :
    GuidedOnlineAddSequenceMove(calculator, taxaToAdd, proposePendantBranchLengths, maxLength, subdivideTop),
    lcfit_failures_(0),
    lcfit_attempts_(0)
{ }

LcfitOnlineAddSequenceMove::~LcfitOnlineAddSequenceMove()
{
    const double lcfit_failure_rate = static_cast<double>(lcfit_failures_) / lcfit_attempts_;
    std::clog << "[LcfitOnlineAddSequenceMove] lcfit failure rate = " << lcfit_failure_rate << '\n';
}

const std::tuple<bpp::Node*, double, double, double>
LcfitOnlineAddSequenceMove::chooseMoveLocation(bpp::TreeTemplate<bpp::Node>& tree,
                                               const std::string& leafName,
                                               smc::rng* rng)
{
  bpp::Node* n = nullptr;
  double edgeLogDensity;
  std::tie(n, edgeLogDensity) = chooseEdge(tree, leafName, rng);

   assert(n);
   assert(std::isfinite(edgeLogDensity));

   const double d = n->getDistanceToFather();
   double distalBranchLength;

   // TODO: hack
   const double mlDistal = d / 4.0;

   assert(std::isfinite(mlDistal));

   do {
       distalBranchLength = gsl_ran_rayleigh(rng->GetRaw(), mlDistal);
   } while (distalBranchLength > d);

   const double distalLogDensity = std::log(gsl_ran_rayleigh_pdf(distalBranchLength, mlDistal) /
                                            gsl_cdf_rayleigh_P(d, mlDistal));

   return std::make_tuple(n, edgeLogDensity, distalBranchLength, distalLogDensity);
}

AttachmentProposal LcfitOnlineAddSequenceMove::propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    TreeParticle* treeParticle = particle.GetValuePointer();

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
    double edgeLogDensity = 0.0;
    double distalBranchLength = 0.0;
    double distalLogDensity = 0.0;

    std::tie(n, edgeLogDensity, distalBranchLength, distalLogDensity) = chooseMoveLocation(*(treeParticle->tree), leafName, rng);

    assert(n);
    assert(std::isfinite(distalBranchLength));
    assert(std::isfinite(distalLogDensity));

    AttachmentLikelihood al(calculator, n, leafName, distalBranchLength);

    using namespace std::placeholders;
    auto pendant_ll = std::bind(&AttachmentLikelihood::operator(), &al, _1);

    lcfit::LCFitResult pendantFit = lcfit::fit_bsm_log_likelihood(pendant_ll, DEFAULT_INIT, {0.1, 0.15, 0.5});

    double pendantBranchLength, pendantLogDensity;
    bool lcfitFailure = false;

    try {
        ++lcfit_attempts_;
        lcfit::rejection_sampler sampler(rng, pendantFit);
        pendantBranchLength = sampler.sample();
        pendantLogDensity = sampler.log_density(pendantBranchLength);
    } catch (const std::exception& e) {
        // std::clog << "** " << e.what() << '\n';

        lcfitFailure = true;
        ++lcfit_failures_;

        // TODO: hack
        const double mlPendant = n->getDistanceToFather() / 4.0;
        pendantBranchLength = gsl_ran_rayleigh(rng->GetRaw(), mlPendant);

        // FIXME: if this is an lcfit failure, this log proposal density isn't
        // scaled the same. this could screw up the distribution of the
        // results, as these may look "better" than the lcfit results and so
        // change the weight scaling
        pendantLogDensity = std::log(gsl_ran_rayleigh_pdf(pendantBranchLength, mlPendant));
    }

    assert(std::isfinite(pendantBranchLength));
    assert(std::isfinite(pendantLogDensity));

    return AttachmentProposal { n, edgeLogDensity, distalBranchLength, distalLogDensity, pendantBranchLength, pendantLogDensity, -1.0, -1.0, lcfitFailure, pendantFit };
}

}} // namespace sts::online
