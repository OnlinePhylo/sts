#include "lcfit_online_add_sequence_move.h"

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gsl/gsl_cdf.h>
#include <lcfit_rejection_sampler.h>
#include <lcfit_select.h>

#include "composite_tree_likelihood.h"
#include "guided_online_add_sequence_move.h"
#include "lcfit_rejection_sampler.h"
#include "tree_particle.h"
#include "online_util.h"
#include "weighted_selector.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace sts { namespace online {

LcfitOnlineAddSequenceMove::LcfitOnlineAddSequenceMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator,
                                                       const std::vector<std::string>& sequenceNames,
                                                       const std::vector<std::string>& taxaToAdd,
                                                       const std::vector<double>& proposePendantBranchLengths,
                                                       const double maxLength,
                                                       const size_t subdivideTop,
                                                       const double expPriorMean) :
    GuidedOnlineAddSequenceMove(calculator, sequenceNames, taxaToAdd, proposePendantBranchLengths, maxLength, subdivideTop),
    expPriorMean_(expPriorMean),
    lcfit_failures_(0),
    lcfit_attempts_(0)
{
    _proposalMethodName = "LcfitOnlineAddSequenceMove";
}

LcfitOnlineAddSequenceMove::~LcfitOnlineAddSequenceMove()
{
    const double lcfit_failure_rate = static_cast<double>(lcfit_failures_) / lcfit_attempts_;
    std::clog << "[LcfitOnlineAddSequenceMove] lcfit failure rate = "
              << lcfit_failures_ << "/" << lcfit_attempts_
              << " (" << lcfit_failure_rate * 100.0 << "%)\n";
}

    struct WrapperFlexibleTreeLikelihood{
        CompositeTreeLikelihood& ctl;
        const bpp::Node& node;
        const std::string& leafName;
        const double distalLength;
        const double proximalLength;
    };
    
double attachment_lnl_callback(double t, void* data)
{
    WrapperFlexibleTreeLikelihood* al = static_cast<WrapperFlexibleTreeLikelihood*>(data);

    return al->ctl(al->node, al->leafName, t, al->distalLength, al->proximalLength);
}

std::pair<double, double> LcfitOnlineAddSequenceMove::proposeDistal(bpp::Node& n, const std::string& leafName, const double mlDistal, const double mlPendant, smc::rng* rng) const
{
	size_t index = 0;
#if defined(_OPENMP)
	index = omp_get_thread_num();
#endif
	
    assert(mlDistal <= n.getDistanceToFather());
    const double edgeLength = n.getDistanceToFather();
    
    double distal = -1;

    double dd1, dd2;
	//calculator(n, leafName, mlPendant, mlDistal, edgeLength-mlDistal);
    calculator[index]->calculateDistalDerivatives(n, leafName, mlPendant, mlDistal, edgeLength-mlDistal, &dd1, &dd2);
    const double sigma = sqrt(fabs(1/dd2));
    
    // Handle very small branch lengths - attach with distal BL of 0
    if(edgeLength < 1e-8)
        distal = 0;
    else {
        do {
            distal = rng->NormalTruncated(mlDistal, sigma, 0.0);
        } while(distal < 0 || distal > edgeLength);
    }
//    if(std::isnan(distal)) distal= mlDistal;
           //std::cout << dd1 << " "<<dd2<<" "<<distal<<" "<<edgeLength << " "<<mlDistal<<" "<<mlPendant<<std::endl;
    assert(!std::isnan(distal));
    
    // Log density: for CDF F(x) and PDF g(x), limited to the interval (a, b]:
    //
    // g'(x) =   g(x) / [F(b) - F(a)]
    //
    // We are limited to (0, d].
    //
    // GSL gaussian CDFs are for mean 0, hence the mlDistal substraction here.
    const double distalLogDensity = std::log(gsl_ran_gaussian_pdf(distal - mlDistal, sigma)) -
    std::log(gsl_cdf_gaussian_P(edgeLength - mlDistal, sigma) - gsl_cdf_gaussian_P(- mlDistal, sigma));
    assert(!std::isnan(distalLogDensity));
//     std::cout << mlPendant << " "<<mlDistal<<" "<<dd1<<" "<<dd2<<" "<<distal<<" "<<distalLogDensity<<std::endl;
    return std::pair<double, double>(distal, distalLogDensity);
}
    
std::pair<double, double> LcfitOnlineAddSequenceMove::proposePendant(bpp::Node& n, const std::string& leafName, const double mlPendant, const double distalBranchLength, smc::rng* rng) const{
	size_t index = 0;
#if defined(_OPENMP)
	index = omp_get_thread_num();
#endif
	
    // FIXME: Are there actual branch length constraints available somewhere?
    const double min_t = 1e-6;
    const double max_t = 20.0;
    bsm_t model = DEFAULT_INIT;
    
//    _al->initialize(&n, leafName, distalBranchLength);
    WrapperFlexibleTreeLikelihood wftl{*(calculator[index]), n, leafName, distalBranchLength, n.getDistanceToFather()-distalBranchLength};
    lcfit_fit_auto(&attachment_lnl_callback, &wftl, &model, min_t, max_t);

//    _al->finalize();

    const double mlPendantBranchLength = lcfit_bsm_ml_t(&model);
    double pendantBranchLength, pendantLogDensity;

    try {
        ++lcfit_attempts_;
        lcfit::rejection_sampler sampler(rng->GetRaw(), model, 1.0 / expPriorMean_);
        pendantBranchLength = sampler.sample();
        pendantLogDensity = sampler.log_density(pendantBranchLength);
    } catch (const std::exception& e) {
        // std::clog << "** " << e.what() << '\n';
        ++lcfit_failures_;

        // Fall back on original proposal
        return GuidedOnlineAddSequenceMove::proposePendant(n, leafName, mlPendant, distalBranchLength, rng);
    }

    assert(std::isfinite(pendantBranchLength));
    assert(std::isfinite(pendantLogDensity));
    return std::pair<double, double>(pendantBranchLength, pendantLogDensity);
}


}} // namespace sts::online
