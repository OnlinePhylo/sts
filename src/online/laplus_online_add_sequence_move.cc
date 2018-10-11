//
//  laplus_online_add_sequence_move.cpp
//  sts
//
//  Created by Mathieu Fourment on 11/10/18.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "laplus_online_add_sequence_move.h"

#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <algorithm>

#include <gsl/gsl_cdf.h>
#include <lcfit_rejection_sampler.h>
#include <lcfit_select.h>

#include "composite_tree_likelihood.h"
#include "guided_online_add_sequence_move.h"
#include "lcfit_rejection_sampler.h"
#include "tree_particle.h"
#include "online_util.h"
#include "weighted_selector.h"
#include "lesplace.h"
#include "gsl.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace sts { namespace online {
	
	LaplusOnlineAddSequenceMove::LaplusOnlineAddSequenceMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator,
														   const std::vector<std::string>& sequenceNames,
														   const std::vector<std::string>& taxaToAdd,
														   const std::vector<double>& proposePendantBranchLengths,
														   const double maxLength,
														   const size_t subdivideTop,
														   const double expPriorMean) :
	GuidedOnlineAddSequenceMove(calculator, sequenceNames, taxaToAdd, proposePendantBranchLengths, maxLength, subdivideTop),
	lcfit_failures_(0),
	lcfit_attempts_(0)
	{
		_proposalMethodName = "LaplusOnlineAddSequenceMove";
	}
	
	LaplusOnlineAddSequenceMove::~LaplusOnlineAddSequenceMove()
	{
		const double lcfit_failure_rate = static_cast<double>(lcfit_failures_) / lcfit_attempts_;
		std::clog << "[LaplusOnlineAddSequenceMove] laplus failure rate = "
		<< lcfit_failures_ << "/" << lcfit_attempts_
		<< " (" << lcfit_failure_rate * 100.0 << "%)\n";
	}
	
	struct WrapperFlexibleTreeLikelihood{
		CompositeTreeLikelihood& ctl;
		const bpp::Node& node;
		const std::string& leafName;
		const double distalLength;
		const double proximalLength;
		const double pendantLength;
	};
	
	double pendant_lnl_callback(void* data, size_t i, double t)
	{
		WrapperFlexibleTreeLikelihood* al = static_cast<WrapperFlexibleTreeLikelihood*>(data);
		return al->ctl(al->node, al->leafName, t, al->distalLength, al->proximalLength);
	}
	
	void pendant_dlnl_callback(void* data, size_t i, double* d1, double* d2)
	{
		WrapperFlexibleTreeLikelihood* al = static_cast<WrapperFlexibleTreeLikelihood*>(data);
		al->ctl.calculatePendantDerivatives(al->node, al->leafName, al->pendantLength, al->distalLength, al->proximalLength, d1, d2);
	}
	
	double distal_lnl_callback(void* data, size_t i, double t)
	{
		WrapperFlexibleTreeLikelihood* al = static_cast<WrapperFlexibleTreeLikelihood*>(data);
		// pendant, distal, proximal
		return al->ctl(al->node, al->leafName, al->pendantLength, t, al->proximalLength-t);
	}
	
	void distal_dlnl_callback(void* data, size_t i, double* d1, double* d2)
	{
		WrapperFlexibleTreeLikelihood* al = static_cast<WrapperFlexibleTreeLikelihood*>(data);
		al->ctl.calculateDistalDerivatives(al->node, al->leafName, al->pendantLength, al->distalLength, al->proximalLength-al->distalLength, d1, d2);
	}
	
	std::pair<double, double> LaplusOnlineAddSequenceMove::proposeDistal(bpp::Node& n, const std::string& leafName, const double mlDistal, const double mlPendant, smc::rng* rng) const
	{
		size_t index = 0;
#if defined(_OPENMP)
		index = omp_get_thread_num();
#endif
		
		assert(mlDistal <= n.getDistanceToFather());
		const double edgeLength = n.getDistanceToFather();
		//CompositeTreeLikelihood* ctl = calculator[index].get();
		//auto fn = [ctl,edgeLength, mlPendant,leafName](double distal) {
		//	return -ctl(n, leafName, mlPendant, distal, edgeLength-distal);
		//};
		//sts::gsl::minimize(fn, mlDistal, 0, 10, 1000);
		WrapperFlexibleTreeLikelihood wftl{*(calculator[index]), n, leafName, mlDistal, edgeLength, mlPendant};
		
		double params[2];
		lesplace_gamma_fit(&distal_lnl_callback, &distal_dlnl_callback, &wftl, &mlDistal, 1, params);
		double distal = rng->Gamma(params[0], params[1]);
		double distalLogDensity = std::log(gsl_ran_gamma_pdf(distal, params[0], params[1]));
		if(distal > edgeLength){
			n.setDistanceToFather(distal);
		}
		return std::pair<double, double>(distal, distalLogDensity);
	}
	
	std::pair<double, double> LaplusOnlineAddSequenceMove::proposePendant(bpp::Node& n, const std::string& leafName, const double mlPendant, const double distalBranchLength, smc::rng* rng) const{
		size_t index = 0;
#if defined(_OPENMP)
		index = omp_get_thread_num();
#endif
		
		WrapperFlexibleTreeLikelihood wftl{*(calculator[index]), n, leafName, distalBranchLength, n.getDistanceToFather()-distalBranchLength, mlPendant};
		
		double params[2];
		lesplace_gamma_fit(&pendant_lnl_callback, &pendant_dlnl_callback, &wftl, &mlPendant, 1, params);
		const double pendantBranchLength = rng->Gamma(params[0], params[1]);
		const double pendantLogDensity = std::log(gsl_ran_gamma_pdf(pendantBranchLength, params[0], params[1]));
		
		assert(std::isfinite(pendantBranchLength));
		assert(std::isfinite(pendantLogDensity));
		return std::pair<double, double>(pendantBranchLength, pendantLogDensity);
	}
	
	
}} // namespace sts::online
