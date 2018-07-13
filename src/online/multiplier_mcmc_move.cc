#include "multiplier_mcmc_move.h"
#include "composite_tree_likelihood.h"
#include "multiplier_proposal.h"
#include "online_util.h"
#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>

using namespace bpp;

namespace sts { namespace online {

MultiplierMCMCMove::MultiplierMCMCMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator, const std::vector<std::string>& parameters,
                                           const double lambda) :
    OnlineMCMCMove(calculator, parameters, lambda)
{}

MultiplierMCMCMove::~MultiplierMCMCMove()
{
    // Debug bits
    if(n_attempted > 0) {
		if(_parameters.size() == 1){
			std::clog << "Multiplier_mcmc_move:"<< _parameters[0] << ": " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << " " << _lambda << std::endl;
		}
		else{
			std::clog << "Multiplier_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << std::endl;
		}
    }
}

int MultiplierMCMCMove::proposeMove(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    // Choose an edge at random
    TreeParticle* value = particle.GetValuePointer();
    return proposeMove(*value, rng);
}
    
int MultiplierMCMCMove::proposeMove(TreeParticle& particle, smc::rng* rng){
	size_t index = 0;
#if defined(_OPENMP)
	index = omp_get_thread_num();
#endif
	
	if(_parameters.size() == 1){
		double orig_value;
		double lower;
		double upper;
	
		
		if (particle.model->hasParameter(_parameters[0])) {
			orig_value = particle.model->getParameterValue(_parameters[0]);
			const IntervalConstraint* constraint = dynamic_cast<const IntervalConstraint*>(particle.model->getParameter(_parameters[0]).getConstraint());
			lower = constraint->getLowerBound();
			upper = constraint->getUpperBound();
		}
		else if (particle.rateDist->hasParameter(_parameters[0])) {
			orig_value = particle.rateDist->getParameterValue(_parameters[0]);
			const IntervalConstraint* constraint = dynamic_cast<const IntervalConstraint*>(particle.rateDist->getParameter(_parameters[0]).getConstraint());
			lower = constraint->getLowerBound();
			upper = constraint->getUpperBound();
		}
		else{
			std::cerr << _parameters[0] << " " << " not found"<< std::endl;
			exit(10);
		}

		_calculator[index]->initialize(*particle.model, *particle.rateDist, *particle.tree);
		
		double orig_ll = particle.logP;
		
		const Proposal p = positive_real_multiplier(orig_value, lower, upper, _lambda, rng);
		if (particle.model->hasParameter(_parameters[0])) {
			particle.model->setParameterValue(_parameters[0], p.value);
		}
		else if (particle.rateDist->hasParameter(_parameters[0])) {
			particle.rateDist->setParameterValue(_parameters[0], p.value);
		}
		
		_calculator[index]->initialize(*particle.model, *particle.rateDist, *particle.tree);
		
		double new_ll = _calculator[index]->operator()();
		
		particle.logP = new_ll;
		double mh_ratio = std::exp(new_ll + std::log(p.hastingsRatio) - orig_ll);
		if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
			return 1;
		} else {
			// Rejected
			particle.logP = orig_ll;
			if (particle.model->hasParameter(_parameters[0])) {
				particle.model->setParameterValue(_parameters[0], orig_value);
			}
			else if (particle.rateDist->hasParameter(_parameters[0])) {
				particle.rateDist->setParameterValue(_parameters[0], orig_value);
			}
			return 0;
		}
	}
	
    std::vector<bpp::Node*> nodes = onlineAvailableEdges(*particle.tree);
    size_t idx = rng->UniformDiscrete(0, nodes.size() - 1);
    
    bpp::Node* n = nodes[idx];
    const double orig_dist = n->getDistanceToFather();
    
    _calculator[index]->initialize(*particle.model, *particle.rateDist, *particle.tree);
    
    double orig_ll = particle.logP;
    
    const Proposal p = positive_real_multiplier(orig_dist, 1e-6, 100.0, _lambda, rng);
    n->setDistanceToFather(p.value);
    double new_ll = _calculator[index]->operator()();
    
    particle.logP = new_ll;
    double mh_ratio = std::exp(new_ll + std::log(p.hastingsRatio) - orig_ll);
    if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
        return 1;
    } else {
        // Rejected
        particle.logP = orig_ll;
        n->setDistanceToFather(orig_dist);
        return 0;
    }
}
    
}}
