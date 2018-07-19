#include "node_sliding_window_mcmc_move.h"
#include "composite_tree_likelihood.h"
#include "multiplier_proposal.h"
#include "online_util.h"
#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace bpp;

namespace sts { namespace online {
    
	SlidingWindowMCMCMove::SlidingWindowMCMCMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator, const std::vector<std::string>& parameters,
                                           const double lambda) : OnlineMCMCMove(calculator, parameters, lambda)
    {}
    
    SlidingWindowMCMCMove::~SlidingWindowMCMCMove()
    {
        // Debug bits
        if(n_attempted > 0) {
            std::clog << "SlidingWindow_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << std::endl;
        }
    }
    
    int SlidingWindowMCMCMove::proposeMove(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
    {
        // Choose an edge at random
        TreeParticle* value = particle.GetValuePointer();
        return proposeMove(*value, rng);
    }
    
    int SlidingWindowMCMCMove::proposeMove(TreeParticle& particle, smc::rng* rng)
    {
		size_t index = 0;
#if defined(_OPENMP)
		index = omp_get_thread_num();
#endif
        //std::cout << value->particleID <<std::endl;
        std::vector<bpp::Node*> nodes = onlineAvailableEdges(*particle.tree);
        size_t idx = rng->UniformDiscrete(0, nodes.size() - 1);
        
        bpp::Node* n = nodes[idx];
        const double orig_dist = n->getDistanceToFather();

		std::unique_ptr<CompositeTreeLikelihood>& calculator = _calculator[index];
        calculator->initialize(*particle.model, *particle.rateDist, *particle.tree);
        
        double orig_ll = calculator->operator()();
        const double max_bl = 100.0;
        const double min_bl = 1e-6;
        double new_dist = orig_dist + (rng->UniformS() - 0.5)*_lambda;
        if ( new_dist > max_bl ) {
            new_dist = 2.0*max_bl-new_dist;
        }
        else if ( new_dist < 1e-6 ) {
            new_dist = 2.0*min_bl -new_dist;
        }
        
        //const Proposal p = positive_real_multiplier(orig_dist, 1e-6, 100.0, lambda, rng);
        n->setDistanceToFather(new_dist);
        double new_ll = calculator->operator()();
        particle.logP = new_ll;
        
        double mh_ratio = std::exp(new_ll - orig_ll);
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
