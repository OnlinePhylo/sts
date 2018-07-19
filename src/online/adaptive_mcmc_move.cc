//
//  adaptive_mcmc_move.c
//  sts
//
//  Created by Mathieu Fourment on 5/06/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "adaptive_mcmc_move.h"
#include "composite_tree_likelihood.h"
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
	
	AdaptiveMCMCMove::AdaptiveMCMCMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator,
									   gsl_matrix& L,
									   gsl_vector& mu,
									   std::vector<std::unique_ptr<Transform>>& transforms,
									   const double lambda) :
			OnlineMCMCMove(calculator, {}, lambda),
			_L(*gsl_matrix_alloc(L.size1, L.size2)),
			_mu(*gsl_vector_alloc(mu.size)),
			_transforms(transforms){
		for(int i = 0; i < calculator.size(); i++){
			_result.push_back(gsl_vector_alloc(mu.size));
			_work.push_back(gsl_vector_alloc(mu.size));
		}
		gsl_vector_memcpy(&_mu, &mu);
		gsl_matrix_memcpy(&_L, &L);
	}
	
	AdaptiveMCMCMove::AdaptiveMCMCMove(const AdaptiveMCMCMove& adapt):
			OnlineMCMCMove(adapt._calculator, adapt._parameters, adapt._lambda),
			_L(*gsl_matrix_alloc(adapt._L.size1, adapt._L.size2)),
			_mu(*gsl_vector_alloc(adapt._mu.size)),
			_transforms(adapt._transforms){
		for(int i = 0; i < _calculator.size(); i++){
			_result.push_back(gsl_vector_alloc(_mu.size));
			_work.push_back(gsl_vector_alloc(_mu.size));
		}
		gsl_vector_memcpy(&_mu, &adapt._mu);
		gsl_matrix_memcpy(&_L, &adapt._L);
	}
	
	AdaptiveMCMCMove::~AdaptiveMCMCMove()
	{
		// Debug bits
		if(n_attempted > 0) {
			std::clog << "Adaptive_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << " " <<_lambda << std::endl;
		}
		for(int i = 0; i < _result.size(); i++){
			gsl_vector_free(_result[i]);
			gsl_vector_free(_work[i]);
		}
		gsl_vector_free(&_mu);
		gsl_matrix_free(&_L);
	}
	
	int AdaptiveMCMCMove::proposeMove(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
	{
		// Choose an edge at random
		TreeParticle* value = particle.GetValuePointer();
		return proposeMove(*value, rng);
	}
	
	int AdaptiveMCMCMove::proposeMove(TreeParticle& particle, smc::rng* rng){
		
		size_t indexCalculator = 0;
#if defined(_OPENMP)
		indexCalculator = omp_get_thread_num();
#endif
		gsl_vector* work = _work[indexCalculator];
		gsl_vector* result = _result[indexCalculator];
		double orig_ll = particle.logP;
		
		size_t dim = _mu.size;
		gsl_ran_multivariate_gaussian(rng->GetRaw(), &_mu, &_L, result);
		std::map<std::string, double> x;
		double logJacobian = 0;
		
		double new_logQ = 0;
		gsl_ran_multivariate_gaussian_log_pdf(_result[indexCalculator], &_mu, &_L, &new_logQ, work);
		
		int index = 0;
		bool frequencies = false;
		for (std::unique_ptr<Transform>& transform : _transforms) {
			if(dynamic_cast<const SimplexTransform*>(transform.get()) != nullptr){
				const Vdouble& orig_values = particle.model->getFrequencies();
				std::vector<double> transformed_orig_values = transform->transform(orig_values);
				
				std::vector<double> transformed_values;
				size_t temp = index+3;
				for(; index < temp; index++){
					transformed_values.push_back(gsl_vector_get(result, index));
				}
				gsl_vector_set(result, index-3, transformed_orig_values[0]);
				gsl_vector_set(result, index-2, transformed_orig_values[1]);
				gsl_vector_set(result, index-1, transformed_orig_values[2]);
				
				double logJac = 0;
				std::vector<double> values = transform->inverse_transform(transformed_values, &logJac);
				logJacobian -= logJac;
				x["theta"] = values[1] + values[2];
				x["theta1"] = values[0]/(values[0] + values[3]);
				x["theta2"] = values[2]/(values[1] + values[2]);
				
				logJacobian += transform->logJacobian(orig_values);
			}
			else{
				std::string name = transform->getNames()[0];
				double transformed_value = gsl_vector_get(result, index);
				double value = transform->inverse_transform(transformed_value);
				logJacobian -= transform->logJacobian(value);
				x[name] = value;
				
				double orig_value;
				if (particle.model->hasParameter(name)) {
					orig_value = particle.model->getParameterValue(name);
				}
				else if (particle.rateDist->hasParameter(name)) {
					orig_value = particle.rateDist->getParameterValue(name);
				}
				logJacobian += transform->logJacobian(orig_value);
				const double transformed_orig_value = transform->transform(orig_value);
				gsl_vector_set(result, index, transformed_orig_value);

				index++;
			}
		}
		
		double logQ = 0;
		gsl_ran_multivariate_gaussian_log_pdf(result, &_mu, &_L, &logQ, work);
		
		std::vector<double> backup;
		for(const auto& val : x){
			if (particle.model->hasParameter(val.first)) {
				double value = particle.model->getParameterValue(val.first);
				backup.push_back(value);
				particle.model->setParameterValue(val.first, val.second);
			}
			else if (particle.rateDist->hasParameter(val.first)) {
				double value = particle.rateDist->getParameterValue(val.first);
				backup.push_back(value);
				particle.rateDist->setParameterValue(val.first, val.second);
			}
			else{
				std::cerr << val.first << " " << " not found"<< std::endl;
				exit(10);
			}
		}
		
		
		_calculator[indexCalculator]->initialize(*particle.model, *particle.rateDist, *particle.tree);
		
		double new_ll = _calculator[indexCalculator]->operator()();
		
		particle.logP = new_ll;
		
		double mh_ratio = std::exp(new_ll + logQ - orig_ll - new_logQ + logJacobian);
		if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
			return 1;
		} else {
			// Rejected
			particle.logP = orig_ll;
			size_t index = 0;
			for(const auto& val : x){
				if (particle.model->hasParameter(val.first)) {
					particle.model->setParameterValue(val.first, backup[index]);
				}
				else if (particle.rateDist->hasParameter(val.first)) {
					particle.rateDist->setParameterValue(val.first, backup[index]);
				}
				index++;
			}
			return 0;
		}
		
	}
	
}}
