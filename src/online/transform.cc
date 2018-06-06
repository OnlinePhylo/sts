//
//  Transform.cc
//  sts
//
//  Created by Mathieu Fourment on 5/06/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "transform.h"

namespace sts {
	namespace online {
		
	
		double LogTransform::logJacobian(double value){
			return -std::log(value);
		}

		double LogTransform::transform(double value){
			return std::log(value);
		}

		double LogTransform::inverse_transform(double value){
			return std::exp(value);
		}
		
		
		// values: N+1 elements (constrained)
		double SimplexTransform::logJacobian(std::vector<double> values){
			size_t tranformed_size = values.size()-1;
			double jacobian = 0;
			for(size_t i = 0; i < tranformed_size; i++){
				double sum = 0.0;
				for(size_t j = 0; j < i; j++){
					sum += values[j];
				}
				jacobian += log(-1.0/((1.0 - sum)*(values[i]*values[i] - values[i])));
			}
			return jacobian;
		}
		
		// values: N+1 elements (constrained)
		// return N elements
		std::vector<double> SimplexTransform::transform(std::vector<double> values){
			size_t tranformed_size = values.size()-1;
			std::vector<double> xx(tranformed_size, 0);
			for (size_t i = 0; i < tranformed_size; i++) {
				double sum = 1;
				for (size_t j = 0; j < i; j++) {
					sum -= values[j];
				}
				xx[i] = sts::util::logit(values[i]/sum) - std::log(1.0/(tranformed_size-i));
			}
			return xx;
		}
		
		// values: N elements
		// return N+1 elements (constrained)
		std::vector<double> SimplexTransform::inverse_transform(const std::vector<double> values, double* logJacobian){
			size_t tranformed_size = values.size();
			std::vector<double> xx(tranformed_size+1, 0);
			*logJacobian = 0;
			xx[tranformed_size] = 1.0;
			for(size_t i = 0; i < tranformed_size; i++){
				double zi = sts::util::logitinv(values[i] + std::log(1.0/(tranformed_size-i)));
				double sum = 0.0;
				for(size_t j = 0; j < i; j++){
					sum += xx[j];
				}
				xx[i] = (1.0-sum)*zi;
				if(logJacobian != nullptr){
					*logJacobian += log(-1.0/((1.0 - sum)*(xx[i]*xx[i] - xx[i])));
				}
				xx[tranformed_size] -= xx[i];
			}
			return xx;
		}
	}
}