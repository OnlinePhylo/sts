//
//  Transform.h
//  sts
//
//  Created by Mathieu Fourment on 5/06/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef Transform_h
#define Transform_h

#include <cmath>
#include <vector>
#include <assert.h>

#include "util.h"

namespace sts {
	namespace online {
	
		class Transform{
		public:
			Transform(std::vector<std::string> names) : _names(names){}
			
			Transform(const Transform& transform) : _names(transform._names){}
			
			virtual ~Transform(){};
			
			virtual Transform * clone () const = 0;
			
			virtual double logJacobian(double value){
				throw std::runtime_error("logJacobian method should not be called on this object");
			}
			
			virtual double logJacobian(std::vector<double> values){
				throw std::runtime_error("logJacobian method should not be called on this object");
			}
			
			virtual double transform(double value){
				throw std::runtime_error("transform method should not be called on this object");
			}
			
			virtual double inverse_transform(double value){
				throw std::runtime_error("inverse_transform method should not be called on this object");
			}
			
			virtual std::vector<double> transform(std::vector<double> values){
				throw std::runtime_error("transform method should not be called on this object");
			}
			
			virtual std::vector<double> inverse_transform(const std::vector<double> values, double* logJacobian){
				throw std::runtime_error("inverse_transform method should not be called on this object");
			}
			
			const std::vector<std::string>& getNames() const{
				return _names;
			}
		private:
			std::vector<std::string> _names;
		};
		
		class LogTransform : public Transform{
		public:
			LogTransform(std::vector<std::string> names) : Transform(names){}
			
			LogTransform(const LogTransform& transform) : Transform(transform){}
			
			virtual ~LogTransform(){};
			
			LogTransform * clone () const {return new LogTransform(*this);}
			
			virtual double logJacobian(double value);
			
			virtual double transform(double value);
			
			virtual double inverse_transform(double value);
		};
		
		
		class SimplexTransform : public Transform{
		public:
			
			SimplexTransform(std::vector<std::string> names) : Transform(names){}
			
			SimplexTransform(const SimplexTransform& transform) : Transform(transform){}
			
			virtual ~SimplexTransform(){};
			
			SimplexTransform * clone () const {return new SimplexTransform(*this);}
			
			// values: N+1 elements (constrained)
			virtual double logJacobian(std::vector<double> values);
			
			// values: N+1 elements (constrained)
			// return N elements
			virtual std::vector<double> transform(std::vector<double> values);
			
			// values: N elements
			// return N+1 elements (constrained)
			virtual std::vector<double> inverse_transform(const std::vector<double> values, double* logJacobian);
		};
		
	}
}

#endif /* Transform_h */
