//
//  laplus_online_add_sequence_move.hpp
//  sts
//
//  Created by Mathieu Fourment on 11/10/18.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef laplus_online_add_sequence_move_hpp
#define laplus_online_add_sequence_move_hpp

#include "guided_online_add_sequence_move.h"

namespace sts { namespace online {
	
	class LaplusOnlineAddSequenceMove : public GuidedOnlineAddSequenceMove
	{
	public:
		/// Constructor
		///
		/// \param calculator Likelihood calculator
		/// \param taxaToAdd Names of sequences to add, in order
		/// \param proposePendantBranchLengths pendant branch lenghts to attempt attachment with.
		LaplusOnlineAddSequenceMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator,
									const std::vector<std::string>& sequenceNames,
									const std::vector<std::string>& taxaToAdd,
									const std::vector<double>& proposePendantBranchLengths,
									const double maxLength,
									const size_t subdivideTop,
									const double expPriorMean);
		
		virtual ~LaplusOnlineAddSequenceMove();
		
	protected:
		std::pair<double, double> proposeDistal(bpp::Node& n, const std::string& leafName, const double mlDistal, const double mlPendant, smc::rng* rng) const;
		
		std::pair<double, double> proposePendant(bpp::Node& n, const std::string& leafName, const double mlPendant, const double distalBranchLength, smc::rng* rng) const;
		
		
	private:
		const std::tuple<bpp::Node*, double, double, double> chooseMoveLocation(bpp::TreeTemplate<bpp::Node>& tree,
																				const std::string& leafName,
																				smc::rng* rng,
																				size_t particleID);
		
		mutable size_t lcfit_failures_;
		mutable size_t lcfit_attempts_;
	};
	
}} // namespaces

#endif /* laplus_online_add_sequence_move_hpp */
