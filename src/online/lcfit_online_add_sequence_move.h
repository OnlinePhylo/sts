#ifndef STS_MOVES_LCFIT_ADD_SEQUENCE_MOVE_H
#define STS_MOVES_LCFIT_ADD_SEQUENCE_MOVE_H

#include "guided_online_add_sequence_move.h"


namespace sts { namespace online {

/// \brief Adds a taxon to a tree.
///
/// Adds a taxon, \f$s\f$ to a particle using a pplacer-style proposal density on edges:
///
/// - Two passes are made over the tree; likelihood vectors are cached for both sides of every edge
/// - We calculate the log-likelihood of attaching \f$s\f$ at the midpoint of every edge (linear in number of edges,
///   with cached vectors), for every branch length in #proposePendantBranchLengths.
/// - Mid-edge LnLs form proposal density on edges
/// - Some ML-optimization is performed for attachment location on edge, pendant branch lengths; a branch length is
///   proposed from a distribution around these lengths.
///
/// For each time in <c>[1,taxa_to_add.size()]</c>, adds <c>taxa_to_add[time-1]</c> to the tree.
///
/// **NOTE**: All of the sequences in taxa_to_add must have partials registered in calculator.
///
/// The particle weight is updated using:
/// \f[
/// w_{r,k} = \frac{\gamma*(s_{r,k})}{\gamma*(s_{r-1,k})q(s_{r-1,k}\rightarrow s_{r,k})}
/// \f]
/// Where \f$\gamma*\f$ is the log-likelihood of the tree with \f$n-1\f$ and \f$n\f$ taxa, and the proposal density
/// \f$q(s_{r-1,k}\rightarrow s_{r,k})\f$ uses the mid-edge log-likelihood
class LcfitOnlineAddSequenceMove : public GuidedOnlineAddSequenceMove
{
public:
    /// Constructor
    ///
    /// \param calculator Likelihood calculator
    /// \param taxaToAdd Names of sequences to add, in order
    /// \param proposePendantBranchLengths pendant branch lenghts to attempt attachment with.
	LcfitOnlineAddSequenceMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator,
                               const std::vector<std::string>& sequenceNames,
                               const std::vector<std::string>& taxaToAdd,
                               const std::vector<double>& proposePendantBranchLengths,
                               const double maxLength,
                               const size_t subdivideTop,
                               const double expPriorMean);

    virtual ~LcfitOnlineAddSequenceMove();

protected:
    std::pair<double, double> proposeDistal(bpp::Node& n, const std::string& leafName, const double mlDistal, const double mlPendant, smc::rng* rng) const;
    
    std::pair<double, double> proposePendant(bpp::Node& n, const std::string& leafName, const double mlPendant, const double distalBranchLength, smc::rng* rng) const;


private:
    const std::tuple<bpp::Node*, double, double, double> chooseMoveLocation(bpp::TreeTemplate<bpp::Node>& tree,
                                                                           const std::string& leafName,
                                                                           smc::rng* rng,
                                                                           size_t particleID);

    /// Branch lengths to propose from
    std::vector<double> proposePendantBranchLengths;

    double expPriorMean_;
    mutable size_t lcfit_failures_;
    mutable size_t lcfit_attempts_;
};

}} // namespaces

#endif // STS_MOVES_LCFIT_ADD_SEQUENCE_MOVE_H
