#ifndef STS_MOVES_GUIDED_ADD_SEQUENCE_MOVE_H
#define STS_MOVES_GUIDED_ADD_SEQUENCE_MOVE_H

#include "online_add_sequence_move.h"
#include "tripod_optimizer.h"

#include <limits>

namespace sts {
namespace online {

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
class GuidedOnlineAddSequenceMove : public OnlineAddSequenceMove
{
public:
    /// Constructor
    ///
    /// \param calculator Likelihood calculator
    /// \param taxaToAdd Names of sequences to add, in order
    /// \param proposePendantBranchLengths pendant branch lenghts to attempt attachment with.
    /// \param maxLength Maximum length to allow prior to proposing an additional attachment location
    /// \param subdivideTop Further subdivide the top `N` edges before proposing
    GuidedOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                const std::vector<std::string>& taxaToAdd,
                                const std::vector<double>& proposePendantBranchLengths = std::vector<double>(1, 0.0),
                                const double maxLength = std::numeric_limits<double>::max(),
                                const size_t subdivideTop = 0);
protected:
    /// \brief Propose a distal branch length around the ML value
    ///
    /// \param node Chosen node
    /// \param mlDistal Maximum-likelihood estimate of distal branch length
    /// \return A pair of (distal branch length, log density)
    virtual std::pair<double, double> proposeDistal(const double edgeLength, const double mlDistal, smc::rng* rng) const;

    /// \brief Propose a pendant branch length around the ML value
    ///
    /// \return A pair of (pendant branch length, log density)
    virtual std::pair<double, double> proposePendant(const double mlPendant, smc::rng* rng) const;
    virtual AttachmentProposal propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng);

    /// Choose edge on which to insert sequence \c leaf_name
    ///
    /// \param tree
    /// \param leafName Name of the new taxon, already registered with the calculator.
    /// \param rng Random number generator
    /// \returns a pair consisting of the node to insert above, and an unnormalized log-likelihood of proposing the node
    /// (forward proposal density)
    virtual const std::pair<bpp::Node*, double> chooseEdge(bpp::TreeTemplate<bpp::Node>& tree,
                                                           const std::string& leafName,
                                                           smc::rng* rng, size_t particleID);

    virtual TripodOptimizer optimizeBranchLengths(const bpp::Node* insertEdge, const std::string& newLeafName,
                                                  double& distalBranchLength, double& pendantBranchLength);
private:
    std::vector<std::pair<bpp::Node*, double> > subdivideTopN(std::vector<AttachmentLocation> locs,
                                                         const std::vector<double>& logWeights,
                                                         const std::string& leafName);
    /// Branch lengths to propose from
    std::vector<double> proposePendantBranchLengths;
    double maxLength;
    size_t subdivideTop;
};

}
} // namespaces

#endif // STS_MOVES_GUIDED_ADD_SEQUENCE_MOVE_H
