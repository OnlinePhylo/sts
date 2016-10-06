#ifndef STS_ONLINE_UNIFORM_LENGTH_ONLINE_ADD_SEQUENCE_MOVE_H
#define STS_ONLINE_UNIFORM_LENGTH_ONLINE_ADD_SEQUENCE_MOVE_H

#include "online_add_sequence_move.h"

#include <functional>

namespace sts { namespace online {

/// \brief Adds a taxon to a tree.
///
/// Adds a taxon, \f$s\f$ to a random edge, with edges weighted based on length.
/// \f[
/// w_{r,k} = \frac{P(T|D)P(T)}{P(T \setminus t_0|D \setminus t_0)P(T \setminus t_0) q(s_{r-1,k}\rightarrow s_{r,k})}
/// \f]
class UniformLengthOnlineAddSequenceMove : public OnlineAddSequenceMove
{
public:
    /// Constructor
    ///
    /// \param calculator Likelihood calculator
    /// \param branchLengthProposer Branch length proposer - should return (branch length, log density)
    /// \param taxaToAdd Names of sequences to add, in order
    UniformLengthOnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                                       const std::vector<std::string>& taxaToAdd,
                                       std::function<std::pair<double,double>(smc::rng*)> branchLengthProposer);
protected:
    virtual AttachmentProposal propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng);
    std::function<std::pair<double,double>(smc::rng*)> branchLengthProposer;
};

}} // namespaces

#endif // STS_ONLINE_UNIFORM_LENGTH_ONLINE_ADD_SEQUENCE_MOVE_H
