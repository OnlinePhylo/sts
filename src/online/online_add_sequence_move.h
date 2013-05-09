#ifndef STS_MOVES_ADD_SEQUENCE_MOVE_H
#define STS_MOVES_ADD_SEQUENCE_MOVE_H

#include <smctc.hh>

#include <forward_list>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <Bpp/Phyl/TreeTemplate.h>

namespace sts { namespace online {

// Forwards
class TreeParticle;
class CompositeTreeLikelihood;

/// \brief Adds a taxon to a tree.
///
/// Adds a taxon, \f$s\f$ to a random edge, drawing a branch length from the proposal distribution
class OnlineAddSequenceMove
{
public:
    /// Constructor
    ///
    /// \param calculator Likelihood calculator
    /// \param taxaToAdd Names of sequences to add, in order
    OnlineAddSequenceMove(CompositeTreeLikelihood& calculator,
                          std::function<std::pair<double,double>(smc::rng*)> branchLengthProposer,
                          const std::vector<std::string>& taxaToAdd);

    void operator()(long, smc::particle<TreeParticle>&, smc::rng*);
protected:
    CompositeTreeLikelihood& calculator;
    std::forward_list<std::string> taxaToAdd;
    std::function<std::pair<double,double>(smc::rng*)> branchLengthProposer;

    long lastTime;
};

}} // namespaces

#endif // STS_MOVES_ADD_SEQUENCE_MOVE_H
