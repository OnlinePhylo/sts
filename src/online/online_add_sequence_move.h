#ifndef STS_MOVES_ADD_SEQUENCE_MOVE_H
#define STS_MOVES_ADD_SEQUENCE_MOVE_H

#include <lcfit_cpp.h>
#include <smctc.hh>

#include <forward_list>
#include <string>
#include <utility>
#include <vector>

#include <Bpp/Phyl/TreeTemplate.h>

namespace sts { namespace online {

// Forwards
class TreeParticle;
class CompositeTreeLikelihood;

struct AttachmentProposal
{
    bpp::Node* edge;
    double edgeLogProposalDensity;
    double distalBranchLength;
    double distalLogProposalDensity;
    double pendantBranchLength;
    double pendantLogProposalDensity;

    double mlDistalBranchLength;
    double mlPendantBranchLength;

    bool lcfitFailure;
    lcfit::LCFitResult lcfitResult;

    double logProposalDensity() const { return edgeLogProposalDensity + distalLogProposalDensity + pendantLogProposalDensity; };
};

struct ProposalRecord
{
    long int T;
    double originalLogLike;
    double newLogLike;
    double originalLogWeight;
    double newLogWeight;
    AttachmentProposal proposal;
};

/// \brief Adds a taxon to a tree.
///
/// Adds a taxon, \f$s\f$ to a random edge, drawing a branch length from the proposal distribution
/// \f[
/// w_{r,k} = \frac{P(T|D)P(T)}{P(T \setminus t_0|D \setminus t_0)P(T \setminus t_0) q(s_{r-1,k}\rightarrow s_{r,k})}
/// \f]
class OnlineAddSequenceMove
{
public:
    OnlineAddSequenceMove(CompositeTreeLikelihood& treeLikelihood,
                          const std::vector<std::string>& taxaToAdd);
    virtual ~OnlineAddSequenceMove() {};

    void operator()(long, smc::particle<TreeParticle>&, smc::rng*);

    void addProposalRecord(const ProposalRecord& proposalRecord);
    const std::vector<ProposalRecord> getProposalRecords() const;

protected:
    virtual AttachmentProposal propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng) = 0;

    CompositeTreeLikelihood& calculator;

private:
    std::forward_list<std::string> taxaToAdd;
    long lastTime;
    std::vector<ProposalRecord> proposalRecords_;
};

}} // namespaces

#endif // STS_MOVES_ADD_SEQUENCE_MOVE_H
