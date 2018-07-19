#ifndef STS_MOVES_ADD_SEQUENCE_MOVE_H
#define STS_MOVES_ADD_SEQUENCE_MOVE_H

#include <lcfit_cpp.h>
#include <smctc.hh>

#include <forward_list>
#include <string>
#include <utility>
#include <vector>

#include <Bpp/Phyl/TreeTemplate.h>
#include <lcfit_cpp.h>

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

    std::string proposalMethodName;

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
	OnlineAddSequenceMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& treeLikelihood,
                          const std::vector<std::string>& sequenceNames,
                          const std::vector<std::string>& taxaToAdd);
    virtual ~OnlineAddSequenceMove() {};

    void operator()(long, smc::particle<TreeParticle>&, smc::rng*);

	void addProposalRecord(const ProposalRecord& proposalRecord);
	const std::vector<ProposalRecord> getProposalRecords() const;
	void enableProposalRecords(bool enable){_recordProposals = enable;}
	
    void addTaxa(const std::vector<std::string>& taxaToAdd);
    
protected:
    virtual AttachmentProposal propose(const std::string& leafName, smc::particle<TreeParticle>& particle, smc::rng* rng) = 0;

    std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator;
    
    std::vector<std::string> _sequenceNames;
    std::forward_list<std::string> taxaToAdd;
    
    std::map<size_t, std::vector<std::pair<size_t, double>>> _probs;
    std::unordered_map<size_t, std::unordered_map<size_t, std::pair<double, double>>> _mles;
    size_t _toAddCount;
    size_t _counter;
    std::string _proposalMethodName;
    
private:
    long lastTime;
    std::vector<ProposalRecord> proposalRecords_;
	bool _recordProposals;
	
};

}} // namespaces

#endif // STS_MOVES_ADD_SEQUENCE_MOVE_H
