#ifndef __PHYLOFUNC_H__
#define __PHYLOFUNC_H__

#include "smctc.hh"
#include <memory>
#include <string>

/// \class phylo_node
/// Represents the merge of two trees in a forest.
class phylo_node
{
public:
    phylo_node();
    ~phylo_node();

    std::shared_ptr< phylo_node > child1; // Am I nuts, or would it be nice to have a type that stores the child and the distance both to keep them tidy? Note that the "height" may no longer be needed soon.
    std::shared_ptr< phylo_node > child2;
    double dist1, dist2;
    double height;	// convenience for proposals, height must always increase
    int id;	// node id (1..n-1) for leaf nodes, corresponds to index in alignment. n..2n-1 for internal nodes.
    // XXX AD shouldn't this be (0..n-1) or (1..n)?
    bool is_leaf();
};

/// \class phylo_particle
/// A forest in the SMC.
///
/// This class stores the SMC forest implicitly, by specifying the collections
/// of mergers that must be made in order to get the forest from \perp, the
/// completely un-merged state.
class phylo_particle
{
public:
    // The merge novel to this particle. If NULL then the particle is \perp.
    std::shared_ptr< phylo_node > node;
    // The predecessor particles, which specify the rest of the merges for this particle.
    std::shared_ptr< phylo_particle > predecessor;
};

/// \class phylo_particle
/// A particle in the SMC.
class particle
{
public:
    std::shared_ptr< phylo_particle > pp;
};


double logLikelihood(long lTime, const particle & X);

smc::particle<particle> fInitialise(smc::rng *pRng);
long fSelect(long lTime, const smc::particle<particle>& p, smc::rng *pRng);
void fMove(long lTime, smc::particle<particle>& pFrom, smc::rng *pRng);

int fMoveNodeAgeMCMC(long lTime, smc::particle<particle>& pFrom, smc::rng *pRng);

extern std::vector< std::shared_ptr< phylo_node > > leaf_nodes;
extern std::vector< std::pair< std::string, std::string > > aln;

class mcmc_move {
public:
    mcmc_move() : attempted(0), accepted(0) {};
    /// Number of attempted moves
    unsigned int attempted;
    /// Number of accepted moves
    unsigned int accepted;

    int operator()(long time, smc::particle<particle>& from, smc::rng *rng) {
        attempted++;
        int result = do_move(time, from, rng);
        if(result) accepted++;
        return result;
    }

    virtual int do_move(long, smc::particle<particle>&, smc::rng*) const = 0;
};

/// An MCMC move which perturbs branch lengths uniformly from -amount to amount
class uniform_bl_mcmc_move : public mcmc_move {
public:
    uniform_bl_mcmc_move() : amount(0.1) {};
    uniform_bl_mcmc_move(double amount) : amount(amount) {};

    int do_move(long, smc::particle<particle>&, smc::rng*) const;

    double amount;
};

#endif // __PHYLOFUNC_H__

