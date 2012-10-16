/// \file phylomoves.hh
/// \author metatangle, inc.
/// \brief Particle move functors for use with SMCTC.
#ifndef __PHYLOMOVES_HH__
#define __PHYLOMOVES_HH__

#include "smctc.hh"
#include "phylofunc.hh"
#include "hmsbeagle.hh"

/// Class to calculate the likelihood of a forest
class forest_likelihood
{
public:
    /// Constructor

    ///  \param calc Initialized likelihood calculator
    ///  \param leaf_nodes Vector representing \\perp
    explicit forest_likelihood(std::shared_ptr<OnlineCalculator> calc,
            std::vector<std::shared_ptr<phylo_node>> leaf_nodes) : calc(calc), leaf_nodes(leaf_nodes) {};
    /// Copy constructor
    explicit forest_likelihood(const forest_likelihood &other) : calc(other.calc), leaf_nodes(other.leaf_nodes) {};

    double operator()(const particle&) const;

    const std::vector<std::shared_ptr<phylo_node>> get_leaves() const;
    std::shared_ptr<OnlineCalculator> get_calculator() const;
private:
    std::shared_ptr<OnlineCalculator> calc;
    std::vector<std::shared_ptr<phylo_node>> leaf_nodes;
};

/// Abstract class for MCMC moves in smctc.
class mcmc_move
{
public:
    /// Create an mcmc_move
    ///  \param log_likelihood forest_likelihood to use for likelihood calculations
    explicit mcmc_move(forest_likelihood& log_likelihood) : log_likelihood(log_likelihood), attempted(0), accepted(0) {};
    /// Number of attempted moves
    unsigned int attempted;
    /// Number of accepted moves
    unsigned int accepted;

    int operator()(long time, smc::particle<particle>& from, smc::rng *rng);

    /// Override in subclass with MCMC move
    virtual int do_move(long, smc::particle<particle>&, smc::rng*) const = 0;

protected:
    forest_likelihood log_likelihood;
};

/// An MCMC move which perturbs branch lengths uniformly from -amount to amount
class uniform_bl_mcmc_move : public mcmc_move
{
public:
    uniform_bl_mcmc_move(forest_likelihood& log_likelihood) : mcmc_move(log_likelihood), amount(0.1) {};
    uniform_bl_mcmc_move(forest_likelihood& log_likelihood, double amount) : mcmc_move(log_likelihood), amount(amount) {};

    int do_move(long, smc::particle<particle>&, smc::rng*) const;

    /// Amount to perturb branch lengths
    double amount;
};

/// A sequential monte carlo move
class smc_move
{
public:
    explicit smc_move(forest_likelihood& log_likelihood) : log_likelihood(log_likelihood) {};

    /// Function call for use with smctc - calls user-defined do_move, tracks result.
    int operator()(long, smc::particle<particle>&, smc::rng*);

    /// Override in subclass with SMC move
    virtual int do_move(long, smc::particle<particle>&, smc::rng*) const = 0;

    /// Number of times this move has been performed
    int call_count;
protected:
    forest_likelihood log_likelihood;
};

/// Merge of two nodes, with exponential branch length proposal
class rooted_merge: public smc_move
{
public:
    explicit rooted_merge(forest_likelihood& log_likelihood) : smc_move(log_likelihood) {};
    int do_move(long, smc::particle<particle>&, smc::rng*) const;
};

/// Particle initializer.

/// Uses as log_likelihood the \\perp ll
class smc_init
{
public:
    explicit smc_init(forest_likelihood& log_likelihood) : log_likelihood(log_likelihood) {};
    virtual smc::particle<particle> operator()(smc::rng*);
protected:
    forest_likelihood log_likelihood;
};
#endif // __PHYLOMOVES_HH__
