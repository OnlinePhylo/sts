#ifndef __PHYLOMOVES_HH__
#define __PHYLOMOVES_HH__

#include "smctc.hh"
#include "phylofunc.hh"
#include "hmsbeagle.hh"

double log_likelihood(OnlineCalculator, long, const particle&);

/// Abstract class for MCMC moves in smctc.
class mcmc_move
{
public:
    /// Create an mcmc_move
    ///  \param calc OnlineCalculator to use for likelihood calculations
    explicit mcmc_move(std::shared_ptr<OnlineCalculator> calc) : calc(calc), attempted(0), accepted(0) {};
    /// Number of attempted moves
    unsigned int attempted;
    /// Number of accepted moves
    unsigned int accepted;

    int operator()(long time, smc::particle<particle>& from, smc::rng *rng);

    /// Override in subclass with MCMC move
    virtual int do_move(long, smc::particle<particle>&, smc::rng*) const = 0;

protected:
    std::shared_ptr<OnlineCalculator> calc;
};

/// An MCMC move which perturbs branch lengths uniformly from -amount to amount
class uniform_bl_mcmc_move : public mcmc_move
{
public:
    uniform_bl_mcmc_move(std::shared_ptr<OnlineCalculator> calc) : mcmc_move(calc), amount(0.1) {};
    uniform_bl_mcmc_move(std::shared_ptr<OnlineCalculator> calc, double amount) : mcmc_move(calc), amount(amount) {};

    int do_move(long, smc::particle<particle>&, smc::rng*) const;

    double amount;
};

class smc_move
{
public:
    explicit smc_move(std::shared_ptr<OnlineCalculator> calc) : calc(calc) {};

    /// Function call for use with smctc - calls user-defined do_move, tracks result.
    int operator()(long, smc::particle<particle>&, smc::rng*);

    /// Override in subclass with SMC move
    virtual int do_move(long, smc::particle<particle>&, smc::rng*) const = 0;

    /// Number of times this move has been performed
    int call_count;

protected:
    std::shared_ptr<OnlineCalculator> calc;
};

class rooted_merge: public smc_move
{
public:
    int do_move(long, smc::particle<particle>&, smc::rng*) const;
};
#endif // __PHYLOMOVES_HH__
