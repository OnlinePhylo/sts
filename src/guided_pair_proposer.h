#ifndef STS_GUIDED_PAIR_PROPOSER_H
#define STS_GUIDED_PAIR_PROPOSER_H

#include <string>
#include <vector>
#include "smctc.hh"
#include <forest_likelihood.h>
#include <uniform_pair_proposer.h>

namespace sts
{
namespace moves
{

/// \class Guided_pair_proposer
/// \brief "Propose" pairs of nodes to join using a guide tree.
class Guided_pair_proposer
{
public:
    /// Instantiate 
    /// \c tree_file the file from which to read the guide tree.
    explicit Guided_pair_proposer(int max_tries, sts::likelihood::Forest_likelihood&);
    ~Guided_pair_proposer();
    void initialize(const std::string& tree_file);
    void operator()(particle::Particle, smc::rng*, particle::Node_ptr& a, particle::Node_ptr& b, double& fwd_density, double& back_density);    

private:
    int max_tries;
    std::string file;
    std::vector< std::pair< double, std::pair< int, int > > > sampling_dist;
    double cumulative;
    gsl_rng* r;
    sts::likelihood::Forest_likelihood& log_likelihood;
    Uniform_pair_proposer upp;

    void propose(smc::rng *rng, int& leaf1, int& leaf2, double& density);
};


} // namespace moves
} // namespace sts

#endif // STS_GUIDED_PAIR_PROPOSER_H
