#ifndef STS_GUIDED_PAIR_PROPOSER_H
#define STS_GUIDED_PAIR_PROPOSER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <Bpp/Seq/DistanceMatrix.h>
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
    explicit Guided_pair_proposer(sts::likelihood::Forest_likelihood&);
    ~Guided_pair_proposer();
    void initialize(const std::string& tree_file, const std::unordered_map<particle::Node_ptr, std::string>& node_name_map);
    void operator()(particle::Particle, smc::rng*, particle::Node_ptr& a, particle::Node_ptr& b, double& fwd_density, double& back_density);    

private:
    std::string file;
    double strength;
    gsl_rng* r;
    sts::likelihood::Forest_likelihood& log_likelihood;
    Uniform_pair_proposer upp;
    std::unordered_map<particle::Node_ptr, int> node_dm_id_map;
    std::unordered_map<int, particle::Node_ptr> dm_id_node_map;
    bpp::DistanceMatrix* dm;

    double get_weight(int i, int j) const;
    void propose(smc::rng *rng, int& leaf1, int& leaf2, double& density, std::vector< std::pair< double, std::pair< int, int > > >& distribution);
};


} // namespace moves
} // namespace sts

#endif // STS_GUIDED_PAIR_PROPOSER_H
