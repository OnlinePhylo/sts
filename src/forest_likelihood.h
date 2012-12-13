#ifndef STS_LIKELIHOOD_FOREST_LIKELIHOOD_H
#define STS_LIKELIHOOD_FOREST_LIKELIHOOD_H
#include "node_ptr.h"
#include "particle.h"

#include <memory>
#include <vector>

namespace sts
{
namespace likelihood
{

class Online_calculator;

/// Class to calculate the likelihood of a forest
class Forest_likelihood
{
public:
    /// Constructor

    ///  \param calc Initialized likelihood calculator
    ///  \param leaf_nodes Vector representing \\perp
    explicit Forest_likelihood(std::shared_ptr<Online_calculator> calc,
                               std::vector<particle::Node_ptr> leaf_nodes) : calc(calc), leaf_nodes(leaf_nodes) {};
    /// Copy constructor
    explicit Forest_likelihood(const Forest_likelihood &other) : calc(other.calc), leaf_nodes(other.leaf_nodes) {};

    double operator()(const particle::Particle&) const;
    double calculate_log_likelihood(const particle::Particle&) const;

    const std::vector<particle::Node_ptr> get_leaves() const;
    std::shared_ptr<Online_calculator> get_calculator() const;
private:
    std::shared_ptr<Online_calculator> calc;
    std::vector<particle::Node_ptr> leaf_nodes;
};
} // namespace likelihood
} // namespace sts

#endif // STS_LIKELIHOOD_FOREST_LIKELIHOOD_H
