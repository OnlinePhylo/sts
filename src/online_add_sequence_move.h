#ifndef STS_MOVES_ROOTED_MERGE_H
#define STS_MOVES_ROOTED_MERGE_H

#include <smctc.hh>

#include <string>
#include <vector>

// Forwards
namespace sts {
class Tree_particle;
namespace likelihood {
class Beagle_tree_likelihood;
}
}

namespace sts { namespace moves {

/// \brief Adds a taxon to a tree.
///
/// For each time in <c>[1,taxa_to_add.size()]</c>, adds <c>taxa_to_add[time-1]</c> to the tree.

/// **NOTE**: All of the sequences in taxa_to_add must have partials registered in calculator.
class Online_add_sequence_move
{
public:
    /// Constructor
    Online_add_sequence_move(sts::likelihood::Beagle_tree_likelihood& calculator,
                             const std::vector<std::string>& taxa_to_add);

    int operator()(long, smc::particle<sts::Tree_particle>&, smc::rng*) const;

protected:
    sts::likelihood::Beagle_tree_likelihood& calculator;
    std::vector<std::string> taxa_to_add;
};
}} // namespaces

#endif // STS_MOVES_ROOTED_MERGE_H
