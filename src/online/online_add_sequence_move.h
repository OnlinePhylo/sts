#ifndef STS_MOVES_ADD_SEQUENCE_MOVE
#define STS_MOVES_ADD_SEQUENCE_MOVE

#include <smctc.hh>

#include <string>
#include <vector>
#include <Bpp/Phyl/TreeTemplate.h>

namespace sts { namespace online {

// Forwards
class Tree_particle;
class Beagle_tree_likelihood;

struct Branch_lengths
{
    double distal_bl;
    double pendant_bl;
};

/// \brief Adds a taxon to a tree.
///
/// For each time in <c>[1,taxa_to_add.size()]</c>, adds <c>taxa_to_add[time-1]</c> to the tree.

/// **NOTE**: All of the sequences in taxa_to_add must have partials registered in calculator.
class Online_add_sequence_move
{
public:
    /// Constructor
    Online_add_sequence_move(Beagle_tree_likelihood& calculator,
                             const std::vector<std::string>& taxa_to_add);

    /// Choose edge on which to insert sequence \c leaf_name
    ///
    /// \param tree
    /// \param leaf_name Name of the new taxon, already registered with the calculator.
    /// \param rng Random number generator
    /// \returns a pair consisting of the node to insert above, and an unnormalized log-likelihood of proposing the node
    /// (forward proposal density)
    std::pair<bpp::Node*, double> choose_edge(bpp::TreeTemplate<bpp::Node>& tree,
                                              const std::string& leaf_name,
                                              smc::rng* rng);
    int operator()(long, smc::particle<Tree_particle>&, smc::rng*);


    Branch_lengths propose_branch_lengths(const bpp::Node* insert_edge, const std::string& new_leaf_name);

protected:
    Beagle_tree_likelihood& calculator;
    std::vector<std::string> taxa_to_add;
};
}} // namespaces

#endif // STS_MOVES_ADD_SEQUENCE_MOVE
