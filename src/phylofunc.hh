#ifndef __PHYLOFUNC_H__
#define __PHYLOFUNC_H__

#include "smctc.hh"
#include <memory>
#include <string>
#include <vector>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Tree.h>

/// \class phylo_node
/// Represents the merge of two trees in a forest.
class edge;
class online_calculator;

class phylo_node
{
public:
    explicit phylo_node(std::shared_ptr<online_calculator> calc);
    explicit phylo_node(const phylo_node &other);
    ~phylo_node();

    std::shared_ptr<edge> child1;
    std::shared_ptr<edge> child2;

    // convenience for proposals, height must always increase.
    // In the non-clock case, height is the diameter (2 * distance to closest leaf)
    double height;
    int id;	// node id (1..n-1) for leaf nodes, corresponds to index in alignment. n..2n-1 for internal nodes.
    // XXX AD shouldn't this be (0..n-1) or (1..n)?
    bool is_leaf();

    /// Calculate the height once children have been set
    void calc_height();

    /// Make a phylo_node from a bpp Tree
    static std::shared_ptr< phylo_node >
      of_tree(std::shared_ptr< online_calculator >, bpp::TreeTemplate<bpp::Node> &);

    /// Make a phylo_node from a bpp Tree and node number
    static std::shared_ptr< phylo_node >
      of_tree(std::shared_ptr< online_calculator >, bpp::TreeTemplate<bpp::Node> &, int);

private:
    std::weak_ptr<online_calculator> calc;
};

/// An edge
class edge
{
public:
    /// Initialize with a node and distance.
    edge(std::shared_ptr<phylo_node>, double);

    double length;
    std::shared_ptr<phylo_node> node;

    /// Make an edge from a bpp Tree and node number
    static std::shared_ptr< edge > of_tree(std::shared_ptr< online_calculator >, bpp::TreeTemplate<bpp::Node> &, int);
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

    /// Make a phylo_particle from a bpp Tree
    static std::shared_ptr< phylo_particle >
      of_tree(std::shared_ptr< online_calculator >, bpp::TreeTemplate<bpp::Node> &);

    /// Make a phylo_particle from a Newick tree string
    static std::shared_ptr< phylo_particle >
      of_newick_string(std::shared_ptr< online_calculator >, std::string &);
};

/// \class particle

/// A particle in the SMC
typedef std::shared_ptr<phylo_particle> particle;

int tree_count(const std::vector< std::shared_ptr< phylo_node > > &);
std::vector< std::shared_ptr< phylo_node > > uncoalesced_nodes(std::shared_ptr<phylo_particle> pp, std::vector<std::shared_ptr<phylo_node>> leaf_nodes);

void write_tree(std::ostream &out, const std::shared_ptr< phylo_node > root, const std::vector< std::string > &names);

#endif // __PHYLOFUNC_H__
