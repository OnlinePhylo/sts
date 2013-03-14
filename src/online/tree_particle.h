#ifndef STS_TREE_PARTICLE_H
#define STS_TREE_PARTICLE_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include <memory>

namespace sts { namespace online {

/// \brief A particle representing a fully-specified tree.
class Tree_particle
{
public:
    /// \brief Constructor
    ///
    /// \param model *owned* Substitution model
    /// \param tree *owned* tree
    /// \param rate_dist *owned* rate distribution
    /// \param sites Site container
    Tree_particle(bpp::SubstitutionModel* model,
                  bpp::TreeTemplate<bpp::Node>* tree,
                  bpp::DiscreteDistribution* rate_dist,
                  bpp::SiteContainer const* sites);
    Tree_particle();
    Tree_particle(const Tree_particle& other);
    Tree_particle& operator=(const Tree_particle& other);
    Tree_particle& operator=(Tree_particle&&) & = default;
    Tree_particle(Tree_particle&&) = default;
    virtual ~Tree_particle() {};

    std::unique_ptr<bpp::SubstitutionModel> model;
    std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tree;
    std::unique_ptr<bpp::DiscreteDistribution> rate_dist;
    bpp::SiteContainer const* sites;
private:
};

}} // Namespaces

#endif
