#ifndef STS_TREE_PARTICLE_H
#define STS_TREE_PARTICLE_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include <memory>

namespace sts { namespace online {

/// \brief A particle representing a fully-specified tree.
class TreeParticle
{
public:
    /// \brief Constructor
    ///
    /// \param model *owned* Substitution model
    /// \param tree *owned* tree
    /// \param rateDist *owned* rate distribution
    /// \param sites Site container
    TreeParticle(std::unique_ptr<bpp::SubstitutionModel> model,
                 std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tree,
                 std::unique_ptr<bpp::DiscreteDistribution> rateDist,
                 bpp::SiteContainer const* sites);
    TreeParticle();
    TreeParticle(const TreeParticle& other);
    TreeParticle& operator=(const TreeParticle& other);
    TreeParticle& operator=(TreeParticle&&);
    TreeParticle(TreeParticle&&);
    virtual ~TreeParticle() {};

    std::unique_ptr<bpp::SubstitutionModel> model;
    std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tree;
    std::unique_ptr<bpp::DiscreteDistribution> rateDist;
    bpp::SiteContainer const* sites;
    size_t particleID;
    double logP;
private:
};

}} // Namespaces

#endif
