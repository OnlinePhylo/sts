#ifndef STS_TREE_PARTICLE_H
#define STS_TREE_PARTICLE_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include <memory>

namespace sts {

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
                  const bpp::SiteContainer& sites) :
        model(model),
        tree(tree),
        rate_dist(rate_dist),
        sites(sites) {};
private:
    std::unique_ptr<bpp::SubstitutionModel> model;
    std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tree;
    std::unique_ptr<bpp::DiscreteDistribution> rate_dist;
    const bpp::SiteContainer& sites;
};

}

#endif
