#include "tree_particle.h"

using namespace bpp;

namespace sts { namespace online {

TreeParticle::TreeParticle() :
    model(nullptr),
    tree(nullptr),
    rate_dist(nullptr),
    sites(nullptr)
{}

TreeParticle::TreeParticle(SubstitutionModel* model,
                           TreeTemplate<bpp::Node>* tree,
                           DiscreteDistribution* rate_dist,
                           SiteContainer const* sites) :
    model(model),
    tree(tree),
    rate_dist(rate_dist),
    sites(sites)
{}

TreeParticle::TreeParticle(const TreeParticle& other) :
    model(other.model->clone()),
    tree(other.tree->clone()),
    rate_dist(other.rate_dist->clone()),
    sites(other.sites)
{
}

TreeParticle& TreeParticle::operator=(const TreeParticle& other)
{
    sites = other.sites;
    model.reset(other.model->clone());
    tree.reset(other.tree->clone());
    rate_dist.reset(other.rate_dist->clone());
    return *this;
}

TreeParticle& TreeParticle::operator=(TreeParticle&& other)
{
    sites = other.sites;
    model = std::move(other.model);
    tree = std::move(other.tree);
    rate_dist = std::move(other.rate_dist);
    return *this;
}

TreeParticle::TreeParticle(TreeParticle&& other) :
    model(std::move(other.model)),
    tree(std::move(other.tree)),
    rate_dist(std::move(other.rate_dist)),
    sites(other.sites)
{}


}} // Namespaces
