#include "tree_particle.h"

using namespace bpp;

namespace sts
{

Tree_particle::Tree_particle(SubstitutionModel* model,
                             TreeTemplate<bpp::Node>* tree,
                             DiscreteDistribution* rate_dist,
                             SiteContainer const* sites) :
    model(model),
    tree(tree),
    rate_dist(rate_dist),
    sites(sites)
{}

Tree_particle::Tree_particle(const Tree_particle& other) :
    model(other.model->clone()),
    tree(other.tree->clone()),
    rate_dist(other.rate_dist->clone()),
    sites(other.sites)
{
}

Tree_particle& Tree_particle::operator=(const Tree_particle& other)
{
    sites = other.sites;
    model.reset(other.model->clone());
    tree.reset(other.tree->clone());
    rate_dist.reset(other.rate_dist->clone());
    return *this;
}

}
