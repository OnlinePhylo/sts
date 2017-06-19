#include "tree_particle.h"

namespace sts { namespace online {

TreeParticle::TreeParticle() :
    model(nullptr),
    tree(nullptr),
    rateDist(nullptr),
    sites(nullptr)
{}

TreeParticle::TreeParticle(std::unique_ptr<bpp::SubstitutionModel> model,
                           std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tree,
                           std::unique_ptr<bpp::DiscreteDistribution> rateDist,
                           bpp::SiteContainer const* sites) :
    model(std::move(model)),
    tree(std::move(tree)),
    rateDist(std::move(rateDist)),
    sites(sites),
    logP(std::numeric_limits<double>::infinity())
{}

TreeParticle::TreeParticle(const TreeParticle& other) :
    model(other.model->clone()),
    tree(other.tree->clone()),
    rateDist(other.rateDist->clone()),
    sites(other.sites),
    particleID(other.particleID),
    logP(std::numeric_limits<double>::infinity())
{
}

TreeParticle& TreeParticle::operator=(const TreeParticle& other)
{
    sites = other.sites;
    model.reset(other.model->clone());
    tree.reset(other.tree->clone());
    rateDist.reset(other.rateDist->clone());
    particleID = other.particleID;
    logP = other.logP;
    return *this;
}

TreeParticle& TreeParticle::operator=(TreeParticle&& other)
{
    sites = other.sites;
    model = std::move(other.model);
    tree = std::move(other.tree);
    rateDist = std::move(other.rateDist);
    particleID = other.particleID;
    logP = other.logP;
    return *this;
}

TreeParticle::TreeParticle(TreeParticle&& other) :
    model(std::move(other.model)),
    tree(std::move(other.tree)),
    rateDist(std::move(other.rateDist)),
    sites(other.sites),
    particleID(other.particleID),
    logP(std::numeric_limits<double>::infinity())
{}


}} // Namespaces
