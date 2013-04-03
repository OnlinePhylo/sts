#include "composite_tree_likelihood.h"
#include "beagle_tree_likelihood.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>

using namespace std;
using namespace bpp;

namespace sts { namespace online {

CompositeTreeLikelihood::CompositeTreeLikelihood(std::shared_ptr<BeagleTreeLikelihood> calculator) :
    calculator_(calculator),
    tree(nullptr)
{};

CompositeTreeLikelihood::CompositeTreeLikelihood(std::shared_ptr<BeagleTreeLikelihood> calculator,
                                                 std::vector<TreeLogLikelihood> additionalLogLikes) :
    calculator_(calculator),
    additionalLogLikes(additionalLogLikes),
    tree(nullptr)
{};

double CompositeTreeLikelihood::operator()()
{
    assert(tree != nullptr);
    double tree_likelihood = calculator_->calculateLogLikelihood();

    for(TreeLogLikelihood& like : additionalLogLikes)
        tree_likelihood += like(*tree);

    return tree_likelihood;
}

double CompositeTreeLikelihood::logLikelihood()
{
    return this->operator()();
}

void CompositeTreeLikelihood::initialize(const SubstitutionModel& model,
                                         const DiscreteDistribution& rate_dist,
                                         TreeTemplate<Node>& tree)
{
    calculator_->initialize(model, rate_dist, tree);
    this->tree = &tree;
}

void CompositeTreeLikelihood::add(TreeLogLikelihood like)
{
    this->additionalLogLikes.push_back(like);
}

}} // Namespaces
