#include "composite_tree_likelihood.h"
#include "beagle_tree_likelihood.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>

using namespace std;
using namespace bpp;

namespace sts { namespace online {

CompositeTreeLikelihood::CompositeTreeLikelihood(shared_ptr<BeagleTreeLikelihood> calculator) :
    calculator_(calculator),
    tree(nullptr)
{};

CompositeTreeLikelihood::CompositeTreeLikelihood(shared_ptr<BeagleTreeLikelihood> calculator,
                                                     vector<TreeLogLikelihood> additional_log_likes) :
    calculator_(calculator),
    additional_log_likes(additional_log_likes),
    tree(nullptr)
{};

double CompositeTreeLikelihood::operator()()
{
    assert(tree != nullptr);
    double tree_likelihood = calculator_->calculate_log_likelihood();

    for(TreeLogLikelihood& like : additional_log_likes)
        tree_likelihood += like(*tree);

    return tree_likelihood;
}

double CompositeTreeLikelihood::log_likelihood()
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
    this->additional_log_likes.push_back(like);
}

}} // Namespaces
