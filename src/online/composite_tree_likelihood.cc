#include "composite_tree_likelihood.h"
#include "beagle_tree_likelihood.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>

using namespace std;
using namespace bpp;

namespace sts { namespace online {

Composite_tree_likelihood::Composite_tree_likelihood(shared_ptr<Beagle_tree_likelihood> calculator) :
    calculator_(calculator),
    tree(nullptr)
{};

Composite_tree_likelihood::Composite_tree_likelihood(shared_ptr<Beagle_tree_likelihood> calculator,
                                                     vector<Tree_log_likelihood> additional_log_likes) :
    calculator_(calculator),
    additional_log_likes(additional_log_likes),
    tree(nullptr)
{};

double Composite_tree_likelihood::operator()()
{
    assert(tree != nullptr);
    double tree_likelihood = calculator_->calculate_log_likelihood();

    for(Tree_log_likelihood& like : additional_log_likes)
        tree_likelihood += like(*tree);

    return tree_likelihood;
}

double Composite_tree_likelihood::log_likelihood()
{
    return this->operator()();
}

void Composite_tree_likelihood::initialize(const SubstitutionModel& model,
                                           const DiscreteDistribution& rate_dist,
                                           TreeTemplate<Node>& tree)
{
    calculator_->initialize(model, rate_dist, tree);
    this->tree = &tree;
}

void Composite_tree_likelihood::add(Tree_log_likelihood like)
{
    this->additional_log_likes.push_back(like);
}

}} // Namespaces
