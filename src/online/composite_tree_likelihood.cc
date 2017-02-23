#include "composite_tree_likelihood.h"
#include "beagle_tree_likelihood.h"
#include "tripod_optimizer.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>

using namespace std;
using namespace bpp;

namespace sts { namespace online {

CompositeTreeLikelihood::CompositeTreeLikelihood(std::shared_ptr<BeagleTreeLikelihood> calculator) :
    calculator_(calculator),
    tree_(nullptr)
{}

CompositeTreeLikelihood::CompositeTreeLikelihood(std::shared_ptr<BeagleTreeLikelihood> calculator,
                                                 std::vector<TreeLogLikelihood> additionalLogLikes) :
    calculator_(calculator),
    additionalLogLikes_(additionalLogLikes),
    tree_(nullptr)
{}

double CompositeTreeLikelihood::operator()()
{
    assert(tree_ != nullptr && "Uninitialized tree!");
    return calculator_->calculateLogLikelihood() + sumAdditionalLogLikes();
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
    this->tree_ = &tree;
}

void CompositeTreeLikelihood::add(TreeLogLikelihood like)
{
    this->additionalLogLikes_.push_back(like);
}

double CompositeTreeLikelihood::sumAdditionalLogLikes() const
{
    double sum = 0.0;
    for (const TreeLogLikelihood& like : additionalLogLikes_) {
        sum += like(*tree_);
    }
    return sum;
}

std::vector<std::vector<double>> CompositeTreeLikelihood::calculateAttachmentLikelihoods(const std::string& leafName,
                                                                                         const std::vector<AttachmentLocation>& attachmentLocations,
                                                                                         const std::vector<double> pendantBranchLengths)
{
    return calculator_->calculateAttachmentLikelihoods(leafName, attachmentLocations, pendantBranchLengths);
}

std::vector<double> CompositeTreeLikelihood::calculateAttachmentLikelihood(const std::string& leafName,
                                                                           const bpp::Node* node,
                                                                           const double distal,
                                                                           const std::vector<double> pendantBranchLengths)
{
    return calculator_->calculateAttachmentLikelihood(leafName, node, distal, pendantBranchLengths);
}

}} // Namespaces
