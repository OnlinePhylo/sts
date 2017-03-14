#include "composite_tree_likelihood.h"
#include "tripod_optimizer.h"

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>

using namespace std;
using namespace bpp;

namespace sts { namespace online {

CompositeTreeLikelihood::CompositeTreeLikelihood(std::shared_ptr<FlexibleTreeLikelihood> calculator) :
    calculator_(calculator),
    tree_(nullptr)
{}

CompositeTreeLikelihood::CompositeTreeLikelihood(std::shared_ptr<FlexibleTreeLikelihood> calculator,
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

double CompositeTreeLikelihood::operator()(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength)
{
    assert(tree_ != nullptr && "Uninitialized tree!");
    return calculator_->calculateLogLikelihood(distal, taxonName, pendantLength, distalLength, proximalLength);// + sumAdditionalLogLikes();
}
    
double CompositeTreeLikelihood::logLikelihood()
{
    return this->operator()();
}

double CompositeTreeLikelihood::logLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength)
{
    return this->operator()(distal, taxonName, pendantLength, distalLength, proximalLength);
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

void CompositeTreeLikelihood::calculatePendantDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2){
    return calculator_->calculatePendantDerivatives(distal, taxonName, pendantLength, distalLength, proximalLength, d1, d2);
}

void CompositeTreeLikelihood::calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2){
    return calculator_->calculateDistalDerivatives(distal, taxonName, pendantLength, distalLength, proximalLength, d1, d2);
}
}} // Namespaces
