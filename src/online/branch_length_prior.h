#ifndef STS_ONLINE_EXPONENTIAL_BRANCH_LENGTH_PRIOR_H
#define STS_ONLINE_EXPONENTIAL_BRANCH_LENGTH_PRIOR_H

#include <Bpp/Phyl/TreeTemplate.h>
#include <functional>

namespace sts { namespace online {

/// \brief Generic prior density on branch lengths
class BranchLengthPrior
{
public:
    /// \brief Construct with a prior density function
    ///
    /// \param logPriorDensity a function of one parameter, which given a branch length, returns the prior density.
    BranchLengthPrior(std::function<double(double)> logPriorDensity);

    double operator()(const bpp::TreeTemplate<bpp::Node>& tree) const;
private:
    std::function<double(double)> logPriorDensity;
};

}}  // Namespaces

#endif
