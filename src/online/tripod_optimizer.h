#ifndef STS_MOVES_TRIPOD_OPTIMIZER_H
#define STS_MOVES_TRIPOD_OPTIMIZER_H

#include <memory>
#include <string>

//#include "attachment_likelihood.h"
#include "composite_tree_likelihood.h"

namespace sts { namespace online {

class TripodOptimizer
{
public:
    const static double TOLERANCE;

    TripodOptimizer(CompositeTreeLikelihood& ctl, const bpp::Node* insertEdge, const std::string& newLeafName, double d);

    virtual ~TripodOptimizer(){}

    /// Optimize distal branch length, keeping pendant fixed
    double optimizeDistal(const double distal_start, const double pendant, size_t max_iters=10);
    /// Optimize pendant branch length, keeping distal fixed
    double optimizePendant(const double distal, const double pendant_start, size_t max_iters=10);
    double logLike(const double distal, const double pendant, const bool distal_changed=true);

private:
    CompositeTreeLikelihood& _ctl;
    const bpp::Node* _insertEdge;
    const std::string& _newLeafName;
    double d;
};

}} // namespace sts::online

#endif // STS_MOVES_TRIPOD_OPTIMIZER_H
