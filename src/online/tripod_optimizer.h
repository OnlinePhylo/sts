#ifndef STS_MOVES_TRIPOD_OPTIMIZER_H
#define STS_MOVES_TRIPOD_OPTIMIZER_H

#include "beagle_tree_likelihood.h"
#include <memory>
#include <string>

namespace sts { namespace online {

class TripodOptimizer
{
public:
    const static double TOLERANCE;

    TripodOptimizer(std::shared_ptr<BeagleTreeLikelihood> btl, const bpp::Node* insertEdge, const std::string& newLeafName, double d);

    /// Optimize distal branch length, keeping pendant fixed
    double optimizeDistal(const double distal_start, const double pendant, size_t max_iters=10);
    /// Optimize pendant branch length, keeping distal fixed
    double optimizePendant(const double distal, const double pendant_start, size_t max_iters=10);
    double log_like(const double distal, const double pendant, const bool distal_changed=true);

private:
    BeagleBuffer b1, b2;

    int beagleInstance,
        distalBuffer,
        proximalBuffer,
        leafBuffer,
        scratch1,
        scratch2;
    double d;
};

}} // namespace sts::online

#endif // STS_MOVES_TRIPOD_OPTIMIZER_H
