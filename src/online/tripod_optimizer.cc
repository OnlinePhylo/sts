#include "tripod_optimizer.h"

#include <memory>

#include "beagle_tree_likelihood.h"
#include "gsl.h"
#include "util.h"

namespace sts { namespace online {

const double TripodOptimizer::TOLERANCE = 1e-3;

TripodOptimizer::TripodOptimizer(std::shared_ptr<BeagleTreeLikelihood> btl, const bpp::Node* insertEdge, const std::string& newLeafName, double d) :
    b1(btl->borrowBuffer()), b2(btl->borrowBuffer())
{
    this->beagleInstance = btl->beagleInstance();
    this->scratch1 = b1.value();
    this->scratch2 = b2.value();
    this->distalBuffer = btl->getDistalBuffer(insertEdge);
    this->proximalBuffer = btl->getProximalBuffer(insertEdge);
    this->leafBuffer = btl->getLeafBuffer(newLeafName);
    this->d = d;
}

double minimize(std::function<double(double)> fn,
                double rawStart,
                double left,
                double right,
                const size_t maxIters=5)
{
    size_t iter = 0;

    double lefty = fn(left);
    double righty = fn(right);
    double start = rawStart;
    double val;
    double min_x = lefty < righty ? left : right;
    double min_y = std::min(righty, lefty);

    for(iter = 0; iter < maxIters; iter++) {
        val = fn(start);
        if(val < min_y)
            return sts::gsl::minimize(fn, start, left, right, maxIters - iter);

        if(std::abs(start - min_x) < TripodOptimizer::TOLERANCE)
            return start;
        start = (start + min_x) / 2;
    }

    return start;

}

/// Optimize distal branch length, keeping pendant fixed
double TripodOptimizer::optimizeDistal(const double distal_start, const double pendant, size_t max_iters)
{
    auto fn = [&](double distal) {
        return - log_like(distal, pendant, true);
    };
    return minimize(fn, distal_start, 0, d, max_iters);
}

/// Optimize pendant branch length, keeping distal fixed
double TripodOptimizer::optimizePendant(const double distal, const double pendant_start, size_t max_iters)
{
    auto fn = [&](double pendant) {
        return -log_like(distal, pendant, false);
    };

    return minimize(fn, pendant_start, 0, 2.0, max_iters);
}

double TripodOptimizer::log_like(const double distal, const double pendant, const bool distal_changed)
{
    std::vector<BeagleOperation> operations;
    std::vector<double> branch_lengths;
    std::vector<int> node_indices;

    // If distal changed, update partial
    if(distal_changed) {
        operations.push_back(BeagleOperation({scratch1,
                        BEAGLE_OP_NONE,
                        BEAGLE_OP_NONE,
                        distalBuffer,
                        distalBuffer,
                        proximalBuffer,
                        proximalBuffer}));
        branch_lengths.push_back(distal);
        node_indices.push_back(distalBuffer);
        branch_lengths.push_back(d - distal);
        node_indices.push_back(proximalBuffer);
    }
    // Always update root partials
    operations.push_back(BeagleOperation({scratch2,
                    BEAGLE_OP_NONE,
                    BEAGLE_OP_NONE,
                    scratch1,
                    scratch1,
                    leafBuffer,
                    leafBuffer}));
    branch_lengths.push_back(0);
    node_indices.push_back(scratch1);
    branch_lengths.push_back(pendant);
    node_indices.push_back(leafBuffer);

    // Usual thing

    using sts::util::beagle_check;

    beagle_check(beagleUpdateTransitionMatrices(beagleInstance,
                                                0,
                                                node_indices.data(),
                                                NULL,
                                                NULL,
                                                branch_lengths.data(),
                                                node_indices.size()));
    beagle_check(beagleUpdatePartials(beagleInstance, operations.data(), operations.size(), scratch2));

    std::vector<int> scale_indices(operations.size());
    for(size_t i = 0; i < operations.size(); i++)
        scale_indices[i] = operations[i].destinationPartials;

    beagle_check(beagleAccumulateScaleFactors(beagleInstance, scale_indices.data(), scale_indices.size(),
                                              scratch2));
    const int categoryWeightIdx = 0;
    const int stateFreqIdx = 0;
    double logLike;
    beagle_check(beagleCalculateRootLogLikelihoods(beagleInstance,
                                                   &scratch2,
                                                   &categoryWeightIdx,
                                                   &stateFreqIdx,
                                                   &scratch2,
                                                   1,
                                                   &logLike));
    return logLike;
}

}} // namespace sts::online
