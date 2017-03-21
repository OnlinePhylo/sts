#include "tripod_optimizer.h"

#include <memory>


#include "gsl.h"
#include "util.h"

namespace sts { namespace online {

const double TripodOptimizer::TOLERANCE = 1e-3;

TripodOptimizer::TripodOptimizer(CompositeTreeLikelihood& ctl, const bpp::Node* insertEdge, const std::string& newLeafName, double d) :
        _ctl(ctl), _insertEdge(insertEdge), _newLeafName(newLeafName)
{
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
        return - _ctl(*_insertEdge, _newLeafName, pendant, distal, d-distal);
    };
    return minimize(fn, distal_start, 0, d, max_iters);
}

/// Optimize pendant branch length, keeping distal fixed
double TripodOptimizer::optimizePendant(const double distal, const double pendant_start, size_t max_iters)
{
    auto fn = [&](double pendant) {
        return - _ctl(*_insertEdge, _newLeafName, pendant, distal, d-distal);
    };

    return minimize(fn, pendant_start, 0, 2.0, max_iters);
}

}} // namespace sts::online
