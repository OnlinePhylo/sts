#include "multiplier_proposal.h"

namespace sts { namespace online {

Proposal pos_real_multiplier(const double value, const double min_value, const double max_value, const double tuning, smc::rng* rng)
{
    const double random = rng->UniformS();
    const double factor = std::exp(tuning * (random - 0.5));
    double new_value = factor * value;
    bool valid = false;
    do {
        if(new_value < min_value)
            new_value = min_value * min_value / new_value;
        else if(new_value > max_value)
            new_value = max_value * max_value / new_value;
        else
            valid = true;
    } while (!valid);

    return {new_value, new_value / value};
}

}}
