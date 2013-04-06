#include "multiplier_proposal.h"

namespace sts { namespace online {

/// \brief Scale \c value
///
/// Multiplies \c value by multiplier \f$y\f$, using tuning parameter \f$\lambda\f$, where
/// \f{eqnarray*}{
/// y &=& \exp^{\lambda(x - 0.5)}\\
/// x & \tilde{} & Uniform(0, 1)
/// \f}
///
/// The proposal density can be derived from the Jacobian:
///
/// \f{eqnarray*}{
/// y &=& e^{\lambda(x - 0.5)}\\
/// x &=& \frac{\log(y)}{\lambda} - 0.5 \lambda \\
/// \frac{d x}{d y} &=& \frac{1}{\lambda y}\\
/// & & \\
/// p(y) &=& p(x) \left|\frac{d x}{d y}\right|\\
/// p(y) &=& \frac{1}{\lambda y}
/// \f}
Proposal positive_real_multiplier(const double value, const double min_value, const double max_value, const double tuning, smc::rng* rng)
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

    const double forwardDensity = 1 / (tuning * factor);

    return Proposal{new_value, forwardDensity, new_value / value};
}

}}
