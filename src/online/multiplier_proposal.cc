#include "multiplier_proposal.h"

namespace sts { namespace online {

/// \brief Scale \c value
///
/// Multiplies \c value by multiplier \f$y\f$, using tuning parameter \f$\lambda\f$, where
/// \f{eqnarray*}{
/// y &=& \exp^{\lambda(u - 0.5)}\\
/// u & \sim & Uniform(0, 1)
/// \f}
///
/// The proposal density can be derived from the Jacobian:
///
/// \f{eqnarray*}{
/// x^\ast &=& x e^{\lambda(u - 0.5)}\\
/// u &=& \frac{\log(x^\ast/x)}{\lambda} + 0.5 \\
/// \frac{d u}{d x^\ast} &=& \frac{1}{\lambda x^\ast}\\
/// & & \\
/// q(x^\ast|x) &=& p(u(x^\ast)) \left|\frac{d u}{d x^\ast}\right|\\
/// &=& \frac{1}{\lambda x^\ast}\\
/// q(x|x^\ast)/q(x^\ast|x) &=& x^\ast/x
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

    const double forwardDensity = 1 / (tuning * new_value);

    return Proposal{new_value, forwardDensity, new_value / value};
}

}}
