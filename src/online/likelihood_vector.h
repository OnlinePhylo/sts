#ifndef STS_ONLINE_LIKELIHOOD_VECTOR
#define STS_ONLINE_LIKELIHOOD_VECTOR

#include <vector>
#include <cstddef>

namespace sts { namespace online {

/// A container for likelihood vectors
class Likelihood_vector
{
public:
    /// Constructor
    Likelihood_vector(const size_t n_rates, const size_t n_sites, const size_t n_states);

    /// Move constructor
    Likelihood_vector(Likelihood_vector&& other);

    /// Move assignment
    Likelihood_vector& operator=(Likelihood_vector&& other);

    inline const std::vector<double>& get() const { return v; };
    inline std::vector<double>& get() { return v; };

    /// \brief Access a likelihov value for a <c>(rate,site,state)</c> combination
    inline double& operator()(const size_t rate, const size_t site, const size_t state)
    {
        return v[index(rate, site, state)];
    }

    /// \brief Access a likelihov value for a <c>(rate,site,state)</c> combination
    inline double operator()(const size_t rate, const size_t site, const size_t state) const
    {
        return v[index(rate, site, state)];
    }

    /// \brief Vector product with \c other.
    double log_dot(const Likelihood_vector& other) const;

    inline size_t n_rates() const { return n_rates_; }
    inline size_t n_sites() const { return n_sites_; }
    inline size_t n_states() const { return n_states_; }
private:
    size_t n_rates_,
           n_sites_,
           n_states_;

    /// Actual storage
    std::vector<double> v;

    /// Index of <c>(rate, site, state)</c> within \c v
    size_t index(const size_t rate, const size_t site, const size_t state) const;
};

}}

#endif
