#ifndef STS_ONLINE_LIKELIHOOD_VECTOR
#define STS_ONLINE_LIKELIHOOD_VECTOR

#include <cstddef>
#include <vector>

namespace bpp {
class DiscreteDistribution;
class SubstitutionModel;
}

namespace sts { namespace online {

/// A container for likelihood vectors
class LikelihoodVector
{
public:
    /// Constructor
    LikelihoodVector(const size_t n_rates, const size_t n_sites, const size_t n_states);

    /// Move constructor
    LikelihoodVector(LikelihoodVector&& other);

    /// Move assignment
    LikelihoodVector& operator=(LikelihoodVector&& other);

    inline const std::vector<double>& get() const { return v; };
    inline std::vector<double>& get() { return v; };

    /// \brief Access the likelihood value at <c>(rate,site,state)</c>
    inline double& operator()(const size_t rate, const size_t site, const size_t state)
    {
        return v[index(rate, site, state)];
    }

    /// \brief Access the likelihood value at <c>(rate,site,state)</c>
    inline double operator()(const size_t rate, const size_t site, const size_t state) const
    {
        return v[index(rate, site, state)];
    }

    /// Alias for calling \c .data() on #get()
    inline double* data() { return v.data(); }

    /// Alias for calling \c .data() on #get()

    inline size_t nRates() const { return nRates_; }
    inline size_t nSites() const { return nSites_; }
    inline size_t nStates() const { return nStates_; }
private:
    size_t nRates_,
           nSites_,
           nStates_;

    /// Actual storage
    std::vector<double> v;

    /// Index of <c>(rate, site, state)</c> within \c v
    size_t index(const size_t rate, const size_t site, const size_t state) const;

};

}}

#endif
