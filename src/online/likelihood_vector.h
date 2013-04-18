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
    inline const double* data() const { return v.data(); }

    /// \brief Vector product with \c other, for equiprobable rates / frequencies
    ///
    /// Given this vector, \f$x\f$, and \c other \f$y\f$, computes:
    ///
    /// \f[
    ///    \sum_{i \in sites} \log \left(\sum_{j \in rates} \; \sum_{k \in states} x_{ijk} y_{ijk} \right)
    /// \f]
    ///
    double logDot(const LikelihoodVector& other) const;

    /// \brief Vector product with \c other, with given rate probabilities and state frequencies.
    ///
    /// \f[
    ///    \sum_{i \in sites} \log \left(\sum_{j \in rates} w_j \sum_{k \in states} x_{ijk} y_{ijk} \right)
    /// \f]
    ///
    /// \param other Another likelihood vector
    /// \param model Substitution model
    /// \param rateDist Rate distribution
    double logDot(const LikelihoodVector& other, const bpp::SubstitutionModel& model, const bpp::DiscreteDistribution& rateDist) const;
    double logDot(const LikelihoodVector& other, const std::vector<double>& freqs, const std::vector<double>& rateWeights) const;

    /// \brief Log-likelihood of this vector
    double logLikelihood() const;
    /// \brief Log-likelihood of this vector, using rate weights and frequencies from Bio++
    double logLikelihood(const bpp::SubstitutionModel& model, const bpp::DiscreteDistribution& rate_dist) const;
    /// \brief Log-likelihood of this vector, using rate weights and frequencies from vectors
    double logLikelihood(const std::vector<double>& freqs, const std::vector<double>& rateWeights) const;

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
