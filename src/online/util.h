/// \file util.h
#ifndef STS_PARTICLE_UTIL_H
#define STS_PARTICLE_UTIL_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>

namespace sts
{
namespace likelihood
{
class Online_calculator;
}
namespace util
{


bpp::SiteContainer* read_alignment(std::istream &, const bpp::Alphabet *);

bpp::SiteContainer* unique_sites(const bpp::SiteContainer& sites, bool verbose = false);

std::vector<double> compressed_site_weights(const bpp::SiteContainer&, const bpp::SiteContainer&);

/// Print a vector to an output stream, tab delimited.
/// \param vec Vector
/// \param out Destination stream
template <typename T>
void print_vector(const std::vector<T>& vec, std::ostream& out = std::cout)
{
    for(const T & i : vec) {
        out << i << '\t';
    }
    out << std::endl;
}

std::string beagle_errstring(const int beagle_error_code);
void beagle_check(int return_code);

/// Calculate the ESS of a set of log weights.
///
/// \[ ESS = \frac{1}{\sum_i^N w_i^2 \]
double effectiveSampleSize(const std::vector<double>& logWeights);
    
double logit(double x);
double logitinv(double x);


} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_UTIL_H
