/// \file util.h
#ifndef STS_PARTICLE_UTIL_H
#define STS_PARTICLE_UTIL_H

#include "node_ptr.h"
#include "particle.h"
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

int count_uncoalesced_trees(const std::vector<particle::Node_ptr> &);
std::vector<particle::Node_ptr> uncoalesced_nodes(particle::Particle pp, std::vector<particle::Node_ptr> leaf_nodes);

void write_tree(std::ostream &out, const particle::Node_ptr root, const std::unordered_map<particle::Node_ptr, std::string> &names);
bpp::SiteContainer* read_alignment(std::istream &, const bpp::Alphabet *);

void register_nodes(likelihood::Online_calculator&,
                    const particle::Node_ptr,
                    const std::unordered_map<particle::Node_ptr, std::string>&);

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

/// Combine hashed value of \c t with \c seed.
///
/// From Boost
template<typename T>
void hash_combine(size_t& seed, const T& t)
{
    std::hash<T> h;
    seed ^= h(t) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/// Reverse a map
///
/// \pre values in m are unique
template<typename K, typename V>
std::unordered_map<V, K> unordered_invert(const std::unordered_map<K,V>& m)
{
    std::unordered_map<V, K> result;
    for(auto &item : m) {
        if(result.count(item.second))
            throw std::out_of_range("duplicate map value during invert");
        result[item.second] = item.first;
    }
    return result;
}


} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_UTIL_H
