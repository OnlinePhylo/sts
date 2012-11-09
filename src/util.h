#ifndef STS_PARTICLE_UTIL_HPP
#define STS_PARTICLE_UTIL_HPP

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

int uncoalesced_count_trees(const std::vector<particle::Node_ptr> &);
std::vector<particle::Node_ptr> uncoalesced_nodes(particle::Particle pp, std::vector<particle::Node_ptr> leaf_nodes);

void write_tree(std::ostream &out, const particle::Node_ptr root, const std::unordered_map<particle::Node_ptr, std::string> &names);
bpp::SiteContainer* read_alignment(std::istream &, const bpp::Alphabet *);

void register_nodes(likelihood::Online_calculator&,
                    const particle::Node_ptr,
                    const std::unordered_map<particle::Node_ptr, std::string>&);

bpp::SiteContainer* unique_sites(const bpp::SiteContainer& sites, bool verbose = false);
std::vector<double> compressed_site_weights(const bpp::SiteContainer&, const bpp::SiteContainer&);

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_UTIL_HPP
