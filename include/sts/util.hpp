#ifndef STS_PARTICLE_UTIL_HPP
#define STS_PARTICLE_UTIL_HPP

#include <cassert>
#include <iostream>
#include <memory>
#include <numeric>
#include <stack>
#include <string>
#include <unordered_set>
#include <vector>

#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Io/IoSequenceFactory.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/SiteTools.h>

#include "sts/likelihood/online_calculator.hpp"
#include "sts/particle/node.hpp"
#include "sts/particle/state.hpp"

namespace sts
{
namespace util
{

int uncoalesced_count_trees(const std::vector<particle::Node_ptr> &);
std::vector<particle::Node_ptr> uncoalesced_nodes(particle::Particle pp, std::vector<particle::Node_ptr> leaf_nodes);

void write_tree(std::ostream &out, const particle::Node_ptr root, const std::unordered_map<particle::Node_ptr, std::string> &names);
bpp::SiteContainer* read_alignment(std::istream &, const bpp::Alphabet *);

/// Find the number of trees (that is, trees consisting of more than one node) from a collection of uncoalesced nodes.
/// \param uncoalesced The uncoalesced nodes.
/// \return The count.
int uncoalesced_count_trees(const std::vector<particle::Node_ptr> &uncoalesced)
{
    int result = 0;
    for(auto i : uncoalesced) {
        if(!i->is_leaf())
            result++;
    }
    return result;
}

/// Find the uncoalesced nodes for a particle.
/// \param pp Input particle
/// \return vector of uncoalesced phylo_nodes.
std::vector<particle::Node_ptr> uncoalesced_nodes(const particle::Particle pp, const std::vector<particle::Node_ptr> leaf_nodes)
{
    // Our set of phylo nodes that can be used in proposal.
    std::unordered_set<particle::Node_ptr> proposal_set;
    // The nodes that have already been coalesced, to be removed later.
    std::unordered_set<particle::Node_ptr> coalesced;
    // Insert all of the leaf nodes into the proposal set.
    proposal_set.insert(leaf_nodes.begin(), leaf_nodes.end());
    // Walk back to predecessor particles, adding root nodes to
    // proposal_set and collecting coalesced nodes in `coalesced`.
    for(particle::Particle cur = pp->predecessor; cur != NULL; cur = cur->predecessor) {
        // Skip if the particle is \perp.
        if(cur->node == NULL) continue;
        // Skip if we've already processed this subtree, such that it's already found in coalesced.
        if(coalesced.find(cur->node) != coalesced.end()) continue;
        // Insert this active root node to the proposal set.
        proposal_set.insert(cur->node);
        // Recursively add all descendants of the root nodes to the coalesced set using a std::stack.
        std::stack<particle::Node_ptr> s;
        s.push(cur->node);
        while(s.size() > 0) {
            particle::Node_ptr n = s.top();
            s.pop();
            if(n->is_leaf()) continue;	// leaf node, nothing more to do.
            coalesced.insert(n->child1->node);
            coalesced.insert(n->child2->node);
            s.push(n->child1->node);
            s.push(n->child2->node);
        }
    }

    std::vector<particle::Node_ptr> pvec(proposal_set.begin(), proposal_set.end());
    std::vector<particle::Node_ptr> cvec(coalesced.begin(), coalesced.end());
    sort(pvec.begin(), pvec.end());
    sort(cvec.begin(), cvec.end());

    // The set difference of available (i.e. proposal_set) and coalesced nodes yields the final proposal set; store it
    // in prop_vector.
    std::vector<particle::Node_ptr> prop_vector(proposal_set.size() + coalesced.size());
    // UGH: std::set_difference requires an ordered container class
    // AG: that's the only way to do a set difference efficiently, right?
    auto last_ins = set_difference(pvec.begin(), pvec.end(), cvec.begin(),
                                   cvec.end(), prop_vector.begin());
    prop_vector.resize(last_ins - prop_vector.begin());

    return prop_vector;
}

void write_tree(std::ostream &out, const particle::Node_ptr root, const std::unordered_map<particle::Node_ptr, std::string>& names)
{
    std::unordered_set<particle::Node_ptr> visited;
    std::stack<particle::Node_ptr> s;
    s.push(root);
    while(!s.empty()) {
        particle::Node_ptr cur = s.top();
        if(cur->is_leaf()) {
            auto iter = names.find(cur);
            out << iter->second;
            visited.insert(cur);
            s.pop();
            continue;
        }
        if(!visited.count(cur->child1->node)) {
            out << "(";
            s.push(cur->child1->node);
            continue;
        } else if(!visited.count(cur->child2->node)) {
            out << ":" << cur->child1->length << ",";
            s.push(cur->child2->node);
            continue;
        }
        out << ":" << cur->child2->length << ")";
        visited.insert(cur);
        s.pop();
    }
    out << ";\n";
}

/// Read an alignment from a stream
/// \param in Input stream
/// \param alphabet The alphabet to use.
bpp::SiteContainer* read_alignment(std::istream &in, const bpp::Alphabet *alphabet)
{
    // Holy boilerplate - Bio++ won't allow reading FASTA files as alignments
    bpp::IOSequenceFactory fac;
    std::unique_ptr<bpp::ISequence> reader = std::unique_ptr<bpp::ISequence>(
                fac.createReader(bpp::IOSequenceFactory::FASTA_FORMAT));
    std::unique_ptr<bpp::SequenceContainer> seqs = std::unique_ptr<bpp::SequenceContainer>(reader->read(in, alphabet));

    // Have to look up by name
    std::vector<std::string> names = seqs->getSequencesNames();
    bpp::SiteContainer *sequences = new bpp::VectorSiteContainer(alphabet);

    for(auto name : names) {
        sequences->addSequence(seqs->getSequence(name), true);
    }

    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sequences);

    return sequences;
}

/// Determine new site weights after compression.
/// \param orig Original sites
/// \param compressed sites after compression to unique sites
/// \returns A vector, where the value at each position is the appropriate weight for <c>compressed</c> site \c i
std::vector<double> compressed_site_weights(const bpp::SiteContainer& orig, const bpp::SiteContainer& compressed)
{
    assert(compressed.getNumberOfSites() <= orig.getNumberOfSites());

    std::vector<double> result(compressed.getNumberOfSites(), 0);

    // Get the first index at which each site from orig appears in compressed.
    std::vector<int> m = bpp::PatternTools::getIndexes(orig, compressed);
    for(unsigned int i = 0; i < m.size(); ++i) {
        assert(m[i] >= 0);
        ++result[m[i]];
    }

    unsigned int tot = std::accumulate(result.begin(), result.end(), 0);
    assert(tot == orig.getNumberOfSites());

    return result;
}

/// Get the unique sites in an alignment
/// \param sites Original sites
/// \param verbose Print a message if sites are compressed?
/// \returns A sequence container containing only the unique sites from \param sites.
bpp::SiteContainer* unique_sites(const bpp::SiteContainer& sites, bool verbose = false)
{
    bpp::SiteContainer *compressed = bpp::PatternTools::shrinkSiteSet(sites);

    if(verbose && compressed->getNumberOfSites() < sites.getNumberOfSites())
        std::cerr << "Reduced from "
                  << sites.getNumberOfSites()
                  << " to " << compressed->getNumberOfSites()
                  << " sites"
                  << std::endl;

    return compressed;
}

/// Register a tree of nodes with an Online_calculator.
///  \param calc Online_calculator instance
///  \param root Root node to register. Children are navigated.
///  \param names Map from node to taxon names.
void register_nodes(likelihood::Online_calculator& calc,
                    const particle::Node_ptr root,
                    const std::unordered_map<particle::Node_ptr, std::string>& names)
{
    std::stack<particle::Node_ptr> to_register;
    to_register.push(root);

    while(!to_register.empty()) {
        particle::Node_ptr n = to_register.top();
        to_register.pop();
        if(n->is_leaf()) {
            assert(names.count(n) == 1);
            calc.register_leaf(n, names.at(n));
        } else {
            calc.register_node(n);
            to_register.push(n->child1->node);
            to_register.push(n->child2->node);
        }
    }
}

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_UTIL_HPP
