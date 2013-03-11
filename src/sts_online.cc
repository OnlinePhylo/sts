#include <Bpp/Phyl/Io/NexusIoTree.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/UniformDiscreteDistribution.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/TreeTools.h>

#include <Bpp/Phyl/Model/JCnuc.h>

// TEST CODE
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>

#include "tclap/CmdLine.h"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "beagle_tree_likelihood.h"
#include "tree_particle.h"
#include "util.h"

#define _STRINGIFY(s) #s
#define STRINGIFY(s) _STRINGIFY(s)

namespace cl = TCLAP;
using namespace std;
typedef bpp::TreeTemplate<bpp::Node> Tree;

const bpp::DNA DNA;

/// Partition an alignment into reference and query sequences
/// \param all_sequences Site container with all sequences
/// \param taxa_in_tree Names of reference sequences
/// \param ref *out* Reference alignment
/// \param query *out* Query alignment.
void partition_alignment(const bpp::SiteContainer& all_sequences,
                         const vector<string> taxa_in_tree,
                         bpp::SiteContainer& ref,
                         bpp::SiteContainer& query)
{
    unordered_set<string> ref_taxa(begin(taxa_in_tree), end(taxa_in_tree));
    for(size_t i = 0; i < all_sequences.getNumberOfSequences(); i++) {
        const bpp::Sequence& sequence = all_sequences.getSequence(i);
        if(ref_taxa.count(sequence.getName()))
            ref.addSequence(sequence, false);
        else
            query.addSequence(sequence, false);
    }
}

vector<unique_ptr<Tree>> read_trees(bpp::IMultiTree& reader, std::string path)
{
    vector<bpp::Tree*> unmanaged_trees;
    reader.read(path, unmanaged_trees);
    vector<unique_ptr<Tree>> result;
    result.reserve(unmanaged_trees.size());
    for(bpp::Tree* t : unmanaged_trees) {
        Tree *tt = new Tree(*t);
        delete(t);

        // Trees must be bifurcating for use with BEAGLE.
        // Root by making the first leaf an outgroup
        tt->newOutGroup(tt->getLeaves()[0]);
        tt->resetNodesId();
        assert(!tt->isMultifurcating());
        assert(tt->isRooted());
        result.emplace_back(tt);
    }
    return result;
}

int main(int argc, char **argv)
{
    cl::CmdLine cmd("Run STS starting from an extant posterior", ' ',
                   STRINGIFY(STS_VERSION));
    cl::UnlabeledValueArg<string> alignment_path(
        "alignment", "Input fasta alignment.", true, "", "fasta", cmd);
    cl::UnlabeledValueArg<string> tree_posterior(
        "posterior_trees", "Posterior tree file in NEXUS format",
        true, "", "trees.nex", cmd);
    cl::UnlabeledValueArg<string> param_posterior(
        "posterior_params", "Posterior parameter file, tab delimited",
        true, "", "params", cmd);

    try {
        cmd.parse(argc, argv);
    } catch(TCLAP::ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    // Get alignment

    // Read trees
    bpp::NexusIOTree tree_reader;
    vector<unique_ptr<Tree>> trees;
    try {
        trees = read_trees(tree_reader, tree_posterior.getValue());
    } catch(bpp::Exception &e) {
        cerr << "error reading " << tree_posterior.getValue() << ": " << e.what() << endl;
        return 1;
    }
    cout << "read " << trees.size() << " trees" << endl;

    ifstream alignment_fp(alignment_path.getValue());
    unique_ptr<bpp::SiteContainer> sites(sts::util::read_alignment(alignment_fp, &DNA));
    alignment_fp.close();
    bpp::VectorSiteContainer ref(&DNA), query(&DNA);
    partition_alignment(*sites, trees[0]->getLeavesNames(), ref, query);
    cout << ref.getNumberOfSequences() << " reference sequences" << endl;
    cout << query.getNumberOfSequences() << " query sequences" << endl;

    // TODO: allow model specification
    bpp::JCnuc model(&DNA);
    // TODO: Allow rate distribution specification
    //bpp::ConstantDistribution rate_dist(1.0);
    bpp::GammaDiscreteDistribution rate_dist(4, 0.358);

    vector<sts::Tree_particle> particles;
    size_t i = 0;
    for(unique_ptr<Tree>& tree : trees) {
        bpp::DRHomogeneousTreeLikelihood like(*tree, &model, &rate_dist, false, false);
        like.setData(ref);
        like.initialize();
        const double bpp_ll = -like.getValue();

        sts::likelihood::Beagle_tree_likelihood btl(ref, model, rate_dist);
        btl.load_substitution_model(model);
        btl.load_rate_distribution(rate_dist);
        const double beagle_ll = btl.calculate_log_likelihood(*tree);

        //std::vector<double> site_log_likes(ref.getNumberOfSites(), -900000);
        //beagleGetSiteLogLikelihoods(btl.get_beagle_instance(), site_log_likes.data());
        //for(size_t i = 0; i < 15; i++) std::cout << site_log_likes[i] << '\t';
        //std::cout << '\n';
        //site_log_likes = like.getLogLikelihoodForEachSite();
        //for(size_t i = 0; i < 15; i++) std::cout << site_log_likes[i] << '\t';
        //std::cout << '\n';

        //std::cout << bpp::TreeTemplateTools::treeToParenthesis(*tree);
        std::cout << "tree likelihood[" << i++ << "] BEAGLE = " << beagle_ll;
        std::cout << "\tBio++  = " << bpp_ll << "\tdelta=" << std::abs(beagle_ll - bpp_ll) << std::endl;

        //particles.emplace_back(model.clone(), tree.release(), rate_dist.clone(), &ref);
    }
}
