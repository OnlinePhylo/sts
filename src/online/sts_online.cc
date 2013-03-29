#include <Bpp/Phyl/Io/NexusIoTree.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/UniformDiscreteDistribution.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

#include <Bpp/Phyl/Model/JCnuc.h>

// TEST CODE
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>

#include "tclap/CmdLine.h"

#include "smctc.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "config.h"
#include "beagle_tree_likelihood.h"
#include "online_add_sequence_move.h"
#include "online_smc_init.h"
#include "multiplier_mcmc_move.h"
#include "node_slider_mcmc_move.h"
#include "tree_particle.h"
#include "util.h"


namespace cl = TCLAP;
using namespace std;
typedef bpp::TreeTemplate<bpp::Node> Tree;

using namespace sts::online;

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
                    STS_VERSION);
    cl::ValueArg<int> burnin("b", "burnin-count", "Number of trees to discard as burnin", false, 0, "#", cmd);
    cl::ValueArg<int> particle_factor("p", "particle-factor", "Multiple of number of trees to determine particle count",
                                      false, 1, "#", cmd);
    cl::ValueArg<int> mcmc_count("m", "mcmc-moves", "Number of MCMC moves per-particle",
                                 false, 5, "#", cmd);
    cl::UnlabeledValueArg<string> alignment_path(
        "alignment", "Input fasta alignment.", true, "", "fasta", cmd);
    cl::UnlabeledValueArg<string> tree_posterior(
        "posterior_trees", "Posterior tree file in NEXUS format",
        true, "", "trees.nex", cmd);
    //cl::UnlabeledValueArg<string> param_posterior(
        //"posterior_params", "Posterior parameter file, tab delimited",
        //true, "", "params", cmd);

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
    // Discard burnin
    if(burnin.getValue() > 0) {
        if(burnin.getValue() >= trees.size()) {
            cerr << "Burnin (" << burnin.getValue() << ") exceeds number of trees (" << trees.size() << ")\n";
            return 1;
        }
        trees.erase(trees.begin(), trees.begin() + burnin.getValue());
    }
    cerr << "read " << trees.size() << " trees" << endl;

    ifstream alignment_fp(alignment_path.getValue());
    unique_ptr<bpp::SiteContainer> sites(sts::util::read_alignment(alignment_fp, &DNA));
    alignment_fp.close();
    bpp::VectorSiteContainer ref(&DNA), query(&DNA);
    partition_alignment(*sites, trees[0]->getLeavesNames(), ref, query);
    cerr << ref.getNumberOfSequences() << " reference sequences" << endl;
    cerr << query.getNumberOfSequences() << " query sequences" << endl;

    // TODO: allow model specification
    bpp::JCnuc model(&DNA);
    // TODO: Allow rate distribution specification
    bpp::ConstantDistribution rate_dist(1.0);
    //bpp::GammaDiscreteDistribution rate_dist(4, 0.358);

    vector<Tree_particle> particles;
    particles.reserve(trees.size());
    for(unique_ptr<Tree>& tree : trees) {
        particles.emplace_back(model.clone(), tree.release(), rate_dist.clone(), &ref);
    }

    Beagle_tree_likelihood calculator(*sites, model, rate_dist);

    // SMC
    Online_smc_init init_fn(particles);
    Online_add_sequence_move move_fn(calculator, query.getSequencesNames());

    smc::sampler<Tree_particle> sampler(particle_factor.getValue() * trees.size(), SMC_HISTORY_NONE);
    smc::mcmc_moves<Tree_particle> mcmc_moves;
    mcmc_moves.AddMove(Multiplier_mcmc_move(calculator), 4.0);
    mcmc_moves.AddMove(Node_slider_mcmc_move(calculator), 1.0);
    smc::moveset<Tree_particle> moveset(init_fn, move_fn);
    moveset.SetMCMCSelector(mcmc_moves);
    moveset.SetNumberOfMCMCMoves(mcmc_count.getValue());

    sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED, 0.99);
    sampler.SetMoveSet(moveset);
    sampler.Initialise();
    size_t n_query = query.getNumberOfSequences();
    vector<string> sequence_names =query.getSequencesNames();
    for(size_t n = 0; n < n_query; n++) {
        const double ess = sampler.IterateEss();
        cerr << "Iter " << n << ": ESS=" << ess << " sequence=" << sequence_names[n] << endl;
    }

    for(size_t i = 0; i < sampler.GetNumber(); i++) {
        const Tree_particle& p = sampler.GetParticleValue(i);
        calculator.initialize(*p.model, *p.rate_dist, *p.tree);
        double log_weight = calculator.calculate_log_likelihood();
        string s = bpp::TreeTemplateTools::treeToParenthesis(*p.tree);
        cout << log_weight << '\t' << s;
    }
}
