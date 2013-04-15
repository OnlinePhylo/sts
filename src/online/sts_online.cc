#include <Bpp/Phyl/Io/NexusIoTree.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/UniformDiscreteDistribution.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

#include <Bpp/Phyl/Model/JCnuc.h>

#include <gsl/gsl_randist.h>

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
#include "branch_length_prior.h"
#include "beagle_tree_likelihood.h"
#include "composite_tree_likelihood.h"
#include "online_add_sequence_move.h"
#include "online_smc_init.h"
#include "multiplier_mcmc_move.h"
#include "node_slider_mcmc_move.h"
#include "multiplier_smc_move.h"
#include "node_slider_smc_move.h"
#include "tree_particle.h"
#include "weighted_selector.h"
#include "util.h"


namespace cl = TCLAP;
using namespace std;
typedef bpp::TreeTemplate<bpp::Node> Tree;

using namespace sts::online;

const bpp::DNA DNA;

/// Partition an alignment into reference and query sequences
/// \param allSequences Site container with all sequences
/// \param taxaInTree Names of reference sequences
/// \param ref *out* Reference alignment
/// \param query *out* Query alignment.
void partitionAlignment(const bpp::SiteContainer& allSequences,
                         const vector<string> taxaInTree,
                         bpp::SiteContainer& ref,
                         bpp::SiteContainer& query)
{
    unordered_set<string> ref_taxa(begin(taxaInTree), end(taxaInTree));
    for(size_t i = 0; i < allSequences.getNumberOfSequences(); i++) {
        const bpp::Sequence& sequence = allSequences.getSequence(i);
        if(ref_taxa.count(sequence.getName()))
            ref.addSequence(sequence, false);
        else
            query.addSequence(sequence, false);
    }
}

vector<unique_ptr<Tree>> readTrees(bpp::IMultiTree& reader, std::string path)
{
    vector<bpp::Tree*> unmanagedTrees;
    reader.read(path, unmanagedTrees);
    vector<unique_ptr<Tree>> result;
    result.reserve(unmanagedTrees.size());
    for(bpp::Tree* t : unmanagedTrees) {
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

template<typename T>
class RangeConstraint : public TCLAP::Constraint<T>
{
public:
    RangeConstraint(T minValue, T maxValue, bool inclusive=true) :
        minValue(minValue),
        maxValue(maxValue),
        inclusive(inclusive)
    {};

    std::string shortID() const
    {
        const char start = inclusive ? '[' : '(',
                   end = inclusive ? ']' : ')';
        return start + std::to_string(minValue) + ',' + std::to_string(maxValue) + end;
    };
    std::string description() const { return shortID(); };

    bool check(const T& val) const
    {
        if(inclusive)
            return val >= minValue && val <= maxValue;
        else
            return val > minValue && val < maxValue;
    }
private:
    T minValue, maxValue;
    bool inclusive;
};

int main(int argc, char **argv)
{
    cl::CmdLine cmd("Run STS starting from an extant posterior", ' ',
                    sts::STS_VERSION);
    cl::ValueArg<int> burnin("b", "burnin-count", "Number of trees to discard as burnin", false, 0, "#", cmd);

    RangeConstraint<double> resample_range(0.0, 1.0, false);
    cl::ValueArg<double> resample_threshold("", "resample-threshold", "Resample when the ESS falls below T * n_particles",
                                            false, 0.99, &resample_range, cmd);
    cl::ValueArg<int> particleFactor("p", "particle-factor", "Multiple of number of trees to determine particle count",
                                      false, 1, "#", cmd);
    cl::ValueArg<int> mcmcCount("m", "mcmc-moves", "Number of MCMC moves per-particle",
                                 false, 0, "#", cmd);
    cl::ValueArg<int> treeSmcCount("", "tree-moves", "Number of additional tree-altering SMC moves per added sequence",
                                 false, 0, "#", cmd);
    cl::ValueArg<double> blPriorExpMean("", "edge-prior-exp-mean", "Mean of exponential prior on edges",
                                           false, 0.1, "float", cmd);

    cl::UnlabeledValueArg<string> alignmentPath(
        "alignment", "Input fasta alignment.", true, "", "fasta", cmd);
    cl::UnlabeledValueArg<string> treePosterior(
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
    bpp::NexusIOTree treeReader;
    vector<unique_ptr<Tree>> trees;
    try {
        trees = readTrees(treeReader, treePosterior.getValue());
    } catch(bpp::Exception &e) {
        cerr << "error reading " << treePosterior.getValue() << ": " << e.what() << endl;
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

    ifstream alignment_fp(alignmentPath.getValue());
    unique_ptr<bpp::SiteContainer> sites(sts::util::read_alignment(alignment_fp, &DNA));
    alignment_fp.close();
    bpp::VectorSiteContainer ref(&DNA), query(&DNA);
    partitionAlignment(*sites, trees[0]->getLeavesNames(), ref, query);
    cerr << ref.getNumberOfSequences() << " reference sequences" << endl;
    cerr << query.getNumberOfSequences() << " query sequences" << endl;

    // TODO: allow model specification
    bpp::JCnuc model(&DNA);
    // TODO: Allow rate distribution specification
    bpp::ConstantDistribution rate_dist(1.0);
    //bpp::GammaDiscreteDistribution rate_dist(4, 0.358);

    // TODO: Other distributions
    const double expPriorMean = blPriorExpMean.getValue();
    std::function<double(double)> exponentialPrior = [expPriorMean](const double d) {
        return std::log(gsl_ran_exponential_pdf(d, expPriorMean));
    };

    vector<TreeParticle> particles;
    particles.reserve(trees.size());
    for(unique_ptr<Tree>& tree : trees) {
        particles.emplace_back(model.clone(), tree.release(), rate_dist.clone(), &ref);
    }

    std::shared_ptr<BeagleTreeLikelihood> beagleLike =
        make_shared<BeagleTreeLikelihood>(*sites, model, rate_dist);
    CompositeTreeLikelihood tree_like(beagleLike);
    tree_like.add(BranchLengthPrior(exponentialPrior));

    const int treeMoveCount = treeSmcCount.getValue();
    // move selection
    std::vector<smc::moveset<TreeParticle>::move_fn> smcMoves;
    smcMoves.push_back(OnlineAddSequenceMove(tree_like, query.getSequencesNames()));

    WeightedSelector<size_t> additionalSMCMoves;
    if(treeMoveCount) {
        smcMoves.push_back(MultiplierSMCMove(tree_like));
        smcMoves.push_back(NodeSliderSMCMove(tree_like));

        // Twice as many multipliers
        additionalSMCMoves.push_back(1, 20);
        additionalSMCMoves.push_back(2, 5);
    }

    std::function<long(long, const smc::particle<TreeParticle>&, smc::rng*)> moveSelector =
        [treeMoveCount,&query,&smcMoves,&additionalSMCMoves](long time, const smc::particle<TreeParticle>&, smc::rng* rng) -> long {
       const size_t blockSize = 1 + treeMoveCount;

       // Add a sequence, followed by treeMoveCount randomly selected moves
       const bool addSequenceStep = (time - 1) % blockSize == 0;
       if(addSequenceStep)
           return 0;
       return additionalSMCMoves.choice();
    };

    // SMC
    OnlineSMCInit particleInitializer(particles);

    smc::sampler<TreeParticle> sampler(particleFactor.getValue() * trees.size(), SMC_HISTORY_NONE);
    smc::mcmc_moves<TreeParticle> mcmcMoves;
    mcmcMoves.AddMove(MultiplierMCMCMove(tree_like), 4.0);
    mcmcMoves.AddMove(NodeSliderMCMCMove(tree_like), 1.0);
    smc::moveset<TreeParticle> moveSet(particleInitializer, moveSelector, smcMoves, mcmcMoves);
    moveSet.SetNumberOfMCMCMoves(mcmcCount.getValue());

    sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED, resample_threshold.getValue());
    sampler.SetMoveSet(moveSet);
    sampler.Initialise();
    const size_t nIters = (1 + treeMoveCount) * query.getNumberOfSequences();
    vector<string> sequenceNames = query.getSequencesNames();
    for(size_t n = 0; n < nIters; n++) {
        const double ess = sampler.IterateEss();
        cerr << "Iter " << n << ": ESS=" << ess << " sequence=" << sequenceNames[n / (1 + treeMoveCount)] << endl;
    }

    double maxWeight = -1e10;
    for(size_t i = 0; i < sampler.GetNumber(); i++) {
        const TreeParticle& p = sampler.GetParticleValue(i);
        beagleLike->initialize(*p.model, *p.rateDist, *p.tree);
        const double logWeight = beagleLike->calculateLogLikelihood();
        maxWeight = std::max(logWeight, maxWeight);
        string s = bpp::TreeTemplateTools::treeToParenthesis(*p.tree);
        cout << logWeight << '\t' << s;
    }

    clog << "Maximum LL: " << maxWeight << '\n';
}
