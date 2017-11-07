#include <Bpp/Phyl/Io/NexusIoTree.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>

#include "tclap/CmdLine.h"

#include "smctc.hh"
#include "json/json.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "sts_config.h"
#include "branch_length_prior.h"
#include "simple_flexible_tree_likelihood.h"
#ifndef NO_BEAGLE
#include "beagle_flexible_tree_likelihood.h"
#endif
#include "flexible_tree_likelihood.h"
#include "composite_tree_likelihood.h"
#include "uniform_online_add_sequence_move.h"
#include "uniform_length_online_add_sequence_move.h"
#include "gsl.h"
#include "guided_online_add_sequence_move.h"
#include "lcfit_online_add_sequence_move.h"
#include "online_smc_init.h"
#include "lazarus_mcmc_move.h"
#include "multiplier_mcmc_move.h"
#include "node_slider_mcmc_move.h"
#include "multiplier_smc_move.h"
#include "node_slider_smc_move.h"
#include "tree_particle.h"
#include "weighted_selector.h"
#include "util.h"

#include "proposal_guided_parsimony.h"
#include "flexible_parsimony.h"
#include "node_sliding_window_mcmc_move.h"

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
/// \param leaf_ids *out* IDs for the leaves in order of addition
void partitionAlignment(const bpp::SiteContainer& allSequences,
                        const vector<string> taxaInTree,
                        bpp::SiteContainer& ref,
                        bpp::SiteContainer& query,
                        unordered_map<string, size_t>& leaf_ids)
{
    unordered_set<string> ref_taxa(begin(taxaInTree), end(taxaInTree));
    size_t ref_i=0;
    size_t query_i=taxaInTree.size();
    for(size_t i = 0; i < allSequences.getNumberOfSequences(); i++) {
        const bpp::Sequence& sequence = allSequences.getSequence(i);
        if(ref_taxa.count(sequence.getName())){
            ref.addSequence(sequence, false);
            leaf_ids[sequence.getName()]=ref_i++;
        }else{
            query.addSequence(sequence, false);
            leaf_ids[sequence.getName()]=query_i++;
        }
    }
}

vector<unique_ptr<Tree>> readTrees(bpp::IMultiTree& reader, std::string path)
{
    vector<bpp::Tree*> unmanagedTrees;
    reader.read(path, unmanagedTrees);
    vector<unique_ptr<Tree>> result;
    result.reserve(unmanagedTrees.size());
    string leaf_root = "";
    for(bpp::Tree* t : unmanagedTrees) {
        Tree *tt = new Tree(*t);
        delete(t);
        if(result.size() == 0){
            leaf_root = tt->getLeaves()[0]->getName();
        }

        // Trees must be bifurcating for use with BEAGLE.
        // Every tree is rooted using the same leaf as outgroup (first leaf of first tree)
        bpp::Node* leaf = tt->getNode(leaf_root);
        tt->newOutGroup(leaf);
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
    cl::ValueArg<int> treeSmcCount("", "tree-moves",
                                   "Number of additional tree-altering SMC moves per added sequence",
                                   false, 0, "#", cmd);
    cl::ValueArg<string> particleGraphPath("g", "particle-graph",
                                           "Path to write particle graph in graphviz format",
                                           false, "", "path", cmd);
    cl::ValueArg<double> blPriorExpMean("", "edge-prior-exp-mean", "Mean of exponential prior on edges",
                                           false, 0.1, "float", cmd);
    std::vector<std::string> methodNames { "uniform-edge", "uniform-length", "guided", "lcfit", "guided-parsimony" };
    cl::ValuesConstraint<std::string> allowedProposalMethods(methodNames);
    cl::ValueArg<std::string> proposalMethod("", "proposal-method", "Proposal mechanism to use", false,
                                             "lcfit", &allowedProposalMethods, cmd);
    cl::ValueArg<double> maxLength("", "max-length", "When discretizing the tree for guided moves, "
                                   "divide edges into lengths no greater than <length>",
                                   false, std::numeric_limits<double>::max(), "length", cmd);
    cl::ValueArg<size_t> subdivideTop("", "divide-top", "Subdivide the top <N> edges to bits of no longer than max-length.",
                                   false, 0, "N", cmd);
    cl::SwitchArg fribbleResampling("", "fribble", "Use fribblebits resampling method", cmd, false);
    cl::MultiArg<double> pendantBranchLengths("", "pendant-bl", "Guided move: attempt attachment with pendant bl X", false, "X", cmd);

    cl::UnlabeledValueArg<string> alignmentPath(
        "alignment", "Input fasta alignment.", true, "", "fasta", cmd);
    cl::UnlabeledValueArg<string> treePosterior(
        "posterior_trees", "Posterior tree file in NEXUS format",
        true, "", "trees.nex", cmd);

    cl::UnlabeledValueArg<string> jsonOutputPath("json_path", "JSON output path", true, "", "path", cmd);
    
    cl::ValueArg<long> seedCmd("s", "seed", "Seed for random number generator", false, -1, "#", cmd);
    
    cl::ValueArg<double> exponent("e", "exponent", "Exponent for calculting propbablity attachments in step 1", false, 0.05, "#", cmd);

    //cl::UnlabeledValueArg<string> param_posterior(
        //"posterior_params", "Posterior parameter file, tab delimited",
        //true, "", "params", cmd);

    try {
        cmd.parse(argc, argv);
    } catch(TCLAP::ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    // Register a GSL error handler that throws exceptions instead of aborting.
    gsl_set_error_handler(&sts_gsl_error_handler);
    
    long seed;
    if(seedCmd.getValue() >= 0 ){
        seed = seedCmd.getValue();
    }
    else{
        seed = time(NULL);
    }
    cout << "Seed: " << seed <<endl;
    
    gsl_rng *rng = gsl_rng_alloc (gsl_rng_rand48);
    gsl_rng_set (rng, seed);
    
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
    clog << "read " << trees.size() << " trees" << endl;
    
    vector<double> bls = trees[0]->getBranchLengths();
    double sum = std::accumulate(bls.begin(), bls.end(), 0.0);
    double mean = sum / bls.size();
    sort(bls.begin(),bls.end());
    double median = bls[bls.size()/2];
    clog << "Mean branch length: " << mean <<endl;
    clog << "Median branch length: " << median <<endl;

    ifstream alignment_fp(alignmentPath.getValue());
    unique_ptr<bpp::SiteContainer> sites(sts::util::read_alignment(alignment_fp, &DNA));
    alignment_fp.close();
    bpp::VectorSiteContainer ref(&DNA), query(&DNA);
    unordered_map<string, size_t> leaf_ids;
    partitionAlignment(*sites, trees[0]->getLeavesNames(), ref, query,leaf_ids);
    cerr << ref.getNumberOfSequences() << " reference sequences" << endl;
    cerr << query.getNumberOfSequences() << " query sequences" << endl;

    if(query.getNumberOfSequences() == 0)
        throw std::runtime_error("No query sequences!");
    
    
    // Leaf IDs reflect their ranking in the alignment (starting from 0)
    // Internal nodes ID are greater or equal than the total number of sequences
    std::vector<string> names = sites->getSequencesNames();
    for(unique_ptr<Tree>& tree : trees){
        size_t nameCounter = names.size();
        for(bpp::Node* node : tree->getNodes()){
            if(node->isLeaf()){
                size_t pos = find(names.begin(), names.end(), node->getName()) - names.begin();
                node->setId(static_cast<int>(pos));
            }
            else {
                node->setId(static_cast<int>(nameCounter));
                nameCounter++;
            }
        }
    }
    

    // TODO: allow model specification
    bpp::JCnuc model(&DNA);
    // TODO: Allow rate distribution specification
    bpp::ConstantRateDistribution rate_dist;
    //bpp::GammaDiscreteDistribution rate_dist(4, 0.358);

    // TODO: Other distributions
    const double expPriorMean = blPriorExpMean.getValue();
    std::function<double(double)> exponentialPrior = [expPriorMean](const double d) {
        return std::log(gsl_ran_exponential_pdf(d, expPriorMean));
    };

    vector<TreeParticle> particles;
    particles.reserve(trees.size());
    for(unique_ptr<Tree>& tree : trees) {
        tree->getRootNode()->getSon(0)->setDistanceToFather(tree->getRootNode()->getSon(0)->getDistanceToFather() +
                                                           tree->getRootNode()->getSon(1)->getDistanceToFather());
        tree->getRootNode()->getSon(1)->setDistanceToFather(0.0);
        
        particles.emplace_back(std::unique_ptr<bpp::SubstitutionModel>(model.clone()),
                               std::unique_ptr<bpp::TreeTemplate<bpp::Node>>(tree.release()),
                               std::unique_ptr<bpp::DiscreteDistribution>(rate_dist.clone()),
                               &ref);
    }
    for(int i = 0; i < particles.size(); i++){
        particles[i].particleID = i;
    }

    std::unique_ptr<bpp::SitePatterns> _patterns(new bpp::SitePatterns(sites.get()));
    
#ifndef NO_BEAGLE
    shared_ptr<FlexibleTreeLikelihood> beagleLike(new BeagleFlexibleTreeLikelihood(*_patterns.get(), model, rate_dist));
#else
    shared_ptr<FlexibleTreeLikelihood> beagleLike(new SimpleFlexibleTreeLikelihood(*_patterns.get(), model, rate_dist));
#endif
    
    CompositeTreeLikelihood treeLike(beagleLike);
    treeLike.add(BranchLengthPrior(exponentialPrior));

    const int treeMoveCount = treeSmcCount.getValue();
    // move selection
    std::vector<smc::moveset<TreeParticle>::move_fn> smcMoves;

    std::vector<double> pbl = pendantBranchLengths.getValue();
    if(pbl.empty())
        pbl = {0.0, median};
    
    std::unique_ptr<OnlineAddSequenceMove> onlineAddSequenceMove;
    const string& name = proposalMethod.getValue();
    if(name == "uniform-length" || name == "uniform-edge") {
        auto branchLengthProposer = [expPriorMean](smc::rng* rng) -> std::pair<double, double> {
            const double v = rng->Exponential(expPriorMean);
            const double logDensity = std::log(gsl_ran_exponential_pdf(v, expPriorMean));
            return {v, logDensity};
        };
        if(name == "uniform-length") {
            onlineAddSequenceMove.reset(new UniformLengthOnlineAddSequenceMove(treeLike, sites->getSequencesNames(), query.getSequencesNames(), branchLengthProposer));
        } else {
            onlineAddSequenceMove.reset(new UniformOnlineAddSequenceMove(treeLike, sites->getSequencesNames(), query.getSequencesNames(), branchLengthProposer));
        }
    } else{
        GuidedOnlineAddSequenceMove* p = nullptr;
        
        if(name == "guided") {
            p = new GuidedOnlineAddSequenceMove(treeLike, sites->getSequencesNames(), query.getSequencesNames(), pbl, maxLength.getValue(), subdivideTop.getValue());
        } else if(name == "lcfit") {
            p = new LcfitOnlineAddSequenceMove(treeLike, sites->getSequencesNames(), query.getSequencesNames(), pbl, maxLength.getValue(), subdivideTop.getValue(), expPriorMean);
        } else if(name == "guided-parsimony") {
            std::shared_ptr<FlexibleParsimony> pars = make_shared<FlexibleParsimony>(*_patterns.get(), DNA);
            p = new ProposalGuidedParsimony(pars, treeLike, sites->getSequencesNames(), query.getSequencesNames(), expPriorMean);
        }
        else{
            throw std::runtime_error("Unknown sequence addition method: " + name);
        }
        p->_heating = exponent.getValue();
        onlineAddSequenceMove.reset(p);
    }


    {
        using namespace std::placeholders;
        auto wrapper = std::bind(&OnlineAddSequenceMove::operator(), std::ref(*onlineAddSequenceMove), _1, _2, _3);
        smcMoves.push_back(wrapper);
    }

    if(treeMoveCount) {
        smcMoves.push_back(MultiplierSMCMove(treeLike));
        smcMoves.push_back(NodeSliderSMCMove(treeLike));
    }

    std::function<long(long, const smc::particle<TreeParticle>&, smc::rng*)> moveSelector =
        [treeMoveCount,&query,&smcMoves](long time, const smc::particle<TreeParticle>&, smc::rng* rng) -> long {
       const size_t blockSize = 1 + treeMoveCount;

       // Add a sequence, followed by treeMoveCount randomly selected moves
       const bool addSequenceStep = (time - 1) % blockSize == 0;
       if(addSequenceStep)
           return 0;
            WeightedSelector<size_t> additionalSMCMoves{*rng};
       // Twice as many multipliers
       additionalSMCMoves.push_back(1, 20);
       additionalSMCMoves.push_back(2, 5);
       return additionalSMCMoves.choice();
    };

    // SMC
    OnlineSMCInit particleInitializer(particles);

    smc::sampler<TreeParticle> sampler(particleFactor.getValue() * trees.size(), SMC_HISTORY_RAM, gsl_rng_default, seed);
    smc::mcmc_moves<TreeParticle> mcmcMoves;
//    mcmcMoves.AddMove(MultiplierMCMCMove(treeLike), 4.0);
//    mcmcMoves.AddMove(NodeSliderMCMCMove(treeLike), 1.0);
//    mcmcMoves.AddMove(SlidingWindowMCMCMove(treeLike), 1.0);
    mcmcMoves.AddMove(LazarusMCMCMove(treeLike, sampler,leaf_ids), 3.0);
    
    smc::moveset<TreeParticle> moveSet(particleInitializer, moveSelector, smcMoves, mcmcMoves);
    moveSet.SetNumberOfMCMCMoves(mcmcCount.getValue());

    // Output
    Json::Value jsonRoot;
    Json::Value& jsonTrees = jsonRoot["trees"];
    Json::Value& jsonIters = jsonRoot["generations"];
    if(jsonOutputPath.isSet()) {
        Json::Value& v = jsonRoot["run"];
        v["nQuerySeqs"] = static_cast<unsigned int>(query.getNumberOfSequences());
        v["nParticles"] = static_cast<unsigned int>(sampler.GetNumber());
        for(size_t i = 0; i < argc; i++)
            v["args"][i] = argv[i];
        v["version"] = sts::STS_VERSION;
        if(v["seed"].isNull()) v["seed"] = static_cast<unsigned int>(seed);
    }

    sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED, resample_threshold.getValue());
    sampler.SetMoveSet(moveSet);
    sampler.Initialise();
    const size_t nIters = (1 + treeMoveCount) * query.getNumberOfSequences();
    vector<string> sequenceNames = query.getSequencesNames();

    smc::DatabaseHistory database_history;

    for(size_t n = 0; n < nIters; n++) {
        double ess = 0.0;

        if (fribbleResampling.getValue()) {
            ess = sampler.IterateEssVariable(&database_history);
        } else {
            ess = sampler.IterateEss();
        }

        cerr << "Iter " << n << ": ESS=" << ess << " sequence=" << sequenceNames[n / (1 + treeMoveCount)] << endl;
        if(jsonOutputPath.isSet()) {
            Json::Value& v = jsonIters[n];
            v["T"] = static_cast<unsigned int>(n + 1);
            v["ess"] = ess;
            v["sequence"] = sequenceNames[n / (1 + treeMoveCount)];
            //v["totalUpdatePartialsCalls"] = static_cast<unsigned int>(BeagleTreeLikelihood::totalBeagleUpdateTransitionsCalls());
            v["totalUpdatePartialsCalls"] = static_cast<unsigned int>(AbstractFlexibleTreeLikelihood::operationCallCount);
            if (fribbleResampling.getValue()) {
                Json::Value ess_array;
                for (size_t i = 0; i < database_history.ess.size(); ++i)
                    ess_array.append(database_history.ess[i]);
                v["essHistory"] = ess_array;
            }
            
            std::set<size_t> set;
            for(size_t i = 0; i < sampler.GetNumber(); i++) {
                const TreeParticle& p = sampler.GetParticleValue(i);
                set.insert(p.particleID);
            }
            v["uniqueParticles"] = static_cast<unsigned int>(set.size());
        }
    }

    double maxLogLike = -std::numeric_limits<double>::max();
    for(size_t i = 0; i < sampler.GetNumber(); i++) {
        const TreeParticle& p = sampler.GetParticleValue(i);
//        treeLike.initialize(*p.model, *p.rateDist, *p.tree);
//        const double logLike = beagleLike->calculateLogLikelihood();
//        maxLogLike = std::max(logLike, maxLogLike);
        string s = bpp::TreeTemplateTools::treeToParenthesis(*p.tree);
        if(jsonOutputPath.isSet()) {
            Json::Value& v = jsonTrees[jsonTrees.size()];
//            v["treeLogLikelihood"] = logLike;
//            v["totalLikelihood"] = treeLike();
            v["particleID"] = static_cast<unsigned int>(p.particleID);
            v["newickString"] = s;
            v["logWeight"] = sampler.GetParticleLogWeight(i);
            v["treeLength"] = p.tree->getTotalLength();
        }
    }

    std::vector<ProposalRecord> proposalRecords = onlineAddSequenceMove->getProposalRecords();
    Json::Value& jsonProposals = jsonRoot["proposals"];
    for (size_t i = 0; i < proposalRecords.size(); ++i) {
        const auto& pr = proposalRecords[i];
        Json::Value& v = jsonProposals[i];
        v["T"] = static_cast<unsigned int>(pr.T);
        v["originalLogLike"] = pr.originalLogLike;
        v["newLogLike"] = pr.newLogLike;
        v["originalLogWeight"] = pr.originalLogWeight;
        v["newLogWeight"] = pr.newLogWeight;
        v["distalBranchLength"] = pr.proposal.distalBranchLength;
        v["distalLogProposalDensity"] = pr.proposal.distalLogProposalDensity;
        v["pendantBranchLength"] = pr.proposal.pendantBranchLength;
        v["pendantLogProposalDensity"] = pr.proposal.pendantLogProposalDensity;
        v["edgeLogProposalDensity"] = pr.proposal.edgeLogProposalDensity;
        v["logProposalDensity"] = pr.proposal.logProposalDensity();
        v["mlDistalBranchLength"] = pr.proposal.mlDistalBranchLength;
        v["mlPendantBranchLength"] = pr.proposal.mlPendantBranchLength;

        v["proposalMethodName"] = pr.proposal.proposalMethodName;
        
        if(pr.T == nIters){
            maxLogLike = std::max(pr.newLogLike, maxLogLike);
        }
    }

    if(jsonOutputPath.isSet()) {
        ofstream jsonOutput(jsonOutputPath.getValue());
        Json::StyledWriter writer;
        jsonOutput << writer.write(jsonRoot);
    }

    if(particleGraphPath.isSet()) {
        ofstream gOut(particleGraphPath.getValue());
        sampler.StreamParticleGraph(gOut);
    }

    clog << "Maximum LL: " << maxLogLike << '\n';

    gsl_rng_free(rng);
}
