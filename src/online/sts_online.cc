#include <Bpp/Phyl/Io/NexusIoTree.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/HKY85.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/Protein/WAG01.h>
#include <Bpp/Phyl/Model/Protein/LG08.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/DNA.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>

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
#include "multiplier_mcmc_move.h"
#include "node_slider_mcmc_move.h"
#include "delta_exchange_mcmc_move.h"
#include "local_mcmc_move.h"
#include "snni_mcmc_move.h"
#include "multiplier_smc_move.h"
#include "node_slider_smc_move.h"
#include "tree_particle.h"
#include "weighted_selector.h"
#include "util.h"

#include "proposal_guided_parsimony.h"
#include "flexible_parsimony.h"
#include "node_sliding_window_mcmc_move.h"
#include <boost/algorithm/string.hpp>
#include <boost/range/combine.hpp>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

#include "dirichlet_prior.h"

#include "scale_mcmc_move.h"
#include "transform.h"
#include "adaptive_mcmc_move.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace cl = TCLAP;
using namespace std;
typedef bpp::TreeTemplate<bpp::Node> Tree;

using namespace sts::online;


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

std::vector<std::vector<double>> readParameters(std::string path, const std::vector<string>& p, unsigned burnin){
    ifstream inputStream;
    string buffer;
    inputStream.open(path);
    getline(inputStream, buffer, '\n');
    
    getline(inputStream, buffer, '\n');
    std::vector<std::string> strs;
    boost::split(strs, buffer, boost::is_any_of("\t ,"));

    std::vector<int> map;
    map.resize(p.size());
    std::transform(p.cbegin(), p.cend(), map.begin(),
                   [strs](const std::string& name){return static_cast<int>(find(strs.begin(), strs.end(), name) - strs.begin());});
    
    std::vector<std::vector<double>> params;
    int count = 0;
    while ( !inputStream.eof() ) {
        getline(inputStream, buffer, '\n');
        if(buffer.size() == 0) continue;
        if(count >= burnin){
            boost::split(strs, buffer, boost::is_any_of("\t ,"));
            std::vector<double> sample;
            for(int index: map){
                sample.push_back(stod(strs[index]));
            }
            params.push_back(sample);
        }
        count++;
    }
    inputStream.close();
    return params;
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
    
    cl::UnlabeledValueArg<string> jsonOutputPath("json_path", "JSON output path", false, "", "path", cmd);
    
    cl::ValueArg<string> paramsPath("P", "parameters", "Parameters input path", false, "", "path", cmd);
	std::vector<std::string> modelNames { "JC69", "K80", "HKY", "GTR", "LG", "WAG" };
	cl::ValuesConstraint<std::string> allowedModels(modelNames);
    cl::ValueArg<string> modelArg("M", "model", "Substitution model", false, "JC69", &allowedModels, cmd);
    cl::ValueArg<int> catCountArg("c", "categories", "Number of categories for Gamma distribution ", false, 1, "#", cmd);
    cl::SwitchArg topoPriorArg("", "uniform-topology", "Use uniform prior on topologies", cmd, false);
    
    cl::ValueArg<long> seedCmd("s", "seed", "Seed for random number generator", false, -1, "#", cmd);
    
    cl::ValueArg<double> exponent("e", "exponent", "Exponent for calculting propbablity attachments in step 1", false, 0.05, "#", cmd);

    cl::ValueArg<string> stemArg("o", "stem", "Stem to output posterior samples", false, "", "stem", cmd);
	
	cl::SwitchArg mvnArg("G", "mvn", "Use a multivariate normal proposal", cmd, false);
	cl::SwitchArg uniProposalArg("u", "uni", "Use a univariate proposal for each parameter", cmd, true);
#if defined(_OPENMP)
	cl::ValueArg<int> threadCountArg("T", "threads", "Number of threads", false, 1, "#", cmd);
#endif
	
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

	std::string modelString = modelArg.getValue();
	unique_ptr<bpp::Alphabet> alphabet;	
	if(modelString == "LG" || modelString == "WAG"){
		alphabet.reset(new bpp::ProteicAlphabet());
	}
	else{
		alphabet.reset(new bpp::DNA());
	}
	
    unique_ptr<bpp::GeneticCode> gCode;
    //    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    //    if (codonAlphabet) {
    //        string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml.getParams(), "Standard", "", true, true);
    //        ApplicationTools::displayResult("Genetic Code", codeDesc);
    //
    //        gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    //    }
    
    ifstream alignment_fp(alignmentPath.getValue());
    unique_ptr<bpp::SiteContainer> sites(sts::util::read_alignment(alignment_fp, alphabet.get()));
    alignment_fp.close();
    bpp::VectorSiteContainer ref(alphabet.get()), query(alphabet.get());
    partitionAlignment(*sites, trees[0]->getLeavesNames(), ref, query);
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
    
    // TODO: Other distributions
    const double expPriorMean = blPriorExpMean.getValue();
    std::function<double(double)> exponentialPrior = [expPriorMean](const double d) {
        return std::log(gsl_ran_exponential_pdf(d, expPriorMean));
    };

    vector<TreeParticle> particles;
    particles.reserve(trees.size());

    int catCount = catCountArg.getValue();
    
    bpp::SubstitutionModel*  model = nullptr;
//    bpp::SubstitutionModelSet* modelSet = nullptr;
    bpp::DiscreteDistribution* rate_dist    = nullptr;
    
    std::vector<std::string> mrbayesNames;
    std::map<std::string, size_t> paramNames;
    std::vector<std::vector<double>> params;
    if(paramsPath.isSet()){
        if(modelString == "K80"){
            paramNames["kappa"] = 0;
            mrbayesNames.push_back("kappa");
        }
        else if(modelString == "HKY"){
            size_t N = 4;
            const char* vinit[] = {"kappa", "pi(A)", "pi(C)", "pi(G)", "pi(T)"};
            for (size_t i = 0; i < N; i++) {
                mrbayesNames.push_back(vinit[i]);
                paramNames[vinit[i]] = i;
            }
        }
        else if(modelString == "GTR"){
            size_t N = 10;
            const char* vinit[] = {"r(C<->T)", "r(A<->T)", "r(G<->T)", "r(A<->C)", "r(C<->G)", "r(A<->G)", "pi(A)", "pi(C)", "pi(G)", "pi(T)"};
            for (size_t i = 0; i < N; i++) {
                mrbayesNames.push_back(vinit[i]);
                paramNames[vinit[i]] = i;
            }
        }
        if(catCount > 1){
			// Do not assign paramNames.size() directly to paramNames["alpha"]
			size_t size = paramNames.size();
			paramNames["alpha"] = size;
			mrbayesNames.push_back("alpha");
        }
    }
    if(paramNames.size() > 0){
        params = readParameters(paramsPath.getValue(), mrbayesNames, burnin.getValue());
        cout << "read "<< params.size()<<" samples"<<endl;
        assert(params.size() == trees.size());
    }
    
    std::vector<std::string> modelParameterNames;
    
    //for (auto tup : boost::combine(params, trees)) {
    for(int i = 0; i < trees.size(); i++){
        unique_ptr<Tree>& tree = trees[i];
        tree->getRootNode()->getSon(0)->setDistanceToFather(tree->getRootNode()->getSon(0)->getDistanceToFather() +
                                                            tree->getRootNode()->getSon(1)->getDistanceToFather());
        tree->getRootNode()->getSon(1)->setDistanceToFather(0.0);
        
        if(paramsPath.isSet()){
			if(dynamic_cast<bpp::NucleicAlphabet*>(alphabet.get()) != nullptr){
				bpp::NucleicAlphabet* nucAlphabet = dynamic_cast<bpp::NucleicAlphabet*>(alphabet.get());
				if(modelString == "K80"){
					model = new bpp::K80(nucAlphabet, params[i][paramNames["kappa"]]);
                    modelParameterNames.push_back("kappa");
				}
				else if(modelString == "HKY"){
					model = new bpp::HKY85(nucAlphabet, params[i][paramNames["kappa"]], params[i][paramNames["pi(A)"]], params[i][paramNames["pi(C)"]], params[i][paramNames["pi(G)"]], params[i][paramNames["pi(T)"]]);
                    const char *vinit[] = {"kappa", "theta", "theta1", "theta2"};
                    modelParameterNames.assign(vinit, vinit+4);
				}
				else if(modelString == "GTR"){
                    double f = params[i][paramNames["r(A<->G)"]];
					double a = params[i][paramNames["r(C<->T)"]]/f;
					double b = params[i][paramNames["r(A<->T)"]]/f;
					double c = params[i][paramNames["r(G<->T)"]]/f;
					double d = params[i][paramNames["r(A<->C)"]]/f;
					double e = params[i][paramNames["r(C<->G)"]]/f;
					
					model = new bpp::GTR(nucAlphabet, a, b, c, d, e, params[i][paramNames["pi(A)"]], params[i][paramNames["pi(C)"]], params[i][paramNames["pi(G)"]], params[i][paramNames["pi(T)"]]);
                    const char *vinit[] = {"a", "b", "c", "d", "e", "theta", "theta1", "theta2"};
                    modelParameterNames.assign(vinit, vinit+8);
				}
				else{
					model = new bpp::JCnuc(nucAlphabet);
				}
			}
			else if(dynamic_cast<bpp::ProteicAlphabet*>(alphabet.get()) != nullptr){
				bpp::ProteicAlphabet* protAlphabet = dynamic_cast<bpp::ProteicAlphabet*>(alphabet.get());
				if(modelString == "LG"){
					model = new bpp::LG08(protAlphabet);
				}
				else if(modelString == "WAG"){
					model = new bpp::WAG01(protAlphabet);
				}
			}
        }
		else{
			bpp::NucleicAlphabet* nucAlphabet = dynamic_cast<bpp::NucleicAlphabet*>(alphabet.get());
            if(nucAlphabet){
                model = new bpp::JCnuc(nucAlphabet);
                rate_dist = new bpp::ConstantRateDistribution();
            }
            else{
                bpp::ProteicAlphabet* protAlphabet = dynamic_cast<bpp::ProteicAlphabet*>(alphabet.get());
                if(modelString == "LG"){
                    model = new bpp::LG08(protAlphabet);
                }
                else if(modelString == "WAG"){
                    model = new bpp::WAG01(protAlphabet);
                }
            }
            
        }
        
        if(catCount == 1){
            rate_dist = new bpp::ConstantRateDistribution();
        }
        else{
            rate_dist = new bpp::GammaDiscreteRateDistribution(catCount, params[i][paramNames["alpha"]]);
        }
        
        particles.emplace_back(std::unique_ptr<bpp::SubstitutionModel>(model),
                               std::unique_ptr<bpp::TreeTemplate<bpp::Node>>(tree.release()),
                               std::unique_ptr<bpp::DiscreteDistribution>(rate_dist),
                               &ref);
    }

    std::unique_ptr<bpp::SitePatterns> _patterns(new bpp::SitePatterns(sites.get()));

	size_t threadCount = 1;
#if defined(_OPENMP)
	threadCount = threadCountArg.getValue();
#endif
	std::vector<std::unique_ptr<CompositeTreeLikelihood>> treeLikes;
	for(size_t i = 0; i < threadCount; i++){
#ifndef NO_BEAGLE
		shared_ptr<FlexibleTreeLikelihood> beagleLike(new BeagleFlexibleTreeLikelihood(*_patterns.get(), *particles[0].model, *particles[0].rateDist));
#else
		shared_ptr<FlexibleTreeLikelihood> beagleLike(new SimpleFlexibleTreeLikelihood(*_patterns.get(), *particles[0].model, *particles[0].rateDist));
#endif
		CompositeTreeLikelihood* treeLike = new CompositeTreeLikelihood(beagleLike);
		treeLike->add(BranchLengthPrior(exponentialPrior));
		treeLikes.push_back(std::unique_ptr<CompositeTreeLikelihood>(treeLike));
	
	// Kappa prior
    if(modelString == "K80" || modelString == "HKY"){
		std::function<double(double)> ratesGammaPrior = [](const double d) {
			return std::log(gsl_ran_gamma_pdf(d, 0.05, 10));
		};
		std::vector<std::string> rr{modelString+".kappa"};
		std::unique_ptr<Prior> ratePrior(new Prior(rr, ratesGammaPrior));
		treeLike->add(ratePrior);
    }
	// Frequencies prior
    if(modelString == "HKY" || modelString == "GTR"){
//        std::function<double()> freqsFlatDirichlet = []() {
//            return gsl_sf_lngamma(4);
//        };
        std::vector<std::string> ff{modelString+".pA",modelString+".pC",modelString+".pG",modelString+".pT"};
        std::unique_ptr<Prior> freqsPrior(new DirichletPrior(ff));
        
        treeLike->add(freqsPrior);
    }
	// GTR relative rate prior
	if (modelString == "GTR") {
		std::function<double(double)> ratesGammaPrior = [](const double d) {
			return std::log(gsl_ran_gamma_pdf(d, 0.05, 10));
		};
		//std::unique_ptr<Prior> prior(new Prior(model->getParameters().getParameterNames(), ratesGammaPrior));
		std::vector<std::string> rr{modelString+".a",modelString+".b",modelString+".c",modelString+".d",modelString+".e"};
		std::unique_ptr<Prior> ratesPrior(new Prior(rr, ratesGammaPrior));
		treeLike->add(ratesPrior);
	}
	// Prior for rate heterogenity across sites (alpha)
	if(catCount > 1){
		std::function<double(double)> alphaGammaPrior = [](const double d) {
			return std::log(gsl_ran_gamma_pdf(d, 0.05, 10));
		};
        std::function<double(double)> alphaExpPrior = [](const double d) {
            return std::log(gsl_ran_exponential_pdf(d, 1.0));
        };
		std::vector<std::string> alpha{"alpha"};
        std::unique_ptr<Prior> alphaPrior(new Prior(alpha, alphaExpPrior));
		treeLike->add(alphaPrior);
	}
    // Uniform prior on topologies
    if(topoPriorArg.getValue()){
        double topologyLogP = -std::log(gsl_sf_doublefact(sites->getNumberOfSequences()*2-5));
        
        std::function<double(double)> uniformTopologyPrior = [topologyLogP](const double d) {
            return topologyLogP;
        };
        std::vector<std::string> empty{};
        std::unique_ptr<Prior> topoPrior(new Prior({}, uniformTopologyPrior));
        treeLike->add(topoPrior);
    }
}

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
            onlineAddSequenceMove.reset(new UniformLengthOnlineAddSequenceMove(treeLikes, sites->getSequencesNames(), query.getSequencesNames(), branchLengthProposer));
        } else {
            onlineAddSequenceMove.reset(new UniformOnlineAddSequenceMove(treeLikes, sites->getSequencesNames(), query.getSequencesNames(), branchLengthProposer));
        }
    } else{
        GuidedOnlineAddSequenceMove* p = nullptr;
        
        if(name == "guided") {
            p = new GuidedOnlineAddSequenceMove(treeLikes, sites->getSequencesNames(), query.getSequencesNames(), pbl, maxLength.getValue(), subdivideTop.getValue());
        } else if(name == "lcfit") {
            p = new LcfitOnlineAddSequenceMove(treeLikes, sites->getSequencesNames(), query.getSequencesNames(), pbl, maxLength.getValue(), subdivideTop.getValue(), expPriorMean);
        } else if(name == "guided-parsimony") {
            std::shared_ptr<FlexibleParsimony> pars = make_shared<FlexibleParsimony>(*_patterns.get(), *alphabet);
            p = new ProposalGuidedParsimony(pars, treeLikes, sites->getSequencesNames(), query.getSequencesNames(), expPriorMean);
        }
        else{
            throw std::runtime_error("Unknown sequence addition method: " + name);
        }
        p->_heating = exponent.getValue();
        onlineAddSequenceMove.reset(p);
    }
	
		
	if(false)
    if(catCount > 1){
        double firstMoment = 0;
        double secondMoment = 0;
        for(auto& p : params){
            double value = p[paramNames["alpha"]];
            firstMoment += value;
            secondMoment += value*value;
        }
        firstMoment /= params.size();
        secondMoment /= params.size();
        double firstMoment2 = firstMoment*firstMoment;
        double shape = firstMoment2/(secondMoment - firstMoment2);
        double scale = (secondMoment - firstMoment2)/firstMoment;
        
        std::cout << "shape: " << shape << " scale: " << scale <<endl;
        
        std::function<std::tuple<double, double>(smc::rng*)> empiricalGammaProposal = [shape, scale](smc::rng* rng) {
            double alpha = gsl_ran_gamma(rng->GetRaw(), shape, scale);
            double alphaLogP = gsl_ran_gamma_pdf(alpha, shape, scale);
            return std::make_tuple(alpha, alphaLogP);
        };
        onlineAddSequenceMove->setGammaProposal(empiricalGammaProposal);
    }
    
    std::unique_ptr<gsl_matrix, decltype(&gsl_matrix_free)> L(gsl_matrix_alloc(1,1), gsl_matrix_free);
    std::unique_ptr<gsl_vector, decltype(&gsl_vector_free)> mu(gsl_vector_alloc(1), gsl_vector_free);
    
    if(modelString == "GTR" || modelString == "HKY"){
        size_t dimRates = 5;
        if(modelString == "HKY"){
            dimRates = 1;
        }
        size_t dim = 3 + dimRates;
        size_t sampleCount = params.size();
        mu = std::unique_ptr<gsl_vector, decltype(&gsl_vector_free)>(gsl_vector_alloc(dim), gsl_vector_free);
        L = std::unique_ptr<gsl_matrix, decltype(&gsl_matrix_free)>(gsl_matrix_alloc(dim, dim), gsl_matrix_free);
//        std::unique_ptr<gsl_matrix, decltype(&gsl_matrix_free)> L(gsl_matrix_alloc(dim, dim), gsl_matrix_free);
//        std::unique_ptr<gsl_vector, decltype(&gsl_vector_free)> mu(gsl_vector_alloc(dim), gsl_vector_free);
        std::vector<std::vector<double>> paramsT;
        paramsT.resize(dim);
        // log transform relative rates or kappa
        if(modelString == "GTR"){
            for(size_t i = 0; i < 5; i++){
                for(size_t j = 0; j < sampleCount; j++){
                    paramsT[i].push_back(std::log(params[j][i]/params[j][5]));
                }
            }
        }
        // log transform kappa
        else{
            for(size_t j = 0; j < sampleCount; j++){
                paramsT[0].push_back(std::log(params[j][0]));
            }
        }
        
        //Transform frequencies
        for (size_t i = 0; i < 3; i++) {
            for(size_t k = 0; k < sampleCount; k++){
                double sum = 1;
                for (size_t j = 0; j < i; j++) {
                    sum -= params[k][dimRates+j+1]; // +1 because there are 6 rate parameters in GTR in Mrbayes
                }
                paramsT[dimRates+i].push_back(sts::util::logit(params[k][dimRates+i+1]/sum) - std::log(1.0/(3-i)));
            }
        }

        // means
        for(size_t i = 0; i < dim; i++){
            double mean = std::accumulate(paramsT[i].cbegin(), paramsT[i].cend(), 0.0);
            gsl_vector_set(mu.get(), i, mean/sampleCount);
        }

        // covariance matrix
        for(size_t i = 0; i < dim; i++){
            double mu1 = gsl_vector_get(mu.get(), i);
            double var = 0;
            for(size_t k = 0; k < sampleCount; k++){
                var += pow(mu1-paramsT[i][k], 2);
            }
            gsl_matrix_set(L.get(), i, i, var/(sampleCount-1));
            
            for(size_t j = i+1; j < dim; j++){
                double mu2 = gsl_vector_get(mu.get(), j);
                double cov = 0;
                for(size_t k = 0; k < sampleCount; k++){
                    cov += (paramsT[i][k]-mu1)*(paramsT[j][k]-mu2);
                }
                gsl_matrix_set(L.get(), i, j, cov/(sampleCount-1));
                gsl_matrix_set(L.get(), j, i, cov/(sampleCount-1));
            }
        }
        gsl_linalg_cholesky_decomp1(L.get());
        
        gsl_vector* muptr = mu.get();
        gsl_matrix* Lptr = L.get();
        
        // C++14
        //std::function<std::tuple<std::vector<double>, double>(smc::rng*)> empiricalMVNProposal = [mu=move(mu), L=move(L)](smc::rng* rng) {
        std::function<std::tuple<std::map<std::string, double>, double>(smc::rng*)> empiricalMVNProposal = [muptr, Lptr, modelParameterNames, dimRates](smc::rng* rng) {
            size_t dim = muptr->size;
            gsl_vector* result = gsl_vector_alloc(dim);
            gsl_vector* work = gsl_vector_alloc(dim);
            gsl_ran_multivariate_gaussian(rng->GetRaw(), muptr, Lptr, result);
            std::map<std::string, double> x;
            double jacobian = 0;
            for(size_t i = 0; i < dimRates; i++){
                x[modelParameterNames[i]] = std::exp(gsl_vector_get(result, i));
                jacobian -= gsl_vector_get(result, i);
            }
            
            vector<double> xx(4, 0);
            xx[3] = 1.0;
            for(size_t i = 0; i < 3; i++){
                double zi = sts::util::logitinv(gsl_vector_get(result, dimRates+i) + std::log(1.0/(3-i)));
                double sum = 0.0;
                for(size_t j = 0; j < i; j++){
                    sum += xx[j];
                }
                xx[i] = (1.0-sum)*zi;
                double logx = log(xx[i]);
                jacobian += log(-1.0/((1.0 - sum)*(xx[i]*xx[i] - xx[i])));
                xx[3] -= xx[i];
            }

            x["theta"] = xx[1] + xx[2];
            x["theta1"] = xx[0]/(xx[0] + xx[3]);
            x["theta2"] = xx[2]/(xx[1] + xx[2]);
            
            double logP = 0;
            gsl_ran_multivariate_gaussian_log_pdf(result, muptr, Lptr, &logP, work);
            
            logP += jacobian;
            gsl_vector_free(work);
            gsl_vector_free(result);
            return std::make_tuple(x, logP);
        };
        
        onlineAddSequenceMove->setMVNProposal(empiricalMVNProposal);
    }
	
    {
        using namespace std::placeholders;
        auto wrapper = std::bind(&OnlineAddSequenceMove::operator(), std::ref(*onlineAddSequenceMove), _1, _2, _3);
        smcMoves.push_back(wrapper);
    }

    if(treeMoveCount) {
        smcMoves.push_back(MultiplierSMCMove(*treeLikes[0]));
        smcMoves.push_back(NodeSliderSMCMove(*treeLikes[0]));
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

    smc::sampler<TreeParticle> sampler(particleFactor.getValue() * trees.size(), SMC_HISTORY_NONE, gsl_rng_default, seed);
    smc::mcmc_moves<TreeParticle> mcmcMoves;
	LocalMCMCMove localMove(treeLikes);
	NNIMCMCMove snniMove(treeLikes);
    MultiplierMCMCMove multMove(treeLikes);
//    NodeSliderMCMCMove sliderMove(treeLike);
//    SlidingWindowMCMCMove slidingMove(treeLike);
    mcmcMoves.AddMove(multMove, 4.0);
//    mcmcMoves.AddMove(sliderMove, 1.0);
//    mcmcMoves.AddMove(slidingMove, 1.0);
	mcmcMoves.AddMove(localMove, 10.0);
	mcmcMoves.AddMove(snniMove, 1.0);
	
	if(uniProposalArg.getValue()){
		if(model->getNumberOfParameters() > 0){
			for(string& p : model->getParameters().getParameterNames()){
				if(model->getParameterNameWithoutNamespace(p).size() == 1){
					MultiplierMCMCMove multMove(treeLikes, {model->getParameterNameWithoutNamespace(p)});
					mcmcMoves.AddMove(multMove, 1.0);
				}
			}
			if(model->getParameters().hasParameter(model->getNamespace() + "theta")){
				DeltaExchangeMCMCMove deltaMove(treeLikes, {"theta", "theta1", "theta2"}, 0.2);
				mcmcMoves.AddMove(deltaMove, 1.0);
			}
		}
		if(rate_dist->getNumberOfParameters() > 0){
			MultiplierMCMCMove multMove(treeLikes, {"alpha"});
			mcmcMoves.AddMove(multMove, 1.0);
		}
	}

	size_t dimRates = 0;
	size_t dimFreqs = 0;
	if(modelString == "GTR" || modelString == "HKY" || modelString == "K80"){
		dimRates = 5;
		if(modelString == "HKY"|| modelString == "K80"){
			dimRates = 1;
		}
		if(modelString == "GTR" || modelString == "HKY"){
			dimFreqs = 3;
		}
	}
	size_t dim = dimRates + dimFreqs;
	if(catCount > 1){
		dim++;
	}
	if(dim  < 2){
		std::cerr << "STS cannot use a multivariate proposal with less than 2 free parameters" << std::endl;
		exit(1);
	}
	
	if(mvnArg.getValue()){
		size_t sampleCount = params.size();
		std::vector<std::vector<double>> paramsT;
		paramsT.resize(dim);
		
		gsl_vector* mu = gsl_vector_alloc(dim);
		gsl_matrix* L = gsl_matrix_alloc(dim, dim);
		
		if(dimRates > 0){
			// log transform relative rates or kappa
			if(modelString == "GTR"){
				vector<double> means(5,0);
				for(size_t i = 0; i < 5; i++){
					for(size_t j = 0; j < sampleCount; j++){
						paramsT[i].push_back(std::log(params[j][i]/params[j][5]));
						means[i] += params[j][i];
					}
				}
			}
			// log transform kappa
			else{
				for(size_t j = 0; j < sampleCount; j++){
					paramsT[0].push_back(std::log(params[j][0]));
				}
			}
		}
		
		if(dimFreqs > 0){
			size_t indexParamFreqs = dimRates;
			 // +1 because there are 6 rate parameters in GTR in Mrbayes
			if(modelString == "GTR") indexParamFreqs++;
			//Transform frequencies
			for (size_t i = 0; i < 3; i++) {
				for(size_t k = 0; k < sampleCount; k++){
					double sum = 1;
					for (size_t j = 0; j < i; j++) {
						sum -= params[k][indexParamFreqs+j];
					}
					paramsT[dimRates+i].push_back(sts::util::logit(params[k][indexParamFreqs+i]/sum) - std::log(1.0/(3-i)));
				}
			}
		}
		
		if(catCount > 1){
			size_t indexParamAlpha = dimRates + dimFreqs;
			if(dimFreqs > 0) indexParamAlpha++;
			if(modelString == "GTR") indexParamAlpha++;
			
			for(size_t j = 0; j < sampleCount; j++){
				paramsT[dimRates+dimFreqs].push_back(std::log(params[j][indexParamAlpha]));
			}
		}
		
		// means
		for(size_t i = 0; i < dim; i++){
			double mean = std::accumulate(paramsT[i].cbegin(), paramsT[i].cend(), 0.0);
			gsl_vector_set(mu, i, mean/sampleCount);
		}
		
		// covariance matrix
		for(size_t i = 0; i < dim; i++){
			double mu1 = gsl_vector_get(mu, i);
			double var = 0;
			for(size_t k = 0; k < sampleCount; k++){
				var += pow(mu1-paramsT[i][k], 2);
			}
			gsl_matrix_set(L, i, i, var/(sampleCount-1));
			
			for(size_t j = i+1; j < dim; j++){
				double mu2 = gsl_vector_get(mu, j);
				double cov = 0;
				for(size_t k = 0; k < sampleCount; k++){
					cov += (paramsT[i][k]-mu1)*(paramsT[j][k]-mu2);
				}
				gsl_matrix_set(L, i, j, cov/(sampleCount-1));
				gsl_matrix_set(L, j, i, cov/(sampleCount-1));
			}
		}
		gsl_linalg_cholesky_decomp1(L);
		
		// Add the transforms in the same order as the L and mu matrices
		// Add rate transforms
		std::vector<Transform*> transforms;
		for(string& p : model->getParameters().getParameterNames()){
			std::string p2 = model->getParameterNameWithoutNamespace(p);
			if(p2 != "theta" && p2 != "theta1" && p2 != "theta2" && p2 != "alpha"){
				LogTransform* transform = new LogTransform({model->getParameterNameWithoutNamespace(p)});
				transforms.push_back(transform);
			}
		}
		// Add frequencies transform
		if(model->getParameters().hasParameter(model->getNamespace() + "theta")){
			SimplexTransform* transform = new SimplexTransform({"theta", "theta1", "theta2"});
			transforms.push_back(transform);
		}
		// Add alpha transform
		if(catCount > 1){
			LogTransform* transform = new LogTransform({"alpha"});
			transforms.push_back(transform);
		}
		AdaptiveMCMCMove adaptMove(treeLikes, *L, *mu, transforms);
		mcmcMoves.AddMove(adaptMove, 1.0);
	}
	
    smc::moveset<TreeParticle> moveSet(particleInitializer, moveSelector, smcMoves, mcmcMoves);
    moveSet.SetNumberOfMCMCMoves(mcmcCount.getValue());

    // Output
    Json::Value jsonRoot;
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
#if defined(_OPENMP)
	sampler.SetNumberOfThreads(threadCount);
#endif
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
	
	if(stemArg.isSet()){
		// output trees
		ofstream treesOutput(stemArg.getValue() + ".trees");
		// output parameters as a csv file
		ofstream logOutput(stemArg.getValue() + ".log");
		logOutput << "particle\tll\tTL";
		
		if(modelString == "GTR"){
			logOutput << "\tr(A<->C)\tr(A<->G)\tr(A<->T)\tr(C<->G)\tr(C<->T)\tr(G<->T)";
		}else
		for (const std::string& paramName: model->getParameters().getParameterNames()) {
			string name = model->getParameterNameWithoutNamespace(paramName);
			if(name != "theta" && name != "theta1" && name != "theta2"){
				logOutput << "\t" << paramName;
			}
		}
		
		for (const std::string& paramName: model->getParameters().getParameterNames()) {
			string name = model->getParameterNameWithoutNamespace(paramName);
			if(name == "theta" || name == "theta1" || name == "theta2"){
				logOutput << "\tpi(A)\tpi(C)\tpi(G)\tpi(T)";
				break;
			}
		}
		if(catCount > 1) logOutput << "\talpha";
		
		logOutput << endl;
		for(size_t i = 0; i < sampler.GetNumber(); i++) {
			const TreeParticle& p = sampler.GetParticleValue(i);
			
			maxLogLike = std::max(p.logP, maxLogLike);
			
			string s = bpp::TreeTemplateTools::treeToParenthesis(*p.tree);
			treesOutput << s;
			
			logOutput << i << "\t" << p.logP << "\t" << p.tree->getTotalLength();
		
			if(modelString == "GTR"){
				double sum = 1;
				sum += p.model->getParameterValue("a");
				sum += p.model->getParameterValue("b");
				sum += p.model->getParameterValue("c");
				sum += p.model->getParameterValue("d");
				sum += p.model->getParameterValue("e");
				logOutput << "\t" << p.model->getParameterValue("d")/sum;
				logOutput << "\t" << 1.0/sum;//f
				logOutput << "\t" << p.model->getParameterValue("b")/sum;
				logOutput << "\t" << p.model->getParameterValue("e")/sum;
				logOutput << "\t" << p.model->getParameterValue("a")/sum;
				logOutput << "\t" << p.model->getParameterValue("c")/sum;
			}
			else{
				for (const std::string& paramName: p.model->getParameters().getParameterNames()) {
					string name = p.model->getParameterNameWithoutNamespace(paramName);
					if(name != "theta" && name != "theta1" && name != "theta2"){
						logOutput << "\t" << p.model->getParameterValue(name);
					}
				}
			}
		
			for (const std::string& paramName: p.model->getParameters().getParameterNames()) {
				string name = p.model->getParameterNameWithoutNamespace(paramName);
				if(name == "theta" || name == "theta1" || name == "theta2"){
					const std::vector<double>& frequencies = p.model->getFrequencies();
					for(auto freq : frequencies){
						logOutput << "\t" << freq;
					}
					break;
				}
			}
		
			if(catCount > 1){
				logOutput << "\t" << p.rateDist->getParameterValue("alpha");
			}
			logOutput << endl;
		}
		logOutput.close();
		treesOutput.close();
	}
	
	if(jsonOutputPath.isSet()) {
		Json::Value& jsonTrees = jsonRoot["trees"];
		for(size_t i = 0; i < sampler.GetNumber(); i++) {
			const TreeParticle& p = sampler.GetParticleValue(i);
			string s = bpp::TreeTemplateTools::treeToParenthesis(*p.tree);
			Json::Value& v = jsonTrees[jsonTrees.size()];
			v["particleID"] = static_cast<unsigned int>(p.particleID);
			v["newickString"] = s;
			v["logWeight"] = sampler.GetParticleLogWeight(i);
			v["treeLength"] = p.tree->getTotalLength();
			if(catCount > 1) v["alpha"] = p.rateDist->getParameterValue("alpha");
			for (const std::string& paramName: p.model->getParameters().getParameterNames()) {
				string name = p.model->getParameterNameWithoutNamespace(paramName);
				if(name != "theta" && name != "theta1" && name != "theta2"){
					v[paramName] = p.model->getParameterValue(name);
				}
			}
			
			for (const std::string& paramName: p.model->getParameters().getParameterNames()) {
				string name = p.model->getParameterNameWithoutNamespace(paramName);
				if(name == "theta" || name == "theta1" || name == "theta2"){
					const std::vector<double>& frequencies = p.model->getFrequencies();
					v["piA"] = frequencies[0];
					v["piC"] = frequencies[1];
					v["piG"] = frequencies[2];
					v["piT"] = frequencies[3];
					break;
				}
			}
		}
	
		Json::Value& jsonProposals = jsonRoot["proposals"];
		std::vector<ProposalRecord> proposalRecords = onlineAddSequenceMove->getProposalRecords();
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
			v["substModelLogProposalDensity"] = pr.proposal.substModelLogProposalDensity;

			v["proposalMethodName"] = pr.proposal.proposalMethodName;
		}

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
