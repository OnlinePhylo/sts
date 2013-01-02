#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/NexusIoTree.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

#include <gsl/gsl_randist.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

using bpp::Node;
using bpp::TreeTemplate;
using bpp::TreeTemplateTools;

typedef vector<TreeTemplate<Node>> Forest;

// Global Random Number Generator
std::mt19937 RNG;

/// Split a tree on a random branch
template<typename Trng>
Forest random_split(const TreeTemplate<Node>& tree, Trng& rng)
{
    vector<const Node*> nodes = tree.getNodes();
    Forest result;

    // Choose a random node to be the root
    // Choose a number between 0 and nodes.size()-2 to skip the trifurcating current root node.
    std::uniform_int_distribution<size_t> idx_dist(0, nodes.size() - 2);
    size_t idx = idx_dist(rng);

    // If the selected index is trifurcating, choose the next node (guaranteed to exist from the bounds above)
    if(nodes[idx]->getNumberOfSons() == 3) idx++;

    // Copy and reroot
    TreeTemplate<Node> rerooted = tree;
    rerooted.rootAt(tree.getRootNode()->getId());
    rerooted.newOutGroup(nodes[idx]->getId());
    const Node* root = rerooted.getRootNode();

    // Should be bifurcating
    assert(root->getNumberOfSons() == 2);
    for(size_t i = 0; i < root->getNumberOfSons(); i++) {
        Node* node = TreeTemplateTools::cloneSubtree<Node>(*rerooted.getRootNode()->getSon(i));
        result.emplace_back(node);
    }

    return result;
}

unsigned int count_merged_leaves(const Forest& forest)
{
    auto f = [] (const unsigned acc, const TreeTemplate<Node>& t) {
        return acc + t.getNumberOfLeaves();
    };
    return accumulate(begin(forest), end(forest), 0, f);
}

/// Decompose a 4-taxon tree into a pair of cherries.
Forest cherries(TreeTemplate<Node> tree)
{
    assert(tree.getLeaves().size() == 4);
    TreeTemplate<Node> orig_tree = tree;

    // Returns whether all nodes immediately below this one are leaves.
    auto all_leaves = [] (const Node& node) -> bool {
        if(node.isLeaf()) return false;
        for(size_t i = 0; i < node.getNumberOfSons(); ++i) {
            if(!node.getSon(i)->isLeaf()) return false;
        }
        return true;
    };
    vector<Node*> nodes = orig_tree.getNodes();

    for(const Node* node : nodes) {
        if(node->getNumberOfSons() == 2 and all_leaves(*node)) {
            Forest result;
            tree.newOutGroup(node->getId());
            for(size_t i = 0; i < tree.getRootNode()->getNumberOfSons(); i++) {
                Node* node = bpp::TreeTemplateTools::cloneSubtree<Node>(*tree.getRootNode()->getSon(i));
                result.emplace_back(node);
            }
            return result;
        }
    }

    throw std::runtime_error("FALSE");
}

/// Get the log-likelihood of a sequence under the stationary frequency distribution
double sequence_log_likelihood(const bpp::Sequence& sequence, const bpp::SubstitutionModel* model, const bpp::DiscreteDistribution* rates)
{
    const unsigned int n_states = model->getNumberOfStates(),
          n_sites = sequence.size(),
          n_rates = rates->getNumberOfCategories();
    const vector<double>& freqs = model->getFrequencies();

    double log_like = 0;
    for(unsigned int rate = 0; rate < n_rates; rate++) {
        for(unsigned int site = 0; site < n_sites; site++) {
            double like = 0;
            for(unsigned int i = 0; i < n_states; i++)
                like += freqs[i] * model->getInitValue(i, sequence.getValue(site)) * rates->getProbability(rate);
            log_like += std::log(like);
        }
    }
    return log_like;
}

double forest_likelihood(const Forest& f, const bpp::SiteContainer& data, bpp::SubstitutionModel* model, bpp::DiscreteDistribution* rates)
{
    double result = 0;

    for (auto &tree : f) {
        if(tree.getNumberOfLeaves() == 1) {
            // Compute background frequencies for a leaf
            const Node* leaf = tree.getLeaves()[0];
            assert(data.hasSequence(leaf->getName()));
            result += sequence_log_likelihood(data.getSequence(data.getSequencePosition(leaf->getName())),
                    model,
                    rates);
        } else {
            bpp::RHomogeneousTreeLikelihood like(tree, data, model, rates, true, false);
            like.initialize();
            result += like.getLogLikelihood();
        }
    }

    return result;
}

double forest_length(Forest& f)
{
    auto tree_length = [](const double accum, TreeTemplate<Node>& tree) {
        return accum + tree.getTotalLength();
    };
    return std::accumulate(begin(f), end(f), 0.0, tree_length);
}

vector<unique_ptr<bpp::Tree>> read_nexus_trees(const string& path, const size_t burnin=0)
{
    bpp::NexusIOTree tree_io;
    vector<bpp::Tree*> raw_trees;
    tree_io.read(path, raw_trees);
    vector<unique_ptr<bpp::Tree>> trees;
    trees.reserve(raw_trees.size());
    for(bpp::Tree* t : raw_trees) trees.emplace_back(t);

    if(burnin) {
        assert(trees.size() > burnin);
        // Prune by moving into a new vector
        std::vector<unique_ptr<bpp::Tree>> result;
        result.reserve(trees.size() - burnin);
        auto e = end(trees);
        for(auto i = begin(trees) + burnin; i != e; ++i)
            result.push_back(std::move(*i));
        return result;
    }

    return trees;
}

struct Exponential_branch_prior
{
    double mu;

    double operator()(const Forest& f)
    {
        // Exponential BL prior on a single node
        auto pr_node = [this](const double accum, const Node* node) {
            return accum + node->hasDistanceToFather() ?
                std::log(gsl_ran_exponential_pdf(node->getDistanceToFather(), mu)) :
                0;
        };
        // Exponential BL prior on tree
        auto pr_tree = [&pr_node](const double accum, const TreeTemplate<Node>& tree) -> double {
            vector<const Node*> nodes = tree.getNodes();
            return accum + std::accumulate(begin(nodes), end(nodes), 0.0, pr_node);
        };
        return std::accumulate(begin(f), end(f), 0.0, pr_tree);
    };
};

int run_main(int argc, char**argv)
{
    bpp::BppApplication natural_extension(argc, argv, "natural-extension");
    std::map<string, string> params = natural_extension.getParams();

    unique_ptr<bpp::Alphabet> alphabet(bpp::SequenceApplicationTools::getAlphabet(params, "", false));

    // Sites
    unique_ptr<bpp::VectorSiteContainer> sites(bpp::SequenceApplicationTools::getSiteContainer(alphabet.get(), params));
    sites.reset(bpp::SequenceApplicationTools::getSitesToAnalyse(*sites, params, "", true, false));
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    // Model, rates
    unique_ptr<bpp::SubstitutionModel> model(bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet.get(), sites.get(), params));
    unique_ptr<bpp::DiscreteDistribution> rate_dist(bpp::PhylogeneticsApplicationTools::getRateDistribution(params));

    // Random seed
    int seed = bpp::ApplicationTools::getIntParameter("natural_extension.seed", params, -1, "", false);
    if(seed != -1)
        RNG.seed(seed);
    // Tree
    int burnin = bpp::ApplicationTools::getIntParameter("natural_extension.burnin", params, 0, "", true);
    std::string nexus_tree_path = bpp::ApplicationTools::getAFilePath("natural_extension.trees_nexus", params, true);
    auto trees = read_nexus_trees(nexus_tree_path, burnin);

    // Parameter of exponential prior
    double exp_mean = bpp::ApplicationTools::getDoubleParameter("natural_extension.exp_mean", params, 10.0, "", false);

    // Output
    std::string output_path = bpp::ApplicationTools::getAFilePath("natural_extension.output_path", params, true, false);
    ofstream output_fp(output_path);

    Exponential_branch_prior forest_prior{exp_mean};

    cerr << "Read " << trees.size() << " trees." << endl;

    size_t i = 0;
    output_fp << "index,forest_length,likelihood,prior,posterior,merged_leaves" << endl;

    for(unique_ptr<bpp::Tree>& t : trees) {
        Forest c = random_split(*t, RNG);
        double fl = forest_likelihood(c, *sites, model.get(), rate_dist.get());
        double pr = forest_prior(c);
        double flen = forest_length(c);
        unsigned int merged = count_merged_leaves(c);

        // No longer invariant - random split
        //assert(std::all_of(begin(c), end(c), [](const TreeTemplate<Node>& t) {
                    //return t.getRootNode()->getNumberOfSons() == 2 &&
                        //t.getRootNode()->getSon(0)->isLeaf() &&
                        //t.getRootNode()->getSon(1)->isLeaf();
        //}));
        output_fp << i++ << "," << flen << "," << fl << "," << pr << "," << fl + pr << "," << merged << endl;
    }

    return 0;
}

int main(int argc, char** argv)
{
    try {
        run_main(argc, argv);
    } catch(bpp::Exception& e) {
        cerr << e.what() << endl;
        return 1;
    }
}
