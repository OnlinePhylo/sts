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
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;
using bpp::TreeTemplate;
using bpp::Node;

typedef vector<TreeTemplate<Node>> forest;

/// Decompose a 4-taxon tree into a pair of cherries.
forest cherries(TreeTemplate<Node> tree)
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
            forest result;
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

double forest_likelihood(const forest& f, const bpp::SiteContainer& data, bpp::SubstitutionModel* model, bpp::DiscreteDistribution* rates)
{
    double result = 0;

    for (auto &tree : f) {
        bpp::RHomogeneousTreeLikelihood like(tree, data, model, rates, true, false);
        like.initialize();
        result += like.getLogLikelihood();
    }
    return result;
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

struct exponential_branch_prior
{
    double mu;

    double operator()(const forest& f)
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

    // Tree
    int burnin = bpp::ApplicationTools::getIntParameter("natural_extension.burnin", params, 0, "", true);
    std::string nexus_tree_path = bpp::ApplicationTools::getAFilePath("natural_extension.trees_nexus", params, true);
    auto trees = read_nexus_trees(nexus_tree_path, burnin);

    // Parameter of exponential prior
    double exp_mean = bpp::ApplicationTools::getDoubleParameter("natural_extension.exp_mean", params, 10.0, "", false);

    // Output
    std::string output_path = bpp::ApplicationTools::getAFilePath("natural_extension.output_path", params, true, false);
    ofstream output_fp(output_path);

    auto forest_prior = exponential_branch_prior{exp_mean};

    cerr << "Read " << trees.size() << " trees." << endl;

    size_t i = 0;
    output_fp << "index,likelihood,prior,posterior" << endl;
    for(unique_ptr<bpp::Tree>& t : trees) {
        forest c = cherries(*t);
        double fl = forest_likelihood(c, *sites, model.get(), rate_dist.get());
        double pr = forest_prior(c);
        assert(std::all_of(begin(c), end(c), [](const TreeTemplate<Node>& t) {
                    return t.getRootNode()->getNumberOfSons() == 2 &&
                        t.getRootNode()->getSon(0)->isLeaf() &&
                        t.getRootNode()->getSon(1)->isLeaf();
        }));
        output_fp << i++ << "," << fl << "," << pr << "," << fl + pr << endl;
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
