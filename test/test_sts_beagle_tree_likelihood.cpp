#include "gtest/gtest.h"

#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/HKY85.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>

#include "util.h"
#include "beagle_tree_likelihood.h"
#include "likelihood_vector.h"
#include "online_util.h"

#include <libhmsbeagle/beagle.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>

namespace sts { namespace test { namespace beagle_tree_likelihood {

const bpp::DNA dna;
constexpr double TOLERANCE = 1e-5;

std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tree_of_path(const std::string& newick_path)
{
    bpp::Newick newick_io;
    std::unique_ptr<bpp::Tree> tree(newick_io.read(newick_path));
    std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tt(new bpp::TreeTemplate<bpp::Node>(*tree));
    return tt;
}

std::unique_ptr<bpp::SiteContainer> alignment_of_fasta_path(const std::string& fasta_path, const bpp::Alphabet& alphabet)
{
    std::ifstream aln_stream(fasta_path);
    return std::unique_ptr<bpp::SiteContainer>(sts::util::read_alignment(aln_stream, &alphabet));
}

void test_known_tree(std::string fasta_path,
                     std::string newick_path,
                     bpp::SubstitutionModel& model,
                     bpp::DiscreteDistribution& rate_dist)
{
    using bpp::Node;
    using bpp::Tree;
    using bpp::TreeTemplate;
    using namespace std;

    unique_ptr<bpp::TreeTemplate<Node>> tt = tree_of_path(newick_path);
    unique_ptr<bpp::SiteContainer> aln = alignment_of_fasta_path(fasta_path, dna);

    std::vector<int> node_ids = tt->getNodesId();
    std::sort(node_ids.begin(), node_ids.end());
    for(size_t i = 1; i < node_ids.size(); i++) {
        ASSERT_EQ(node_ids[i], node_ids[i-1] + 1);
    }

    // BEAGLE
    sts::online::BeagleTreeLikelihood beagle_calculator(*aln, model, rate_dist);
    beagle_calculator.initialize(model, rate_dist, *tt);
    beagle_calculator.initialize(model, rate_dist, *tt);
    const double beagle_ll = beagle_calculator.calculateLogLikelihood();
    const size_t llCalls = beagle_calculator.numberOfBeagleUpdateTransitionsCalls();

    // This should not cause any more beagle operations to be executed.
    beagle_calculator.calculateLogLikelihood();
    ASSERT_EQ(llCalls, beagle_calculator.numberOfBeagleUpdateTransitionsCalls());

    // Make dirty
    beagle_calculator.invalidateAll();
    const double beagle_ll_cached = beagle_calculator.calculateLogLikelihood();
    const size_t llCalls2 = beagle_calculator.numberOfBeagleUpdateTransitionsCalls();
    ASSERT_NEAR(beagle_ll, beagle_ll_cached, TOLERANCE);
    ASSERT_EQ(llCalls * 2, llCalls2);

    // Bio++
    bpp::DRHomogeneousTreeLikelihood like(*tt, &model, &rate_dist, false, false);
    like.setData(*aln);
    like.initialize();
    const double bpp_ll = -like.getValue();

    ASSERT_NEAR(beagle_ll, bpp_ll, TOLERANCE);

    const int b = beagle_calculator.getDistalBuffer(tt->getRootNode());
    ASSERT_NEAR(beagle_calculator.logLikelihood(b), beagle_ll, TOLERANCE);
}

TEST(STSBeagleTreeLikelihood, ThirtyJCConstant)
{
    bpp::JCnuc model(&dna);
    bpp::ConstantRateDistribution rates;
    test_known_tree("data/thirty.ma", "data/thirty.tree", model, rates);
}

TEST(STSBeagleTreeLikelihood, ThirtyJCGamma4)
{
    bpp::JCnuc model(&dna);
    bpp::GammaDiscreteRateDistribution rates(4, 0.234);
    test_known_tree("data/thirty.ma", "data/thirty.tree", model, rates);
}

TEST(STSBeagleTreeLikelihood, ThirtyHKYGamma4)
{
    bpp::HKY85 model(&dna, 2.0, 0.4, 0.2, 0.15, 0.25);
    bpp::GammaDiscreteRateDistribution rates(4, 0.234);
    test_known_tree("data/thirty.ma", "data/thirty.tree", model, rates);
}

/// Test that log-likelihood of each mid-edge vector is identical to root.
void test_mid_edge_likelihood_vectors(const std::string& tree_path, const std::string& fasta_path,
                                      const bpp::SubstitutionModel& model, const bpp::DiscreteDistribution& rates)
{
    using namespace bpp;
    using namespace sts::online;
    using std::vector;
    using std::unique_ptr;

    std::unique_ptr<TreeTemplate<Node>> tree = tree_of_path(tree_path);
    std::unique_ptr<SiteContainer> aln = alignment_of_fasta_path(fasta_path, dna);

    sts::online::BeagleTreeLikelihood beagleCalculator(*aln, model, rates);
    beagleCalculator.initialize(model, rates, *tree);

    const vector<BeagleTreeLikelihood::NodePartials> nps = beagleCalculator.getMidEdgePartials();

    const double rootLogLike = beagleCalculator.calculateLogLikelihood();
    for(const BeagleTreeLikelihood::NodePartials& np : nps) {
        double midEdgeLogLike = beagleCalculator.logLikelihood(np.second);
        ASSERT_NEAR(midEdgeLogLike, rootLogLike, TOLERANCE);
    }
}

/// Prune a leaf, verify that logDot log-likelihood is same as full tree
void test_mid_edge_attachment(const std::string& tree_path, const std::string& fasta_path,
                  const bpp::SubstitutionModel& model, const bpp::DiscreteDistribution& rates)
{
    using namespace bpp;
    using namespace sts::online;
    using std::vector;
    using std::unique_ptr;

    std::unique_ptr<TreeTemplate<Node>> tree = tree_of_path(tree_path);
    std::unique_ptr<SiteContainer> aln = alignment_of_fasta_path(fasta_path, dna);

    Node* leaf = tree->getLeaves()[0];
    const std::string leafName = leaf->getName();
    leaf->setDistanceToFather(0.0);

    // Middle of edge
    const double d = leaf->getFather()->getDistanceToFather() + siblings(leaf)[0]->getDistanceToFather();
    leaf->getFather()->setDistanceToFather(d / 2.0);
    siblings(leaf)[0]->setDistanceToFather(d / 2.0);

    sts::online::BeagleTreeLikelihood beagleCalculator(*aln, model, rates);
    beagleCalculator.initialize(model, rates, *tree);
    const double fullLogLikelihood = beagleCalculator.calculateLogLikelihood();

    bpp::Node* sibling = siblings(leaf)[0];
    TreeTemplateTools::dropLeaf(*tree, leafName);

    const int leafBuffer = beagleCalculator.getLeafBuffer(leafName);
    beagleCalculator.initialize(model, rates, *tree);

    const vector<BeagleTreeLikelihood::NodePartials> nps = beagleCalculator.getMidEdgePartials();

    for(const BeagleTreeLikelihood::NodePartials& np : nps) {
        if(np.first == sibling) {
            const double midEdgeLike = beagleCalculator.logDot(leafBuffer, np.second);
            ASSERT_NEAR(midEdgeLike, fullLogLikelihood, TOLERANCE);
            return;
        }
    }
    FAIL() << "Sibling not found.";
}

void test_attachment_likelihood(const std::string& tree_path, const std::string& fasta_path,
                                const bpp::SubstitutionModel& model, const bpp::DiscreteDistribution& rates)
{
    using namespace bpp;
    using namespace sts::online;
    using std::string;
    using std::unique_ptr;
    using std::vector;

    unique_ptr<TreeTemplate<Node>> tree = tree_of_path(tree_path);
    unique_ptr<SiteContainer> aln = alignment_of_fasta_path(fasta_path, dna);

    sts::online::BeagleTreeLikelihood fullCalculator(*aln, model, rates);
    fullCalculator.initialize(model, rates, *tree);
    const double fullLogLikelihood = fullCalculator.calculateLogLikelihood();

    for(const string& leafName : tree->getLeavesNames()) {
        bpp::TreeTemplate<Node> tmpTree(*tree);
        Node* n = tmpTree.getNode(leafName);
        tmpTree.newOutGroup(n);

        const double pendant = n->getDistanceToFather() + sibling(n)->getDistanceToFather();
        // Sibling becomes the root.
        // Insert on first child of the sibling.
        // When sibling becomes root, this edge has all of the length.
        const bpp::Node* sib = sibling(n);
        ASSERT_EQ(static_cast<size_t>(2), sib->getNumberOfSons()) << "Sibling must be bifurcating (dropped " << leafName << ")";
        const bpp::Node* insertEdge = sib->getSon(0);

        const double origInsertLength = insertEdge->getDistanceToFather();
        const double distal = sib->getSon(1)->getDistanceToFather();

        TreeTemplateTools::dropLeaf(tmpTree, leafName);
        fullCalculator.initialize(model, rates, tmpTree);

        ASSERT_NEAR(insertEdge->getDistanceToFather(),
                    distal + origInsertLength,
                    TOLERANCE);

        const double attLike = fullCalculator.calculateAttachmentLikelihood(leafName, insertEdge, distal, {pendant})[0];
        EXPECT_NEAR(fullLogLikelihood, attLike, TOLERANCE) << "removing " << leafName;
    }
}

TEST(STSBeagleTreeLikelihoodMidEdgeThirty, JukesCantorConstant)
{
    bpp::JCnuc model(&dna);
    bpp::ConstantRateDistribution rates;
    test_mid_edge_likelihood_vectors("data/thirty.tree", "data/thirty.ma", model, rates);
    test_attachment_likelihood("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
}

TEST(STSBeagleTreeLikelihoodMidEdgeThirty, JukesCantorGamma6)
{
    bpp::JCnuc model(&dna);
    bpp::GammaDiscreteRateDistribution rates(6, 0.234);
    test_mid_edge_likelihood_vectors("data/thirty.tree", "data/thirty.ma", model, rates);
    test_mid_edge_attachment("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST(STSBeagleTreeLikelihoodMidEdgeThirty, JukesCantorGamma2)
{
    bpp::JCnuc model(&dna);
    bpp::GammaDiscreteRateDistribution rates(2, 0.234);
    test_mid_edge_likelihood_vectors("data/thirty.tree", "data/thirty.ma", model, rates);
    test_mid_edge_attachment("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST(STSBeagleTreeLikelihoodMidEdgeThirty, HKY85Constant)
{
    bpp::HKY85 model(&dna, 2.0, 0.25, 0.25, 0.3, 0.3);
    bpp::ConstantRateDistribution rates;
    test_mid_edge_likelihood_vectors("data/thirty.tree", "data/thirty.ma", model, rates);
    test_mid_edge_attachment("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST(STSBeagleTreeLikelihoodMidEdgeThirty, HKY85Gamma6)
{
    bpp::HKY85 model(&dna, 2.0, 0.4, 0.2, 0.15, 0.25);
    bpp::GammaDiscreteRateDistribution rates(6, 0.234);
    test_mid_edge_likelihood_vectors("data/thirty.tree", "data/thirty.ma", model, rates);
    test_mid_edge_attachment("data/thirty.tree", "data/thirty.ma", model, rates);
    test_attachment_likelihood("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST(STSBeagleTreeLikelihoodMidEdge5taxon, HKY85Gamma6)
{
    bpp::HKY85 model(&dna, 2.0, 0.4, 0.2, 0.15, 0.25);
    bpp::GammaDiscreteRateDistribution rates(6, 0.234);
    test_mid_edge_likelihood_vectors("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
    test_mid_edge_attachment("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
    test_attachment_likelihood("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
}

TEST(STSBeagleTreeLikelihoodMidEdge5taxon, JukesCantorConstant)
{
    bpp::JCnuc model(&dna);
    bpp::ConstantRateDistribution rates;
    test_mid_edge_likelihood_vectors("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
    test_mid_edge_attachment("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
    test_attachment_likelihood("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
}


}}} // namespaces
