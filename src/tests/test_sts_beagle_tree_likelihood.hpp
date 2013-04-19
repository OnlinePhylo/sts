#include "catch.hpp"

#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/JCnuc.h>
#include <Bpp/Phyl/Model/HKY85.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>

#include "util.h"
#include "beagle_tree_likelihood.h"
#include "likelihood_vector.h"

#include <libhmsbeagle/beagle.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>

namespace sts { namespace test { namespace beagle_tree_likelihood {

const bpp::DNA dna;

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
        REQUIRE(node_ids[i] == node_ids[i-1] + 1);
    }

    // BEAGLE
    sts::online::BeagleTreeLikelihood beagle_calculator(*aln, model, rate_dist);
    beagle_calculator.initialize(model, rate_dist, *tt);
    beagle_calculator.initialize(model, rate_dist, *tt);
    const double beagle_ll = beagle_calculator.calculateLogLikelihood();

    // Make dirty
    beagle_calculator.invalidate(tt->getNodes()[4]);
    const double beagle_ll_cached = beagle_calculator.calculateLogLikelihood();
    CHECK(beagle_ll == Approx(beagle_ll_cached));

    // Bio++
    bpp::DRHomogeneousTreeLikelihood like(*tt, &model, &rate_dist, false, false);
    like.setData(*aln);
    like.initialize();
    const double bpp_ll = -like.getValue();

    CHECK(beagle_ll == Approx(bpp_ll));
}

TEST_CASE("sts/beagle_tree_likelihood/thirty/JC/constant", "thirty.ma, jukes-cantor, constant rates")
{
    bpp::JCnuc model(&dna);
    bpp::ConstantDistribution rates(1.0);
    test_known_tree("data/thirty.ma", "data/thirty.tree", model, rates);
}

TEST_CASE("sts/beagle_tree_likelihood/thirty/JC/gamma4", "thirty.ma, jukes-cantor, gamma4 rates")
{
    bpp::JCnuc model(&dna);
    bpp::GammaDiscreteDistribution rates(4, 0.234);
    test_known_tree("data/thirty.ma", "data/thirty.tree", model, rates);
}

TEST_CASE("sts/beagle_tree_likelihood/thirty/HKY/gamma4", "thirty.ma, HKY, gamma4 rates")
{
    bpp::HKY85 model(&dna, 2.0, 0.4, 0.2, 0.15, 0.25);
    bpp::GammaDiscreteDistribution rates(4, 0.234);
    test_known_tree("data/thirty.ma", "data/thirty.tree", model, rates);
}

void test_mid_edge_likelihood_vectors(const std::string& tree_path, const std::string& fasta_path,
                                      const bpp::SubstitutionModel& model, const bpp::DiscreteDistribution& rates)
{
    using namespace bpp;
    using namespace sts::online;
    using std::vector;
    using std::unique_ptr;

    std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tree = tree_of_path(tree_path);
    std::unique_ptr<bpp::SiteContainer> aln = alignment_of_fasta_path(fasta_path, dna);

    sts::online::BeagleTreeLikelihood beagleCalculator(*aln, model, rates);
    beagleCalculator.initialize(model, rates, *tree);

    const vector<BeagleTreeLikelihood::NodePartials> nps = beagleCalculator.getMidEdgePartials();

    const double rootLogLike = beagleCalculator.calculateLogLikelihood();
    for(const BeagleTreeLikelihood::NodePartials& np : nps) {
        double midEdgeLogLike = beagleCalculator.logLikelihood(np.second);
        REQUIRE(midEdgeLogLike == Approx(rootLogLike));
    }
}

TEST_CASE("sts/beagle_tree_likelihood/mid_edge/thirty/jukes_cantor/constant", "Test mid-edge partials")
{
    bpp::JCnuc model(&dna);
    bpp::ConstantDistribution rates(1.0);
    test_mid_edge_likelihood_vectors("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST_CASE("sts/beagle_tree_likelihood/mid_edge/thirty/jukes_cantor/gamma6", "Test mid-edge partials")
{
    bpp::JCnuc model(&dna);
    bpp::GammaDiscreteDistribution rates(6, 0.234);
    test_mid_edge_likelihood_vectors("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST_CASE("sts/beagle_tree_likelihood/mid_edge/thirty/jukes_cantor/gamma2", "Test mid-edge partials")
{
    bpp::JCnuc model(&dna);
    bpp::GammaDiscreteDistribution rates(2, 0.234);
    test_mid_edge_likelihood_vectors("data/thirty.tree", "data/thirty.ma", model, rates);
}


TEST_CASE("sts/beagle_tree_likelihood/mid_edge/thirty/hky85/constant", "Test mid-edge partials")
{
    bpp::HKY85 model(&dna, 2.0, 0.25, 0.25, 0.3, 0.3);
    bpp::ConstantDistribution rates(1.0);
    test_mid_edge_likelihood_vectors("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST_CASE("sts/beagle_tree_likelihood/mid_edge/thirty/hky85/gamma6", "Test mid-edge partials")
{
    bpp::HKY85 model(&dna, 2.0, 0.4, 0.2, 0.15, 0.25);
    bpp::GammaDiscreteDistribution rates(6, 0.234);
    test_mid_edge_likelihood_vectors("data/thirty.tree", "data/thirty.ma", model, rates);
}

TEST_CASE("sts/beagle_tree_likelihood/mid_edge/5taxon/hky85/gamma6", "Test mid-edge partials")
{
    bpp::HKY85 model(&dna, 2.0, 0.4, 0.2, 0.15, 0.25);
    bpp::GammaDiscreteDistribution rates(6, 0.234);
    test_mid_edge_likelihood_vectors("data/5taxon/5taxon.tre", "data/5taxon/5taxon.fasta", model, rates);
}


}}} // namespaces
