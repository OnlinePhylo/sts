/// Tests for STS likelihood calculation

#include "gtest/gtest.h"

#include <cmath>
#include <fstream>
#include <memory>
#include <streambuf>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>

#include "online_calculator.h"
#include "particle.h"
#include "state.h"
#include "util.h"

namespace sts
{
namespace test
{
namespace likelihood
{

constexpr double TOL = 1e-2;

// file to string
std::string slurp(const std::string file_name)
{
    std::ifstream s(file_name);
    return std::string((std::istreambuf_iterator<char>(s)),
                       std::istreambuf_iterator<char>());
}

void test_known_tree_jc69(std::string fasta_path, std::string newick_path, double log_likelihood, bool compress)
{
    const bpp::DNA dna;
    std::ifstream aln_stream(fasta_path);
    std::string nwk_string = slurp(newick_path);
    auto aln = std::shared_ptr<bpp::SiteContainer>(sts::util::read_alignment(aln_stream, &dna));
    auto compressed_aln = std::shared_ptr<bpp::SiteContainer>(util::unique_sites(*aln));
    auto weights = util::compressed_site_weights(*aln, *compressed_aln);

    ASSERT_LE(compressed_aln->getNumberOfSites(), aln->getNumberOfSites());

    auto model = std::shared_ptr<bpp::SubstitutionModel>(new bpp::JCnuc(&dna));
    auto calc = std::make_shared<sts::likelihood::Online_calculator>();
    calc->initialize(compress ? compressed_aln : aln, model);
    if(compress)
        calc->set_weights(weights);
    std::unordered_map<sts::particle::Node_ptr, std::string> names;
    auto root = sts::particle::State::of_newick_string(calc, nwk_string, names);
    // Register
    sts::util::register_nodes(*calc, root->node, names);
    std::unordered_set<sts::particle::Node_ptr> visited;
    double ll = calc->calculate_ll(root->node, visited);

    ASSERT_NEAR(log_likelihood, ll, TOL);
}


TEST(STSLikelihoodKnownTree, Compress)
{
    test_known_tree_jc69("data/bppsim/JC69/JC69.fasta", "data/bppsim/JC69/JC69.dnd",
                         -11745.0178177233, true);
}

TEST(STSLikelihoodKnownTree, NoCompress)
{
    test_known_tree_jc69("data/bppsim/JC69/JC69.fasta", "data/bppsim/JC69/JC69.dnd",
                         -11745.0178177233, false);
}

TEST(STSLikelihoodKnownTree, ThirtyCompress)
{
    test_known_tree_jc69("data/thirty.ma", "data/thirty.tree", -18464.9, true);
}

}
}
}
