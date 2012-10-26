/// Tests for STS likelihood calculation
#ifndef STS_TEST_LIKELIHOOD_HPP
#define STS_TEST_LIKELIHOOD_HPP

#include <cmath>
#include <fstream>
#include <memory>
#include <streambuf>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Phyl/Model/JCnuc.h>
#include "sts/likelihood/online_calculator.hpp"
#include "sts/particle.hpp"
#include "sts/util.hpp"
#include "catch.hpp"

namespace sts
{
namespace test
{
namespace likelihood
{

// file to string
std::string slurp(const std::string file_name)
{
    std::ifstream s(file_name);
    std::string str((std::istreambuf_iterator<char>(s)),
                    std::istreambuf_iterator<char>());
    return str;
}

void test_known_tree_jc69(std::string fasta_path, std::string newick_path, double log_likelihood, bool compress)
{
    const double tol = 1e-5;
    const bpp::DNA dna;
    std::ifstream aln_stream(fasta_path);
    std::string nwk_string = slurp(newick_path);
    auto aln = std::shared_ptr<bpp::SiteContainer>(sts::util::read_alignment(aln_stream, &dna));
    auto compressed_aln = std::shared_ptr<bpp::SiteContainer>(util::unique_sites(*aln));
    auto weights = util::compressed_site_weights(*aln, *compressed_aln);

    REQUIRE(compressed_aln->getNumberOfSites() <= aln->getNumberOfSites());

    auto model = std::shared_ptr<bpp::SubstitutionModel>(new bpp::JCnuc(&dna));
    auto calc = std::make_shared<sts::likelihood::online_calculator>();
    calc->initialize(compress ? compressed_aln : aln, model);
    if(compress)
        calc->set_weights(weights);
    std::unordered_map<sts::particle::node, std::string> names;
    auto root = sts::particle::phylo_particle::of_newick_string(calc, nwk_string, names);
    // Register
    sts::util::register_nodes(*calc, root->node, names);
    std::unordered_set<sts::particle::node> visited;
    double ll = calc->calculate_ll(root->node, visited);

    REQUIRE(std::abs(log_likelihood - ll) < 0.1);
}


TEST_CASE("sts/likelihood/known_tree/compress", "Test calculating the likelihood of a known tree with compressed sites")
{
    test_known_tree_jc69("../data/bppsim/JC69/JC69.fasta", "../data/bppsim/JC69/JC69.dnd",
            -11745.0178177233, true);
}

TEST_CASE("sts/likelihood/known_tree/no_compress", "Test calculating the likelihood of a known tree without compressing sites")
{
    test_known_tree_jc69("../data/bppsim/JC69/JC69.fasta", "../data/bppsim/JC69/JC69.dnd",
            -11745.0178177233, false);
}

TEST_CASE("sts/likelihood/known_tree/thirty/compress", "Test calculating the likelihood of thirty.ma")
{
    test_known_tree_jc69("../data/thirty.ma", "../data/thirty.tree", -18464.9, true);
}

}
}
}

#endif // STS_TEST_LIKELIHOOD_HPP
