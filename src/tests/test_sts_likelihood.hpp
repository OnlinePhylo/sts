/// Tests for STS likelihood calculation
#ifndef STS_TEST_LIKELIHOOD_HPP
#define STS_TEST_LIKELIHOOD_HPP

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
    std::string slurp(const std::string file_name)
    {
        std::ifstream s(file_name);
        std::string str((std::istreambuf_iterator<char>(s)),
                         std::istreambuf_iterator<char>());
        return str;
    }

    TEST_CASE("sts/likelihood/known_tree", "Test calculating the likelihood of a known tree")
    {
        const bpp::DNA dna;
        std::ifstream aln_stream("../data/bppsim/JC69/JC69.fasta");
        //std::ifstream nw_stream("../data/bppsim/JC69/JC69.dnd");
        std::string nwk_string = slurp("../data/bppsim/JC69/JC69.dnd");
        auto aln = std::shared_ptr<bpp::SiteContainer>(sts::util::read_alignment(aln_stream, &dna));
        REQUIRE(aln->getNumberOfSequences() == 79);

        auto model = std::shared_ptr<bpp::SubstitutionModel>(new bpp::JCnuc(&dna));
        auto calc = std::make_shared<sts::likelihood::online_calculator>();
        calc->initialize(aln, model);
        std::unordered_map<sts::particle::node,std::string> names;
        auto root = sts::particle::phylo_particle::of_newick_string(calc, nwk_string, names);
        // Register
        sts::util::register_nodes(*calc, root->node, names);
        std::unordered_set<sts::particle::node> visited;
        double ll = calc->calculate_ll(root->node, visited);
        REQUIRE(-11745.0178177233 == ll);
    }
}
}
}

#endif // STS_TEST_LIKELIHOOD_HPP
