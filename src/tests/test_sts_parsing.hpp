#ifndef TEST_STS_PARSING_HPP
#define TEST_STS_PARSING_HPP

#include "sts/particle.hpp"
#include "sts/likelihood/online_calculator.hpp"
#include "sts/util.hpp"

#include <Bpp/Phyl/TreeTemplateTools.h>

#include "catch.hpp"


namespace sts
{
namespace test
{
namespace parsing
{

using namespace sts::particle;
using sts::likelihood::Online_calculator;

std::shared_ptr<Online_calculator> null_calculator;

TEST_CASE("phylofunc/newick_parsing/one_leaf", "test parsing a newick tree with one leaf")
{
    std::string tree = "A;";
    std::unordered_map<Node_ptr, std::string> names;
    sts::particle::Particle p = State::of_newick_string(null_calculator, tree, names);
    REQUIRE(p->node->is_leaf());
}

TEST_CASE("phylofunc/newick_parsing/two_leaf", "test parsing a newick tree with two leaves")
{
    std::string tree = "(A:2,B:3);";
    std::unordered_map<Node_ptr, std::string> names;
    sts::particle::Particle p = State::of_newick_string(null_calculator, tree, names);
    REQUIRE(!p->node->is_leaf());
    REQUIRE(p->node->child1->length == 2);
    REQUIRE(p->node->child1->node->is_leaf());
    REQUIRE(p->node->child2->length == 3);
    REQUIRE(p->node->child2->node->is_leaf());

    std::set<std::shared_ptr<Node>> node_set;
    node_set.insert(p->node->child1->node);
    node_set.insert(p->node->child2->node);
    unsigned int found = 0;
    for(std::shared_ptr<State> cur = p->predecessor; cur; cur = cur->predecessor) {
        REQUIRE(node_set.count(cur->node) > 0);
        ++found;
    }
    REQUIRE(found == node_set.size());
}

TEST_CASE("phylofunc/newick_parsing/three_leaf", "test parsing a newick tree with three leaves")
{
    std::string tree = "((A:2,B:3):4,C:6);";
    std::unordered_map<Node_ptr, std::string> names;
    sts::particle::Particle p = State::of_newick_string(null_calculator, tree, names);
    REQUIRE(!p->node->is_leaf());
    REQUIRE(p->node->child1->length == 4);
    REQUIRE(!p->node->child1->node->is_leaf());
    REQUIRE(p->node->child1->node->child1->length == 2);
    REQUIRE(p->node->child1->node->child1->node->is_leaf());
    REQUIRE(p->node->child1->node->child2->length == 3);
    REQUIRE(p->node->child1->node->child2->node->is_leaf());
    REQUIRE(p->node->child2->length == 6);
    REQUIRE(p->node->child2->node->is_leaf());

    std::set<std::shared_ptr<Node>> node_set;
    node_set.insert(p->node->child1->node);
    node_set.insert(p->node->child1->node->child1->node);
    node_set.insert(p->node->child1->node->child2->node);
    node_set.insert(p->node->child2->node);
    unsigned int found = 0;
    for(std::shared_ptr<State> cur = p->predecessor; cur; cur = cur->predecessor) {
        REQUIRE(node_set.count(cur->node) > 0);
        ++found;
    }
    REQUIRE(found == node_set.size());
}

TEST_CASE("phylofunc/newick_parsing/four_leaf", "test parsing a newick tree with four leaves")
{
    std::string tree = "((A:2,B:3):4,(C:6,D:7):9);";
    std::unordered_map<Node_ptr, std::string> names;
    sts::particle::Particle p = State::of_newick_string(null_calculator, tree, names);
    REQUIRE(!p->node->is_leaf());
    REQUIRE(p->node->child1->length == 4);
    REQUIRE(!p->node->child1->node->is_leaf());
    REQUIRE(p->node->child1->node->child1->length == 2);
    REQUIRE(p->node->child1->node->child1->node->is_leaf());
    REQUIRE(p->node->child1->node->child2->length == 3);
    REQUIRE(p->node->child1->node->child2->node->is_leaf());
    REQUIRE(p->node->child2->length == 9);
    REQUIRE(!p->node->child2->node->is_leaf());
    REQUIRE(p->node->child2->node->child1->length == 6);
    REQUIRE(p->node->child2->node->child1->node->is_leaf());
    REQUIRE(p->node->child2->node->child2->length == 7);
    REQUIRE(p->node->child2->node->child2->node->is_leaf());

    std::set<std::shared_ptr<Node>> node_set;
    node_set.insert(p->node->child1->node);
    node_set.insert(p->node->child1->node->child1->node);
    node_set.insert(p->node->child1->node->child2->node);
    node_set.insert(p->node->child2->node);
    node_set.insert(p->node->child2->node->child1->node);
    node_set.insert(p->node->child2->node->child2->node);
    int found = 0;
    for(std::shared_ptr<State> cur = p->predecessor; cur; cur = cur->predecessor) {
        REQUIRE(node_set.count(cur->node) > 0);
        ++found;
    }
    REQUIRE(found == node_set.size());
}

static std::string
roundtrip(std::string &tree)
{
    std::unordered_map<Node_ptr, std::string> names;
    sts::particle::Particle p = State::of_newick_string(null_calculator, tree, names);
    std::ostringstream ostream;
    sts::util::write_tree(ostream, p->node, names);
    return ostream.str();
}

TEST_CASE("phylofunc/newick_parsing/round_trip/1", "test a parsed newick tree can round-trip")
{
    std::string tree = "((A:2,B:3):4,(C:6,D:7):9);\n";
    std::vector<std::string> names {"A", "B", "", "C", "D", ""};
//    REQUIRE(roundtrip(tree, names) == tree);
    REQUIRE(roundtrip(tree) == tree);
}

TEST_CASE("phylofunc/newick_parsing/round_trip/2", "test a parsed newick tree can round-trip")
{
    std::string tree = "((A:2.5,((B:3.25,C:4.125):5,D:6):7.5):7.75,((E:8,F:9):10,G:11):9);\n";
//   std::vector<std::string> names {"A", "B", "C", "", "D", "", "", "E", "F", "", "G", ""};
    REQUIRE(roundtrip(tree) == tree);
}

}
}
}

#endif // TEST_STS_PARSING_HPP
