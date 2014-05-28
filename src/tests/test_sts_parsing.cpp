#include "edge.h"
#include "online_calculator.h"
#include "node.h"
#include "particle.h"
#include "state.h"
#include "util.h"

#include <Bpp/Phyl/TreeTemplateTools.h>

#include "gtest/gtest.h"


namespace sts
{
namespace test
{
namespace parsing
{

using namespace sts::particle;
using sts::likelihood::Online_calculator;

std::shared_ptr<Online_calculator> null_calculator;

TEST(phylofunc_newick_parsing, one_leaf)
{
    std::string tree = "A;";
    std::unordered_map<Node_ptr, std::string> names;
    sts::particle::Particle p = State::of_newick_string(null_calculator, tree, names);
    ASSERT_TRUE(p->node->is_leaf());
}

TEST(phylofunc_newick_parsing, two_leaf)
{
    std::string tree = "(A:2,B:3);";
    std::unordered_map<Node_ptr, std::string> names;
    sts::particle::Particle p = State::of_newick_string(null_calculator, tree, names);
    ASSERT_TRUE(!p->node->is_leaf());
    ASSERT_EQ(p->node->child1->length, 2);
    ASSERT_TRUE(p->node->child1->node->is_leaf());
    ASSERT_EQ(p->node->child2->length, 3);
    ASSERT_TRUE(p->node->child2->node->is_leaf());

    std::set<std::shared_ptr<Node>> node_set;
    node_set.insert(p->node->child1->node);
    node_set.insert(p->node->child2->node);
    unsigned int found = 0;
    for(std::shared_ptr<State> cur = p->predecessor; cur; cur = cur->predecessor) {
        ASSERT_GT(node_set.count(cur->node), 0);
        ++found;
    }
    ASSERT_EQ(found, node_set.size());
}

TEST(phylofunc_newick_parsing, three_leaf)
{
    std::string tree = "((A:2,B:3):4,C:6);";
    std::unordered_map<Node_ptr, std::string> names;
    sts::particle::Particle p = State::of_newick_string(null_calculator, tree, names);
    ASSERT_TRUE(!p->node->is_leaf());
    ASSERT_EQ(p->node->child1->length, 4);
    ASSERT_TRUE(!p->node->child1->node->is_leaf());
    ASSERT_EQ(p->node->child1->node->child1->length, 2);
    ASSERT_TRUE(p->node->child1->node->child1->node->is_leaf());
    ASSERT_EQ(p->node->child1->node->child2->length, 3);
    ASSERT_TRUE(p->node->child1->node->child2->node->is_leaf());
    ASSERT_EQ(p->node->child2->length, 6);
    ASSERT_TRUE(p->node->child2->node->is_leaf());

    std::set<std::shared_ptr<Node>> node_set;
    node_set.insert(p->node->child1->node);
    node_set.insert(p->node->child1->node->child1->node);
    node_set.insert(p->node->child1->node->child2->node);
    node_set.insert(p->node->child2->node);
    unsigned int found = 0;
    for(std::shared_ptr<State> cur = p->predecessor; cur; cur = cur->predecessor) {
        ASSERT_GT(node_set.count(cur->node), 0);
        ++found;
    }
    ASSERT_EQ(found, node_set.size());
}

TEST(phylofunc_newick_parsing, four_leaf)
{
    std::string tree = "((A:2,B:3):4,(C:6,D:7):9);";
    std::unordered_map<Node_ptr, std::string> names;
    sts::particle::Particle p = State::of_newick_string(null_calculator, tree, names);
    ASSERT_FALSE(p->node->is_leaf());
    ASSERT_EQ(p->node->child1->length, 4);
    ASSERT_FALSE(p->node->child1->node->is_leaf());
    ASSERT_EQ(p->node->child1->node->child1->length, 2);
    ASSERT_TRUE(p->node->child1->node->child1->node->is_leaf());
    ASSERT_EQ(p->node->child1->node->child2->length, 3);
    ASSERT_TRUE(p->node->child1->node->child2->node->is_leaf());
    ASSERT_EQ(p->node->child2->length, 9);
    ASSERT_FALSE(p->node->child2->node->is_leaf());
    ASSERT_EQ(p->node->child2->node->child1->length, 6);
    ASSERT_TRUE(p->node->child2->node->child1->node->is_leaf());
    ASSERT_EQ(p->node->child2->node->child2->length, 7);
    ASSERT_TRUE(p->node->child2->node->child2->node->is_leaf());

    std::set<std::shared_ptr<Node>> node_set;
    node_set.insert(p->node->child1->node);
    node_set.insert(p->node->child1->node->child1->node);
    node_set.insert(p->node->child1->node->child2->node);
    node_set.insert(p->node->child2->node);
    node_set.insert(p->node->child2->node->child1->node);
    node_set.insert(p->node->child2->node->child2->node);
    unsigned int found = 0;
    for(std::shared_ptr<State> cur = p->predecessor; cur; cur = cur->predecessor) {
        ASSERT_GT(node_set.count(cur->node), 0);
        ++found;
    }
    ASSERT_EQ(found, node_set.size());
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

TEST(phylofunc_newick_parsing, round_trip_1)
{
    std::string tree = "((A:2,B:3):4,(C:6,D:7):9);\n";
    std::vector<std::string> names {"A", "B", "", "C", "D", ""};
    ASSERT_EQ(roundtrip(tree), tree);
}

TEST(phylofunc_newick_parsing, round_trip_2)
{
    std::string tree = "((A:2.5,((B:3.25,C:4.125):5,D:6):7.5):7.75,((E:8,F:9):10,G:11):9);\n";
//   std::vector<std::string> names {"A", "B", "C", "", "D", "", "", "E", "F", "", "G", ""};
    ASSERT_EQ(roundtrip(tree), tree);
}

}
}
}
