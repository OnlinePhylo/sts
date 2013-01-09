#ifndef TEST_STS_GMERGE_HPP
#define TEST_STS_GMERGE_HPP

#include "gm_tree.h"
#include "util.h"

#include "catch.hpp"
#include <iostream>
#include <string>


namespace sts
{
namespace test
{
namespace gmerge
{


TEST_CASE("sts/guidedmerge/test_merge", "")
{
    using sts::particle::Node_ptr;
    using namespace sts::guidedmerge;
    using namespace std;
    const std::string test_topology = "(D,(C,(B,A)));";
    std::unordered_map<string,Node_ptr> names_nodes;
    GM_tree result = GM_tree::of_newick_string(test_topology, names_nodes);

    std::unordered_map<Node_ptr,string> nodes_names = sts::util::unordered_invert(names_nodes);
    REQUIRE(names_nodes.size() == nodes_names.size());

    auto leaves = result.get_leaves();
    REQUIRE(leaves.size() == 4);

    unordered_map<string,GM_node_ptr> names_gmnodes;
    for(auto l : leaves) {
        REQUIRE(static_cast<bool>(l->node));
        REQUIRE(nodes_names.count(l->node) > 0);
        names_gmnodes[nodes_names[l->node]] = l;
    }

    // Check path between "A" and "C"
    auto n1 = names_gmnodes.at("A"), n2 = names_gmnodes.at("C");
    auto p = result.find_path(n1, n2);
    REQUIRE(p.size() == 4);

    // Merge "A" and "C"
    result.merge(n1, n2);

    // Should reduce leaf count by 1
    auto m_leaves = result.get_leaves();
    // "B" and "D" nodes are be unchanged
    REQUIRE(m_leaves.count(names_gmnodes.at("B")) == 1);
    REQUIRE(m_leaves.count(names_gmnodes.at("D")) == 1);
    REQUIRE(result.get_leaves().size() == 3);
    // "A" and "C" nodes are no longer present
    REQUIRE(m_leaves.count(names_gmnodes.at("A")) == 0);
    REQUIRE(m_leaves.count(names_gmnodes.at("C")) == 0);
}

TEST_CASE("sts/guidedmerge/parsing_success_thirty", "")
{
    using namespace sts::guidedmerge;
    using namespace std;
    std::string file_path = "data/thirty.tree";
    std::unordered_map<std::string,sts::particle::Node_ptr> names;
    const GM_tree result = GM_tree::of_newick_path(file_path, names);

    // TODO: Expand
    // Check that merge size counts match python/sts.py for the same tree
    REQUIRE(result.find_k_distance_merges(2).size() == 9);
    REQUIRE(result.find_k_distance_merges(3).size() == 15);
    REQUIRE(result.find_k_distance_merges(4).size() == 16);
}

TEST_CASE("sts/guidedmerge/gmtree_roundtrip_thirty", "")
{
    using namespace sts::guidedmerge;
    using namespace std;
    std::string file_path = "data/thirty.tree";
    std::unordered_map<std::string,sts::particle::Node_ptr> names;
    GM_tree result = GM_tree::of_newick_path(file_path, names);
    std::string nwk = result.to_newick_string(sts::util::unordered_invert(names));

    // TODO: Expand
    REQUIRE(nwk.length() > 0);
}

}
}
}

#endif // TEST_STS_GMERGE_HPP
