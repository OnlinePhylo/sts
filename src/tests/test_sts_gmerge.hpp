#ifndef TEST_STS_GMERGE_HPP
#define TEST_STS_GMERGE_HPP

#include "gm_tree.h"
#include "util.h"

#include "catch.hpp"
#include <iostream>
#include <string>
#include <unordered_map>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/TreeTemplateTools.h>


namespace sts
{
namespace test
{
namespace gmerge
{

std::unordered_map<std::string, sts::particle::Node_ptr> name_map(const bpp::TreeTemplate<bpp::Node>& t)
{
    std::unordered_map<std::string, sts::particle::Node_ptr> result;
    for(const bpp::Node* node : t.getLeaves()) {
        result[node->getName()] = std::make_shared<sts::particle::Node>(nullptr);
    }
    return result;
}

typedef bpp::TreeTemplate<bpp::Node> BppTree;
const std::string test_topology = "(D,(C,(B,A)));";
const std::unique_ptr<BppTree> tree(bpp::TreeTemplateTools::parenthesisToTree(test_topology));
const std::unordered_map<std::string,sts::particle::Node_ptr> names_nodes = name_map(*tree);


TEST_CASE("sts/guidedmerge/test_gmtree_copy", "")
{
    using sts::particle::Node_ptr;
    using namespace sts::guidedmerge;
    using namespace std;
    const GM_tree result = GM_tree::of_newick_string(test_topology, names_nodes);
    auto nodes_names = sts::util::unordered_invert(names_nodes);
    std::unique_ptr<BppTree> t1(result.to_treetemplate(nodes_names));
    bpp::TreeTemplateTools::orderTree(*t1->getRootNode(), true, true);
    GM_tree copy_ctor(result);

    // Check copied tree has same topology
    std::unique_ptr<BppTree> t2(copy_ctor.to_treetemplate(nodes_names));
    bpp::TreeTemplateTools::orderTree(*t2->getRootNode(), true, true);
    REQUIRE(bpp::TreeTemplateTools::haveSameOrderedTopology(*t1->getRootNode(), *t2->getRootNode()));

    // Assignment
    GM_tree assign;
    assign = result;

    // Check original unchanged
    std::unique_ptr<BppTree> t3(result.to_treetemplate(nodes_names));
    bpp::TreeTemplateTools::orderTree(*t3->getRootNode(), true, true);
    REQUIRE(bpp::TreeTemplateTools::haveSameOrderedTopology(*t1->getRootNode(), *t3->getRootNode()));

    // Check same topology
    std::unique_ptr<BppTree> t4(assign.to_treetemplate(nodes_names));
    bpp::TreeTemplateTools::orderTree(*t4->getRootNode(), true, true);
    REQUIRE(bpp::TreeTemplateTools::haveSameOrderedTopology(*t1->getRootNode(), *t4->getRootNode()));
}

TEST_CASE("sts/guidedmerge/test_merge", "")
{
    using sts::particle::Node_ptr;
    using namespace sts::guidedmerge;
    using namespace std;
    std::unordered_map<string,Node_ptr> names_nodes = name_map(*tree);
    GM_tree result = GM_tree::of_treetemplate(*tree, names_nodes);

    std::unordered_map<Node_ptr,string> nodes_names = sts::util::unordered_invert(names_nodes);
    REQUIRE(names_nodes.size() == nodes_names.size());

    REQUIRE(4 == result.get_leaf_count());

    // Check path between "A" and "C"
    auto n1 = names_nodes.at("A"), n2 = names_nodes.at("C");
    REQUIRE(result.rf_distance(n1, n2) == 2);
    REQUIRE(result.path_exists(names_nodes.at("A"), names_nodes.at("B")));

    Node_ptr p = make_shared<sts::particle::Node>(nullptr);
    // Merge "A" and "C"
    result.merge(n1, n2, p);
    // Should reduce leaf count by 1
    REQUIRE(result.get_leaf_count() == 3);

    // Now merge the newly created node with "D"
    REQUIRE(result.path_exists(p, names_nodes.at("D")));
    Node_ptr p2 = make_shared<sts::particle::Node>(nullptr);

    result.merge(names_nodes["D"], p, p2);
    REQUIRE(result.get_leaf_count() == 2);
}

TEST_CASE("sts/guidedmerge/parsing_success_thirty", "")
{
    using namespace sts::guidedmerge;
    using namespace std;
    const std::string file_path = "data/thirty.tree";
    bpp::Newick nwk;
    unique_ptr<BppTree> tree(nwk.read(file_path));
    std::unordered_map<std::string,sts::particle::Node_ptr> names = name_map(*tree);

    const GM_tree result = GM_tree::of_treetemplate(*tree, names);

    // TODO: Expand
    // Check that merge size counts match python/sts.py for the same tree
    // NB: sts.py includes pendant edges in `k`, while this code does not.
    REQUIRE(result.find_k_distance_merges(0).size() == 9);
    REQUIRE(result.find_k_distance_merges(1).size() == 15);
    REQUIRE(result.find_k_distance_merges(2).size() == 16);
}

TEST_CASE("sts/guidedmerge/gmtree_roundtrip_thirty", "")
{
    using namespace sts::guidedmerge;
    using namespace std;
    const std::string file_path = "data/thirty.tree";
    bpp::Newick nwk;
    unique_ptr<BppTree> tree(nwk.read(file_path));
    std::unordered_map<std::string,sts::particle::Node_ptr> names = name_map(*tree);

    const GM_tree result = GM_tree::of_treetemplate(*tree, names);
    std::string n = result.to_newick_string(sts::util::unordered_invert(names));

    // TODO: Expand
    REQUIRE(n.length() > 0);
}

}
}
}

#endif // TEST_STS_GMERGE_HPP
