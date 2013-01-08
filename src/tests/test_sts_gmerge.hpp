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

TEST_CASE("sts/guidedmerge/parsing_success", "")
{
    using namespace sts::guidedmerge;
    using namespace std;
    std::string file_path = "data/thirty.tree";
    std::unordered_map<std::string,sts::particle::Node_ptr> names;
    GM_tree result = GM_tree::of_newick_path(file_path, names);

    // TODO: Expand
    // Check that merge size counts match python/sts.py for the same tree
    REQUIRE(result.find_k_distance_merges(2).size() == 9);
    REQUIRE(result.find_k_distance_merges(3).size() == 15);
    REQUIRE(result.find_k_distance_merges(4).size() == 16);
}

TEST_CASE("sts/guidedmerge/gmtree_roundtrip", "")
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
