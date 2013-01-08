#ifndef TEST_STS_GMERGE_HPP
#define TEST_STS_GMERGE_HPP

#include "gm_tree.h"

#include "catch.hpp"
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
    GM_tree result = GM_tree::of_newick_path(file_path);
    // FALSE!
    REQUIRE(result.find_k_distance_merges(2).size() == 10000);
}


}
}
}

#endif // TEST_STS_GMERGE_HPP
