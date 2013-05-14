/// \file online_util.h
/// \brief Utility functions
/// Utility functions for online inference
#ifndef STS_ONLINE_ONLINE_UTILS_H
#define STS_ONLINE_ONLINE_UTILS_H

#include <forward_list>
#include <stack>
#include <vector>
#include <Bpp/Phyl/TreeTemplate.h>

namespace sts { namespace online {

/// \brief Extract a vector of nodes in postorder from \c root
///
/// \param root Tree root
/// \tparam N Node type
template<typename N>
std::vector<N*> postorder(N* root)
{
    std::stack<N*> to_process;
    std::forward_list<N*> result;
    to_process.push(root);
    while(!to_process.empty()) {
        N* n = to_process.top();
        to_process.pop();
        result.push_front(n);
        for(size_t i = 0; i < n->getNumberOfSons(); i++)
            to_process.push(n->getSon(i));
    }
    return std::vector<N*>(result.begin(), result.end());
}

/// \brief Extract a vector of nodes in preorder from \c root
///
/// \param root Tree root
/// \tparam N Node type
template<typename N>
std::vector<N*> preorder(N* root)
{
    std::vector<N*> result;
    std::stack<N*> to_process;
    to_process.push(root);
    while(!to_process.empty()) {
        N* n = to_process.top();
        to_process.pop();
        result.push_back(n);

        // Add sons in reverse order for left-to-right traversal: FILO
        for(int i = n->getNumberOfSons() - 1; i >= 0; i--)
            to_process.push(n->getSon(i));

    }
    return result;
}

/// \p List siblings of \c node
///
/// \param node Node to examine
/// \tparam N Node type
/// \returns Vector of nodes with same father as \c node. If \c node has no father, an empty vector;
template <typename N>
std::vector<N*> siblings(N* node)
{
    std::vector<N*> result;
    if(!node->hasFather())
        return result;
    N* f = node->getFather();
    for(size_t i = 0; i < f->getNumberOfSons(); i++) {
        N* n = f->getSon(i);
        if(n != node)
            result.push_back(n);
    }
    return result;
}

/// \brief generate a vector of edges which may be placed on / modified
///
/// The online inference algorithm uses a bifurcating tree - the root is placed on an arbitrary edge.
/// To create a uniform proposal density on edges, the edge to the right of the root has branch length 0,
/// and the edge to the left is used for all operations.
/// This function returns a list of all nodes *except* the root and node to the right of the root.
///
/// \param tree Tree
/// \tparam N node type
/// \return Available nodes
template <typename N>
std::vector<N*> onlineAvailableEdges(bpp::TreeTemplate<N>& tree)
{
    N* root = tree.getRootNode();
    assert(root->getNumberOfSons() == 2);
    N* right_of_root = root->getSon(1);
    assert(right_of_root->getDistanceToFather() <= 1e-6);

    std::vector<N*> r = tree.getNodes();
    auto f = [root, right_of_root](const N* n) { return n == root || n == right_of_root; };
    r.erase(std::remove_if(r.begin(), r.end(), f), r.end());

    return r;
}

struct SummaryStatistics
{
    double mean, median, sd, min, max, lower95, upper95;
};
SummaryStatistics summarize(std::vector<double> values);

}}

namespace std {
std::ostream& operator<<(std::ostream& os, const sts::online::SummaryStatistics& value);
}

#endif
