#ifndef STS_LIKELIHOOD_DETAIL_ONLINE_CALCULATOR_HPP
#define STS_LIKELIHOOD_DETAIL_ONLINE_CALCULATOR_HPP

#include <memory>
#include <string>
#include <vector>
#include <stack>
#include <unordered_map>
#include <unordered_set>

#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include "libhmsbeagle/beagle.h"

namespace sts
{
namespace particle
{
class Node;
}

namespace likelihood
{

// XXX should we come up with a convention for method order? I'd be happy with either approximate dependency, split into
// public and private or not.
class Online_calculator
{
public:
    Online_calculator() : verify_cached_ll(false), instance(-1), next_id(0) {};
    ~Online_calculator();
    void initialize(std::shared_ptr<bpp::SiteContainer>, std::shared_ptr<bpp::SubstitutionModel>);
    int get_id();
    void free_id(int id);
    double calculate_ll(std::shared_ptr<sts::particle::Node> node, std::unordered_set<std::shared_ptr<sts::particle::Node>>& visited);
    void invalidate(std::shared_ptr< sts::particle::Node > n);
    void register_node(std::shared_ptr< sts::particle::Node > n);
    void register_leaf(std::shared_ptr< sts::particle::Node > n, const std::string taxon);
    void unregister_node(const sts::particle::Node* n);

    void set_weights(std::vector<double> weights);

    /// If \c true, unconditionally recomputes the log likelihood of \c node.
    /// If the log-likelihood of \c node has previously been calculated, verifies that the new value matches the cached
    /// value.
    bool verify_cached_ll;
private:
    BeagleInstanceDetails instance_details;
    std::shared_ptr<bpp::SiteContainer> sites;
    std::shared_ptr<bpp::SubstitutionModel> model;
    int num_buffers;
    int instance;
    int next_id;
    std::stack<int> free_ids;
    std::unordered_map<const sts::particle::Node*, double> node_ll_map; // caches the root ll at each node
    std::unordered_map<const sts::particle::Node*, int> node_buffer_map; // maps nodes to a beagle buffer ID
    std::unordered_map<std::string, int> taxon_buffer_map; // maps taxon names to beagle buffer id.

    int create_beagle_instance();
    void grow();
    void set_eigen_and_rates_and_weights(int instance);
    void set_eigen_and_rates_and_weights(int, const bpp::SubstitutionModel&);
    int get_buffer(std::shared_ptr< sts::particle::Node > n);
};

}
}

#endif // STS_LIKELIHOOD_DETAIL_ONLINE_CALCULATOR_HPP
