#ifndef STS_LIKELIHOOD_DETAIL_ONLINE_CALCULATOR_FWD_HPP
#define STS_LIKELIHOOD_DETAIL_ONLINE_CALCULATOR_FWD_HPP

#include <boost/shared_ptr.hpp>
#include <memory>
#include <stack>
#include <unordered_map>
#include <unordered_set>

#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include "sts/likelihood/bpp_shim.hpp"

#include "libhmsbeagle/beagle.h"

namespace sts
{
namespace particle
{
class phylo_node;
}

namespace likelihood
{

// XXX should we come up with a convention for method order? I'd be happy with either approximate dependency, split into
// public and private or not.
class online_calculator
{
public:
    online_calculator() : initialized(false), instance(-1), next_id(0)  {};
    ~online_calculator();
    void initialize(boost::shared_ptr<bpp::SiteContainer>, boost::shared_ptr<bpp::SubstitutionModel>);
    int get_id();
    void free_id(int id);
    double calculate_ll(boost::shared_ptr<sts::particle::phylo_node> node, std::unordered_set<boost::shared_ptr<sts::particle::phylo_node>>& visited);
    void register_node(boost::shared_ptr< sts::particle::phylo_node > n);
    void unregister_node(const sts::particle::phylo_node* n);
    bool initialized;

private:
    BeagleInstanceDetails instance_details;
    boost::shared_ptr<bpp::SiteContainer> sites;
    boost::shared_ptr<bpp::SubstitutionModel> model;
    int num_buffers;
    int instance;
    int next_id;
    std::stack<int> free_ids;
    std::unordered_map<const sts::particle::phylo_node*, double> node_ll_map; // caches the root ll at each node
    std::unordered_map<const sts::particle::phylo_node*, int> node_buffer_map; // maps nodes to a beagle buffer ID

    int create_beagle_instance();
    void grow();
    void set_eigen_and_rates_and_weights(int instance);
    void set_eigen_and_rates_and_weights(int, const bpp::SubstitutionModel&);
    int get_buffer( boost::shared_ptr< sts::particle::phylo_node > n );
};

}
}

#endif // STS_LIKELIHOOD_DETAIL_ONLINE_CALCULATOR_FWD_HPP
