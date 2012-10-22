#ifndef STS_LIKELIHOOD_DETAIL_ONLINE_CALCULATOR_FWD_HPP
#define STS_LIKELIHOOD_DETAIL_ONLINE_CALCULATOR_FWD_HPP

#include <memory>
#include <stack>
#include <unordered_map>

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
    void initialize(std::shared_ptr<bpp::SiteContainer>, std::shared_ptr<bpp::SubstitutionModel>);
    int get_id();
    void free_id(int id);
    double calculate_ll(std::shared_ptr<sts::particle::phylo_node> node, std::vector<bool>& visited);
    void invalidate(int);
    bool initialized;

private:
    BeagleInstanceDetails instance_details;
    std::shared_ptr<bpp::SiteContainer> sites;
    std::shared_ptr<bpp::SubstitutionModel> model;
    int num_buffers;
    int instance;
    int next_id;
    std::stack<int> free_ids;
    std::unordered_map<int, double> map_id_ll; // caches the root ll at each ID

    int create_beagle_instance();
    void grow();
    void set_eigen_and_rates_and_weights(int instance);
    void set_eigen_and_rates_and_weights(int, const bpp::SubstitutionModel&);
};

}
}

#endif // STS_LIKELIHOOD_DETAIL_ONLINE_CALCULATOR_FWD_HPP
