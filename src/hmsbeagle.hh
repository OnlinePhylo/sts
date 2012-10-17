/// \file hmsbeagle.hh
/// \author metatangle, inc.
/// \brief Interface between STS and beagle-lib

#ifndef __hmsbeagle__
#define __hmsbeagle__

#include <string>
#include <vector>
#include <stack>
#include <unordered_map>

#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include "libhmsbeagle/beagle.h"
#include "phylofunc.hh"

// XXX should we come up with a convention for method order? I'd be happy with either approximate dependency, split into
// public and private or not.
class online_calculator
{
public:
    online_calculator() : instance(-1), next_id(0), initialized(false) {}; // XXX Fix for transition to AA
    ~online_calculator() {
        if(instance >= 0)
            beagleFinalizeInstance(instance);
    };

    void initialize(std::shared_ptr<bpp::SiteContainer>, std::shared_ptr<bpp::SubstitutionModel>);
    int get_id();
    void free_id(int id);
    double calculate_ll(std::shared_ptr< phylo_node > node, std::vector<bool>& visited);
    bool initialized;

private:
    BeagleInstanceDetails instance_details;
    std::shared_ptr<bpp::SiteContainer> sites;
    std::shared_ptr<bpp::SubstitutionModel> model;
    int num_buffers;
    int instance;
    int next_id;
    std::stack<int> free_ids;
    std::unordered_map< int, double > map_id_ll; // caches the root ll at each ID

    int create_beagle_instance();
    void grow();
    void set_eigen_and_rates_and_weights(int instance);
    void set_eigen_and_rates_and_weights(int, const bpp::SubstitutionModel&);
};

#endif //  __hmsbeagle__
