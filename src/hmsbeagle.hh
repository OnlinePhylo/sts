#ifndef __hmsbeagle__
#define __hmsbeagle__

#include <string>
#include <vector>
#include <stack>
#include <unordered_map>

#include <Bpp/Phyl/Model/AbstractSubstitutionModel.h>

#include "libhmsbeagle/beagle.h"
#include "phylofunc.hh"

std::vector< double > get_partials(const std::string& sequence);

// XXX should we come up with a convention for method order? I'd be happy with either approximate dependency, split into
// public and private or not.
class OnlineCalculator
{
public:
    OnlineCalculator() : stateCount(4), instance(-1), next_id(0) {}; // XXX Fix for transition to AA
    ~OnlineCalculator() {
        beagleFinalizeInstance(instance);
    };

    void initialize(const std::vector<std::string>& seqs);
    int get_id();
    void free_id(int id);
    double calculate_ll(std::shared_ptr< phylo_node > node, std::vector<bool>& visited);

private:
    BeagleInstanceDetails instDetails;
    int stateCount;
    int nPatterns;
    int nPartBuffs;
    int instance;
    int next_id;
    std::stack<int> free_ids;
    // AD: Is this a naming convention of src_dst for maps? Was surprising for me. Could I change to map_id_ll?
    std::unordered_map< int, double > id_ll; // caches the root ll at each ID

    int create_beagle_instance();
    void grow();
    void set_eigen_and_rates_and_weights(int instance);
    void set_eigen_and_rates_and_weights(int instance, const bpp::AbstractReversibleSubstitutionModel& model);
};

#endif //  __hmsbeagle__
