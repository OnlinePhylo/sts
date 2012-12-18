#ifndef STS_MOVES_MCMC_EVENT_H
#define STS_MOVES_MCMC_EVENT_H

#include <string>

namespace sts
{
namespace moves
{

enum class Mcmc_move_result {
    ACCEPTED,
    REJECTED
};

struct Mcmc_event
{
    std::string name;
    long time;
    double size;
    double mh_ratio;
    Mcmc_move_result result;
};

}
}

#endif // STS_MOVES_MCMC_EVENT_H