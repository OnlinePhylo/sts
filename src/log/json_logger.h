#ifndef STS_LOG_JSON_LOGGER_H
#define STS_LOG_JSON_LOGGER_H
#include "metropolis_hastings_move.h"
#include "particle.h"
#include "node_ptr.h"
#include <iostream>
#include <string>
#include <unordered_map>

#include <json/json.h>
#include "smctc.hh"

namespace sts
{

// Forwards
namespace moves { struct Mcmc_event; }

namespace log
{

class Json_logger
{
public:
    explicit Json_logger(std::ostream & out);
    void log(smc::sampler<sts::particle::Particle>& sampler, const std::unordered_map<sts::particle::Node_ptr, std::string>& node_name_map);

    void log_mcmc(const sts::moves::Mcmc_event& event);
    void write();
private:
    std::ostream* out;
    std::unordered_map< sts::particle::Particle, int > particle_id_map; // This is a map that gets filled with particle id's when we log?
    std::unordered_map< sts::particle::Node_ptr, int > node_id_map; // This is a map that gets filled with node id's when we log?
    Json::Value root;
};

}
}
#endif // STS_LOG_JSON_LOGGER_H
