#ifndef STS_LOG_JSON_LOGGER_HPP
#define STS_LOG_JSON_LOGGER_HPP

#include <iostream>
#include <string>
#include <unordered_map>

#include <json/json.h>
#include "smctc.hh"

namespace sts
{

// Forwards
namespace particle
{
class Particle;
class Node_ptr;
}

namespace log
{

class Json_logger
{
public:
    explicit Json_logger(std::ostream & out);
    void log(smc::sampler<sts::particle::Particle>& sampler, const std::unordered_map<sts::particle::Node_ptr, std::string>& node_name_map);
    void write();
private:
    std::ostream* out;
    std::unordered_map< sts::particle::Particle, int > particle_id_map; // This is a map that gets filled with particle id's when we log?
    std::unordered_map< sts::particle::Node_ptr, int > node_id_map; // This is a map that gets filled with node id's when we log?
    Json::Value root;
};

}
}
#endif // STS_LOG_JSON_LOGGER_HPP
