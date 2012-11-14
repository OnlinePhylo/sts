#ifndef STS_LOG_SAMPLER_H
#define STS_LOG_SAMPLER_H


#include <unordered_map>
#include <json/json.h>
#include "smctc.hh"
#include "node.h"
#include "state.h"


namespace sts
{
namespace log
{

// Node information indices within serialized JSON
const unsigned int NODE_ID = 0u,
                   NODE_NAME = 1u,
                   NODE_CHILD1 = 2u,
                   NODE_CHILD2 = 3u,
                   NODE_LENGTH1 = 4u,
                   NODE_LENGTH2 = 5u;

// Particle information indices within serialized JSON.
const unsigned int PARTICLE_ID = 0u,
                   PARTICLE_NAME = 1u,
                   PARTICLE_PREDECESSOR = 2u,
                   PARTICLE_PARTIAL_LL = 3u,
                   PARTICLE_FORWARD_LOG_DENSITY = 4u,
                   PARTICLE_BACKWARD_LOG_DENSITY = 5u;

void to_json(smc::sampler<sts::particle::Particle>& sampler,
             Json::Value& root,
             const std::unordered_map < sts::particle::Node_ptr, std::string > &,
             std::unordered_map< sts::particle::Particle, int >&,
             std::unordered_map< sts::particle::Node_ptr, int >&);
}
}

#endif // STS_LOG_SAMPLER_H
