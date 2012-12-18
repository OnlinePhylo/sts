#include "log/json_logger.h"
#include "log/sampler.h"

#include "mcmc_event.h"

#include <stdexcept>

using namespace sts::moves;

namespace sts
{
namespace log
{

Json::Value to_value(const Mcmc_event& event)
{
    Json::Value result(Json::objectValue);
    result["name"] = event.name;
    result["time"] = static_cast<double>(event.time);
    result["size"] = event.size;
    result["mh_ratio"] = event.mh_ratio;
    result["accepted"] = event.result == Mcmc_move_result::ACCEPTED;
    return result;
}

/// Construct a Json_logger from an output stream
/// \param out Destination for JSON logging.
Json_logger::Json_logger(std::ostream & out) : out(&out)
{
    root["version"] = "0.1";
    Json::Value generations(Json::arrayValue);
    root["generations"] = Json::Value(Json::arrayValue);
    root["mcmc"] = Json::Value(Json::arrayValue);
}

/// Log the current generation of particles.
///
/// This method <b>does not</b> write to the output stream. Writing is delated until \c Json_logger::write() is called.
/// \param sampler
/// \param node_name map Map from nodes to names, for leaves.
void Json_logger::log(smc::sampler<sts::particle::Particle>& sampler, const std::unordered_map<sts::particle::Node_ptr, std::string>& node_name_map)
{
    Json::Value generation_root;
    sts::log::to_json(sampler, generation_root, node_name_map, particle_id_map, node_id_map);
    root["generations"].append(generation_root);
}

void Json_logger::log_mcmc(const Mcmc_event& event)
{
    root["mcmc"].append(to_value(event));
}

/// Write the root node and all logged generations to \c out.
void Json_logger::write()
{
    Json::StyledStreamWriter writer;
    writer.write(*out, root);
}

}
}
