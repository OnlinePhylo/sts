#ifndef STS_LOG_JSON_LOGGER_HPP
#define STS_LOG_JSON_LOGGER_HPP

#include <iostream>
#include <string>
#include <json/json.h>
#include <sts/log/sampler.hpp>

class json_logger
{
public:
    void initialize(std::ostream& out);
    void log(smc::sampler<sts::particle::particle>& sampler, const std::unordered_map<sts::particle::node, std::string>& node_name_map);
private:
    std::ostream* out;
    std::unordered_map< sts::particle::particle, int > particle_id_map; // This is a map that gets filled with particle id's when we log?
    std::unordered_map< sts::particle::node, int > node_id_map; // This is a map that gets filled with node id's when we log?
    void write_object(Json::Value &object);
};

void json_logger::write_object(Json::Value &object)
{
    Json::StyledStreamWriter writer;
    writer.write(*out, object);
}

void json_logger::initialize(std::ostream& out_ref)
{
    out = &out_ref;
    Json::Value root;
    root["version"] = "0.1";
    write_object(root);
}

void json_logger::log(smc::sampler<sts::particle::particle>& sampler, const std::unordered_map<sts::particle::node, std::string>& node_name_map)
{
    Json::Value root;
    sts::log::to_json(sampler, root, node_name_map, particle_id_map, node_id_map);
    write_object(root);
}

#endif // STS_LOG_JSON_LOGGER_HPP
