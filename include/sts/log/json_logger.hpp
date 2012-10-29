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
    std::unordered_map< sts::particle::particle, int > particle_id_map;
    std::unordered_map< sts::particle::node, int > node_id_map;
};

void json_logger::initialize(std::ostream& out) 
{
    this->out = &out;
    Json::Value root;
    root["version"]="0.1";
    Json::StyledStreamWriter writer;
    writer.write(out, root);
    out << std::endl;
}

void json_logger::log(smc::sampler<sts::particle::particle>& sampler, const std::unordered_map<sts::particle::node, std::string>& node_name_map)
{
    if(out==NULL) return;
    Json::Value root;
    Json::StyledWriter writer;
    sts::log::to_json( sampler, root, node_name_map, particle_id_map, node_id_map);
    std::string json_string = writer.write(root);
    (*out) << json_string << std::endl << std::endl;            
}

#endif // STS_LOG_JSON_LOGGER_HPP
