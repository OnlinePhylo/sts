#ifndef STS_LOG_SAMPLER_HPP
#define STS_LOG_SAMPLER_HPP

#include <stack>
#include <memory>
#include <unordered_map>
#include <json/json.h>
#include "sts/particle.hpp"
#include "smctc.hh"

namespace sts {
namespace log {

void to_json(smc::sampler<sts::particle::particle>& sampler, Json::Value& root, const std::unordered_map<sts::particle::node, std::string>& node_name_map, std::unordered_map< sts::particle::particle, int >& particle_id_map, std::unordered_map< sts::particle::node, int >& node_id_map)
{
    root["generation"] = (int)sampler.GetTime();
    Json::Value& states = root["particles"];
    Json::Value& particles = root["states"];
    particles["fields"][0u]="id";
    particles["fields"][1u]="node";
    particles["fields"][2u]="predecessor";
    Json::Value& nodes = root["nodes"];
    nodes["fields"][0u]="id";
    nodes["fields"][1u]="name";
    nodes["fields"][2u]="child1";
    nodes["fields"][3u]="child2";
    nodes["fields"][4u]="length1";
    nodes["fields"][5u]="length2";
    std::unordered_set< sts::particle::particle > particles_visited;
    std::unordered_set< sts::particle::node > nodes_visited;

    int nindex = 0;
    int pindex = 0;

    // Add all the leaf nodes if they haven't already been added.
    for(auto n : node_name_map){
        if(nodes_visited.count(n.first)) continue;
        nodes_visited.insert(n.first);
        if(node_id_map.count(n.first) == 0) node_id_map[n.first] = node_id_map.size();
        int nid = node_id_map[n.first];
        Json::Value& jnode = nodes["data"][nindex++];
        jnode[0u]=nid;
        jnode[1u]=n.second;
    }

    // Traverse the particle system and add particles and any internal tree nodes.
    for(int i=0;i<sampler.GetNumber(); i++){
        // determine whether we've seen this particle previously and if not add it
        sts::particle::particle X = sampler.GetParticleValue(i);
        std::stack<sts::particle::particle> s;
        s.push(X);

        while(s.size()>0){
            sts::particle::particle x = s.top();
            s.pop();
            if(x==NULL) continue;
            if(particles_visited.count(x) != 0) continue;
            particles_visited.insert(x);
            sts::particle::particle pred = x->predecessor;
            if(particle_id_map.count(x) == 0) particle_id_map[x] = particle_id_map.size();
            if(particle_id_map.count(pred) == 0) particle_id_map[pred] = particle_id_map.size();
            int pid = particle_id_map[x];
            Json::Value& jpart = particles["data"][pindex++];
            jpart[0u] = pid;
            if(pred!=NULL) jpart[2u] = particle_id_map[pred];
            if(pred!=NULL) s.push(pred);

            // traverse the nodes below this particle
            std::stack<sts::particle::node> ns;
            if(x->node != NULL) ns.push(x->node);
            while(ns.size()>0){
                sts::particle::node n = ns.top();
                ns.pop();
                if(n==NULL) continue;
                if(nodes_visited.count(n) != 0) continue;
                nodes_visited.insert(n);
                if(node_id_map.count(n) == 0) node_id_map[n] = node_id_map.size();
                // Determine whether we've seen this node previously and if not add it
                int nid = node_id_map[n];
                Json::Value& jnode = nodes["data"][nindex++];

                sts::particle::node c1 = n->child1->node;
                sts::particle::node c2 = n->child2->node;
                ns.push(c1);
                ns.push(c2);
                if(node_id_map.count(c1) == 0) node_id_map[c1] = node_id_map.size();
                if(node_id_map.count(c2) == 0) node_id_map[c2] = node_id_map.size();

                jnode[0u]=nid;
                jnode[2u]=node_id_map[c1];
                jnode[3u]=node_id_map[c2];
                jnode[4u]=n->child1->length;
                jnode[5u]=n->child2->length;
            }

            // Add a node reference to this particle.
            jpart[1u]=node_id_map[x->node];
        }
        states[i] = particle_id_map[X];
    }
}

}
}

#endif // STS_LOG_SAMPLER_HPP
