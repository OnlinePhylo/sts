#include "log/sampler.h"

#include "edge.h"

#include <memory>
#include <stack>
#include <unordered_set>

namespace sts
{
namespace log
{

/// Serialize a particle.
/// \param sampler The sampler.
/// \param root The root of the JSON object.
/// \param node_name_map A map from the nodes to strings;
/// \param particle_id_map This is a map that gets filled with particle id's when we log?
/// \param node_id_map This is a map that gets filled with particle id's when we log?
void to_json(smc::sampler<sts::particle::Particle>& sampler,
             Json::Value& root,
             const std::unordered_map < sts::particle::Node_ptr, std::string > & node_name_map,
             std::unordered_map< sts::particle::Particle, int >& particle_id_map,
             std::unordered_map< sts::particle::Node_ptr, int >& node_id_map)
{

    root["generation"] = (int)sampler.GetTime();
    Json::Value& particles = root["particles"];
    Json::Value& states = root["states"];
    states["fields"][PARTICLE_ID                  ] = "id";
    states["fields"][PARTICLE_NAME                ] = "node";
    states["fields"][PARTICLE_PREDECESSOR         ] = "predecessor";
    states["fields"][PARTICLE_LL                  ] = "log_like";
    states["fields"][PARTICLE_FORWARD_LOG_DENSITY ] = "forward_log_density";
    states["fields"][PARTICLE_BACKWARD_LOG_DENSITY] = "back_log_density";
    Json::Value& nodes = root["nodes"];
    nodes["fields"][NODE_ID     ] = "id";
    nodes["fields"][NODE_NAME   ] = "name";
    nodes["fields"][NODE_CHILD1 ] = "child1";
    nodes["fields"][NODE_CHILD2 ] = "child2";
    nodes["fields"][NODE_LENGTH1] = "length1";
    nodes["fields"][NODE_LENGTH2] = "length2";
    std::unordered_set< sts::particle::Particle > particles_visited;
    std::unordered_set< sts::particle::Node_ptr > nodes_visited;

    int nindex = 0; // node index
    int pindex = 0; // particle index

    // Add all the leaf nodes if they haven't already been added.
    // XXX Are the only named nodes leaves? If so I think we should call this leaf_name_map everywhere.
    for(auto node_name_pair : node_name_map) { // Iterate over (key,value) pairs in node_name_map.
        if(nodes_visited.count(node_name_pair.first)) continue;
        nodes_visited.insert(node_name_pair.first);
        if(node_id_map.count(node_name_pair.first) == 0) node_id_map[node_name_pair.first] = node_id_map.size();
        int nid = node_id_map[node_name_pair.first];
        Json::Value& jnode = nodes["data"][nindex++];
        jnode[NODE_ID] = nid;
        jnode[NODE_NAME] = node_name_pair.second;
    }

    // Traverse the particle system and add particles and any internal tree nodes.
    for(int i = 0; i < sampler.GetNumber(); i++) {
        // determine whether we've seen this particle previously and if not add it
        sts::particle::Particle X = sampler.GetParticleValue(i);
        std::stack<sts::particle::Particle> s;
        s.push(X);

        while(s.size() > 0) {
            sts::particle::Particle p = s.top();
            s.pop();
            if(p == nullptr) continue;
            // Skip if we've already seen this particle.
            if(particles_visited.count(p) != 0) continue;
            particles_visited.insert(p);
            sts::particle::Particle pred = p->predecessor;
            // Insert particle ID's if needed.
            if(particle_id_map.count(p) == 0) particle_id_map[p] = particle_id_map.size();
            if(particle_id_map.count(pred) == 0) particle_id_map[pred] = particle_id_map.size();
            int pid = particle_id_map[p];
            Json::Value& jpart = states["data"][pindex++];
            jpart[PARTICLE_ID] = pid;
            jpart[PARTICLE_LL] = p->log_likelihood;
            jpart[PARTICLE_FORWARD_LOG_DENSITY] = p->forward_log_density;
            jpart[PARTICLE_BACKWARD_LOG_DENSITY] = p->forward_log_density;
            if(pred != nullptr) jpart[PARTICLE_PREDECESSOR] = particle_id_map[pred];
            if(pred != nullptr) s.push(pred);

            // Traverse the nodes below this particle.
            std::stack<sts::particle::Node_ptr> ns;
            if(p->node != nullptr) ns.push(p->node);
            while(ns.size() > 0) {
                sts::particle::Node_ptr n = ns.top();
                ns.pop();
                if(n == nullptr) continue;
                if(nodes_visited.count(n) != 0) continue;
                nodes_visited.insert(n);
                if(node_id_map.count(n) == 0) node_id_map[n] = node_id_map.size();
                // Determine whether we've seen this node previously and if not add it.
                int nid = node_id_map[n];
                Json::Value& jnode = nodes["data"][nindex++];

                sts::particle::Node_ptr c1 = n->child1->node;
                sts::particle::Node_ptr c2 = n->child2->node;
                ns.push(c1);
                ns.push(c2);
                // Insert node ID's if needed.
                if(node_id_map.count(c1) == 0) node_id_map[c1] = node_id_map.size();
                if(node_id_map.count(c2) == 0) node_id_map[c2] = node_id_map.size();

                jnode[NODE_ID     ] = nid;
                jnode[NODE_CHILD1 ] = node_id_map[c1];
                jnode[NODE_CHILD2 ] = node_id_map[c2];
                jnode[NODE_LENGTH1] = n->child1->length;
                jnode[NODE_LENGTH2] = n->child2->length;
            }

            // Add a node reference to this particle.
            // XXX ??? Was 1u. Is PARTICLE_NAME correct?
            jpart[PARTICLE_NAME] = node_id_map[p->node];
        }
        particles[i] = particle_id_map[X];
    }
}

}
}
