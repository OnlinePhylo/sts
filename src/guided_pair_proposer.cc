#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <algorithm>
#include <vector>
#include "guided_pair_proposer.h"
#include "util.h"

/**
Guided merges, at initialization:
1. calc the inverse of all pairwise dists z_{i,j} = \frac{1}{d_{i,j}}
2. construct a sampling distribution from them as p(i,j) = \frac{z_{i,j}}{\sum_{x \neq y}{z_{x,y}}}
3. sort the p(i,j) so that sampling w/replacement can be done in log(n) time with binary search

Proposing a join:
1. sample a pair i,j from the sampling distribution
2. if i,j are part of the same tree already, pick a join among forest roots uniformly at random (or some other move)
3. else if i,j are not part of the same tree, merge the roots of the trees 
4. the forward density is the sum of p(i,j) from all leaf cross-pairs from the merged trees + (prob of picking an invalid merge * fwd prob in other move)
5. the backward density is the sum of p(i,j) from all leaf cross pairs below all forest roots + (prob of picking an invalid merge * backward prob in other move). 

Implementation details:
For each particle, store a union-find structure mapping nodes to their current root. 
This will facilitate efficient lookup for whether leaf pairs are in the same tree.
keep a running total of backward prob at each successive particle

*/

using namespace bpp;
using namespace std;

namespace sts
{
namespace moves
{

Guided_pair_proposer::Guided_pair_proposer( int max_tries, sts::likelihood::Forest_likelihood& ll) :
max_tries(max_tries), strength(2.0), r(NULL), log_likelihood(ll), upp(ll)
{}

void Guided_pair_proposer::initialize(const std::string& tree_file, const std::unordered_map<particle::Node_ptr, std::string>& node_name_map)
{
  assert(max_tries == 1);
  bpp::Newick nwk;
  std::ifstream infile(tree_file.c_str());
  bpp::Tree* tree = nwk.read(infile);
  dm =  TreeTemplateTools::getDistanceMatrix(*tree);
  int leaves = tree->getNumberOfLeaves();
  int pairs = (leaves*(leaves-1))/2;
  sampling_dist.resize(pairs);
  int k=0;
  cumulative = 0;
  for(int i=0; i<leaves; i++){
    for(int j=i+1; j<leaves; j++){
      cumulative += get_weight(i,j);
      sampling_dist[k].first = cumulative;
      sampling_dist[k].second = make_pair(i,j);
      k++;
    }
  }
  // Construct mappings from node pointer to distance matrix entry and vice-versa.
  unordered_map<string, int> name_dm_id_map;
  for(unsigned int i=0; i<dm->size(); i++){
    name_dm_id_map[dm->getName(i)]=i;
  }
  for(auto& n : node_name_map){
    node_dm_id_map[n.first] = name_dm_id_map[n.second];
  }
  for(auto& n : node_dm_id_map){
    dm_id_node_map[n.second] = n.first;
  }

//  r = gsl_rng_alloc(gsl_rng_mt19937);
}

void Guided_pair_proposer::propose(smc::rng *rng, int& leaf1, int& leaf2, double& density)
{
    int c;
    if(rng != NULL) c = rng->Uniform(0, cumulative);
    if(rng == NULL) c = gsl_ran_flat(r, 0, cumulative);
    pair< double, pair< int, int > > cpair;
    cpair.first = c;
    auto iter = std::lower_bound<>(sampling_dist.begin(), sampling_dist.end(), cpair);
    leaf1 = iter->second.first;
    leaf2 = iter->second.second;
    double prev = 0;
    if(iter != sampling_dist.begin()){
      auto p = iter-1;
      prev = p->first;
    }
    density = (iter->first - prev)/cumulative;
}

void leaf_nodes(particle::Node_ptr& node, unordered_set< particle::Node_ptr >& leaves){
  std::stack< particle::Node_ptr > s;
  s.push(node);
  while(s.size()>0){
    particle::Node_ptr cur = s.top();
    s.pop();
    if(cur->child1 != NULL){
      s.push(cur->child1->node);
      s.push(cur->child2->node);
    }else{
      leaves.insert(cur);
    }
  }
}

double Guided_pair_proposer::get_weight(int i, int j) const
{
  return 1.0 / pow((*dm)(i,j), strength);
}

void Guided_pair_proposer::operator()(particle::Particle pp, smc::rng* rng, particle::Node_ptr& a, particle::Node_ptr& b, double& fwd_density, double& back_density)
{
  int a1, a2;
  double d1;
  int i=0;
  double invalid_probability;
  vector<particle::Node_ptr> uncoalesced;
  uncoalesced = sts::util::uncoalesced_nodes(pp, uncoalesced);
  for(; i<max_tries; i++){
    propose(rng, a1, a2, d1);
    particle::Node_ptr n1 = dm_id_node_map[a1];
    particle::Node_ptr n2 = dm_id_node_map[a2];
    // Check whether n1 and n2 are in different trees in the forest.
    // Do this by first accumulating all forest roots that are not leaf nodes
    // then asking whether n1 or n2 are descendants.
    particle::Node_ptr n1_root = n1, n2_root = n2;
    unordered_set<particle::Node_ptr> all_leaves;
    for(auto u : uncoalesced){
      unordered_set<particle::Node_ptr> leaves;
      leaf_nodes(u, leaves);
      all_leaves.insert(leaves.begin(), leaves.end());
      if(leaves.count(n1)) n1_root = u;
      if(leaves.count(n2)) n2_root = u;
    }
    // Calculate the probability of proposing an invalid merge
    invalid_probability = 0;
    for(auto i = all_leaves.begin(); i != all_leaves.end(); i++){
      auto j = i;
      for(j++; j != all_leaves.end(); j++){	
	invalid_probability += get_weight(node_dm_id_map[*i], node_dm_id_map[*j]);
      }
    }
    invalid_probability /= cumulative;

    
    if(n1_root==n2_root)
      continue;
    // for now only permit merges of leaves
    if(n1!=n1_root)
      continue;
    if(n2!=n2_root)
      continue;
    a = n1;
    b = n2;

    // success! we found a pair to merge.
    // Calculate the probability of proposing this merge
    fwd_density = d1; // valid only when both nodes are leaves
    back_density = d1;
    break;
  }
  if(i==max_tries){
    // Failed to propose a join among leaves. 
    // Instead propose a join uniformly at random.
    upp(pp, rng, a, b, fwd_density, back_density);
    fwd_density *= invalid_probability;

    // if this is a leaf pair, add in the probability of proposing this merge via a weighted proposal
    if(node_dm_id_map.count(a) && node_dm_id_map.count(b)){
      fwd_density += get_weight(node_dm_id_map[a], node_dm_id_map[b]) / cumulative;
    }
    // TODO: implement max_tries > 1
    //  fwd_density *= pow(invalid_probability, max_tries);
  }
    
  // Calculate the back-proposal density
  // This is the density of all other ways to propose the same configuration
  // from all possible previous configurations
  
  // This is the weight of all cherries plus the joint probability of invalid merge and the back proposal density from the uniform merge
  for(auto u : uncoalesced){
    unordered_set<particle::Node_ptr> leaves;
    leaf_nodes(u, leaves);
    // TODO: do we need to factor in 1.0 over the number of active roots?
    if(leaves.size() > 2){
      // when the tree root is not a cherry, the only way to recover it is by
      // proposing it after proposing an invalid merge
      back_density += (invalid_probability + d1);
    }else{
      // for cherries, exclude the possibility of the weighted merge proposal
      particle::Node_ptr c1 = u->child1->node;
      particle::Node_ptr c2 = u->child2->node;
      back_density += invalid_probability + d1 - get_weight(node_dm_id_map[c1], node_dm_id_map[c2])/cumulative;
    }
  }
  if(a == nullptr || b == nullptr){
    cerr << "ruh roh\n";
  }

}

Guided_pair_proposer::~Guided_pair_proposer()
{
  if(r!=NULL) gsl_rng_free(r);
}

} // namespace moves
} // namespace sts
