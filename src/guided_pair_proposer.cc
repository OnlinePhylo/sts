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

Guided_pair_proposer::Guided_pair_proposer( sts::likelihood::Forest_likelihood* ll) :
 strength(2.0), r(NULL), log_likelihood(ll)
{}

void Guided_pair_proposer::initialize(const std::string& tree_file, const std::unordered_map<particle::Node_ptr, std::string>& node_name_map)
{
  bpp::Newick nwk;
  std::ifstream infile(tree_file.c_str());
  bpp::Tree* tree = nwk.read(infile);
  dm =  TreeTemplateTools::getDistanceMatrix(*tree);
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
}

void Guided_pair_proposer::propose(smc::rng *rng, int& leaf1, int& leaf2, double& density, std::vector< std::pair< double, std::pair< int, int > > >& distribution)
{
    double c;
    if(rng != NULL) c = rng->Uniform(0, distribution.back().first);
    if(rng == NULL) c = gsl_ran_flat(r, 0, distribution.back().first);
    pair< double, pair< int, int > > cpair;
    cpair.first = c;
    auto iter = std::lower_bound<>(distribution.begin(), distribution.end(), cpair);
    leaf1 = iter->second.first;
    leaf2 = iter->second.second;
    double prev = 0;
    if(iter != distribution.begin()){
      auto p = iter-1;
      prev = p->first;
    }
    density = (iter->first - prev)/distribution.back().first;
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
  vector<particle::Node_ptr> uncoalesced;
  uncoalesced = sts::util::uncoalesced_nodes(pp, log_likelihood->get_leaves());
  // construct a sampling distribution
  std::vector< std::pair< double, std::pair< int, int > > > distribution;
  double mass = 0;
  for(unsigned i=0; i<uncoalesced.size(); i++){
    for(unsigned j=i+1; j<uncoalesced.size(); j++){
      unordered_set<particle::Node_ptr> n1_leaves, n2_leaves;
      leaf_nodes(uncoalesced[i], n1_leaves);
      leaf_nodes(uncoalesced[j], n2_leaves);
      // use the shortest distance between leaves as the sampling weight for this node pair
      // this approach will have some major problems deeper in the tree
      // a smarter refinement might subtract the distance from root to leaf in each subtree from
      // the pairwise distance and take the average of these values among all leaf pairs.
      double max_mass = 0;
      for(auto l1 : n1_leaves){
	for(auto l2 : n2_leaves){
	  double cur_weight = get_weight(node_dm_id_map[l1], node_dm_id_map[l2]);
	  max_mass = max(max_mass, cur_weight);
	}
      }
      mass += max_mass;
      std::pair< int, int > node_pair(i, j);
      distribution.push_back(make_pair(mass, node_pair));
    }
  }
  
  int a1, a2;
  double d1;
  propose(rng, a1, a2, d1, distribution);
  a = uncoalesced[a1];
  b = uncoalesced[a2];
  fwd_density = d1;
    
  // Calculate the back-proposal density
  // This is the density of all other ways to propose the same configuration
  // from all possible previous configurations  
  // This is simply the number of previous states that lead to the current forest
  int tc = util::count_uncoalesced_trees(uncoalesced);
  back_density = tc > 1.0 ? tc : 1.0;
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
