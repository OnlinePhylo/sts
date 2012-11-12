#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Seq/DistanceMatrix.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <algorithm>
#include <vector>
#include "guided_pair_proposer.h"

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
max_tries(max_tries), log_likelihood(ll), upp(ll)
{}

void Guided_pair_proposer::initialize(const std::string& tree_file)
{
  bpp::Newick nwk;
  std::ifstream infile(tree_file.c_str());
  bpp::Tree* tree = nwk.read(infile);
  DistanceMatrix* dm =  TreeTemplateTools::getDistanceMatrix(*tree);
  int leaves = tree->getNumberOfLeaves();
  int pairs = (leaves*(leaves-1))/2;
  sampling_dist.resize(pairs);
  int k=0;
  cumulative = 0;
  for(int i=0; i<leaves; i++){
    for(int j=i+1; j<leaves; j++){
      cumulative += 1.0 / (*dm)(i,j);
      sampling_dist[k].first = cumulative;
      sampling_dist[k].second = make_pair(i,j);
      k++;
    }
  }
  std::sort(sampling_dist.begin(), sampling_dist.end());
  r = gsl_rng_alloc(gsl_rng_mt19937);
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

void Guided_pair_proposer::operator()(particle::Particle pp, smc::rng* rng, particle::Node_ptr& a, particle::Node_ptr& b, double& fwd_density, double& back_density)
{
  int a1, a2;
  double d1;
  int i=0;
  for(; i<max_tries; i++){
    propose(rng, a1, a2, d1);
    // check whether a1 and a2 are unmerged

    break;
    // success!
    return;
  }
  // propose a join uniformly at random
  upp(pp, rng, a, b, fwd_density, back_density);
}

Guided_pair_proposer::~Guided_pair_proposer()
{
  gsl_rng_free(r);
}

} // namespace moves
} // namespace sts
