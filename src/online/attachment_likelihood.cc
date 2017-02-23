#include "attachment_likelihood.h"

#include <vector>

#include "gsl.h"
#include "util.h"


using sts::util::beagle_check;

namespace sts { namespace online {

AttachmentLikelihood::AttachmentLikelihood(const CompositeTreeLikelihood& ctl) : ctl_(ctl), btl_(ctl_.calculator_)
{
  beagle_instance_ = btl_->beagleInstance();
}
    
void AttachmentLikelihood::initialize(const bpp::Node* edge, const std::string& new_leaf_name, double distal_length)
{
  edge_ = edge;
  distal_length_ = distal_length;
  dirty_ = true;
    
  b1_ = btl_->borrowBuffer();
  b2_ = btl_->borrowBuffer();
  b3_ = btl_->borrowBuffer();
  b4_ = btl_->borrowBuffer();
    
  scratch1_ = b1_->value();
  scratch2_ = b2_->value();
  scratch3_ = b3_->value();
  scratch4_ = b4_->value();
    
  distal_buffer_ = btl_->getDistalBuffer(edge_);
  proximal_buffer_ = btl_->getProximalBuffer(edge_);
  leaf_buffer_ = btl_->getLeafBuffer(new_leaf_name);
}

void AttachmentLikelihood::finalize()
{
  b4_.reset();
  b3_.reset();
  b2_.reset();
  b1_.reset();
}

void AttachmentLikelihood::setDistalLength(double distal_length)
{
  distal_length_ = distal_length;
  dirty_ = true;
}

double AttachmentLikelihood::operator()(double pendant_length)
{
  std::vector<BeagleOperation> operations;
  std::vector<double> branch_lengths;
  std::vector<int> node_indices;

  // If distal changed, update partial
  if (dirty_) {
    operations.push_back(BeagleOperation({scratch1_,
              BEAGLE_OP_NONE,
              BEAGLE_OP_NONE,
              distal_buffer_,
              distal_buffer_,
              proximal_buffer_,
              proximal_buffer_}));
    branch_lengths.push_back(distal_length_);
    node_indices.push_back(distal_buffer_);
    branch_lengths.push_back(edge_->getDistanceToFather() - distal_length_);
    node_indices.push_back(proximal_buffer_);
    dirty_ = false;
  }

  // Always update root partials
  branch_lengths.push_back(pendant_length);
  node_indices.push_back(leaf_buffer_);


  beagle_check(beagleUpdateTransitionMatrices(beagle_instance_,
                                              0,
                                              node_indices.data(),
                                              NULL,
                                              NULL,
                                              branch_lengths.data(),
                                              node_indices.size()));
    
    if(operations.size() > 0) beagle_check(beagleUpdatePartials(beagle_instance_, operations.data(), operations.size(), BEAGLE_OP_NONE));
    
  const int categoryWeightIdx = 0;
  const int stateFreqIdx = 0;
  double logLike;
    
    beagle_check(beagleCalculateEdgeLogLikelihoods(beagle_instance_,
                                                   &scratch1_,
                                                   &leaf_buffer_,
                                                   &leaf_buffer_,
                                                   NULL,
                                                   NULL,
                                                   &categoryWeightIdx,
                                                   &stateFreqIdx,
                                                   NULL,
                                                   1,
                                                   &logLike,
                                                   NULL,
                                                   NULL));

  return logLike;
}
    
double AttachmentLikelihood::derivatives_pendant(double pendant_length, double &d1, double &d2)
{
    beagle_check(beagleUpdateTransitionMatrices(beagle_instance_,
                                                0,
                                                &leaf_buffer_,
                                                &scratch3_,
                                                &scratch4_,
                                                &pendant_length,
                                                1));
    
    const int categoryWeightIdx = 0;
    const int stateFreqIdx = 0;
    double logLike;
    
    beagle_check(beagleCalculateEdgeLogLikelihoods(beagle_instance_,
                                                   &scratch1_,
                                                   &leaf_buffer_,
                                                   &leaf_buffer_,
                                                   &scratch3_,
                                                   &scratch4_,
                                                   &categoryWeightIdx,
                                                   &stateFreqIdx,
                                                   NULL,
                                                   1,
                                                   &logLike,
                                                   &d1,
                                                   &d2));

    return logLike;
}

double AttachmentLikelihood::derivatives_distal(double pendant_length, double &d1, double &d2)
{
    std::vector<BeagleOperation> operations;
    std::vector<double> branch_lengths;
    std::vector<int> node_indices;
    
    // calculate partial for distal and pendant in scratch2
    operations.push_back(BeagleOperation({scratch2_,
        BEAGLE_OP_NONE,
        BEAGLE_OP_NONE,
        leaf_buffer_,
        leaf_buffer_,
        proximal_buffer_,
        proximal_buffer_}));
    branch_lengths.push_back(pendant_length);
    node_indices.push_back(leaf_buffer_);
    branch_lengths.push_back(edge_->getDistanceToFather() - distal_length_);
    node_indices.push_back(proximal_buffer_);
    
    // update matrices of distal and pendant
    beagle_check(beagleUpdateTransitionMatrices(beagle_instance_,
                                                0,
                                                node_indices.data(),
                                                NULL,
                                                NULL,
                                                branch_lengths.data(),
                                                node_indices.size()));
    
    // update matrix and derivatives of distal
    beagle_check(beagleUpdateTransitionMatrices(beagle_instance_,
                                                0,
                                                &distal_buffer_,
                                                &scratch3_,
                                                &scratch4_,
                                                &distal_length_,
                                                1));
    
    beagle_check(beagleUpdatePartials(beagle_instance_, operations.data(), operations.size(), NULL));
    
    const int categoryWeightIdx = 0;
    const int stateFreqIdx = 0;
    double logLike;
    
    beagle_check(beagleCalculateEdgeLogLikelihoods(beagle_instance_,
                                                   &scratch2_,
                                                   &distal_buffer_,
                                                   &distal_buffer_,
                                                   &scratch3_,
                                                   &scratch4_,
                                                   &categoryWeightIdx,
                                                   &stateFreqIdx,
                                                   NULL,
                                                   1,
                                                   &logLike,
                                                   &d1,
                                                   &d2));

    return logLike;
}

}} // namespace sts::online
