#include "attachment_likelihood.h"

#include <vector>

#include "gsl.h"
#include "util.h"

namespace sts { namespace online {

AttachmentLikelihood::AttachmentLikelihood(const CompositeTreeLikelihood& ctl,
                                           bpp::Node* edge,
                                           const std::string& new_leaf_name,
                                           double distal_length)
    : ctl_(ctl),
      edge_(edge),
      distal_length_(distal_length),
      dirty_(true),
      btl_(ctl_.calculator_),
      b1_(btl_->borrowBuffer()),
      b2_(btl_->borrowBuffer())
{
  beagle_instance_ = btl_->beagleInstance();
  scratch1_ = b1_.value();
  scratch2_ = b2_.value();
  distal_buffer_ = btl_->getDistalBuffer(edge_);
  proximal_buffer_ = btl_->getProximalBuffer(edge_);
  leaf_buffer_ = btl_->getLeafBuffer(new_leaf_name);
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
  }

  // Always update root partials
  operations.push_back(BeagleOperation({scratch2_,
            BEAGLE_OP_NONE,
            BEAGLE_OP_NONE,
            scratch1_,
            scratch1_,
            leaf_buffer_,
            leaf_buffer_}));
  branch_lengths.push_back(0.0);
  node_indices.push_back(scratch1_);
  branch_lengths.push_back(pendant_length);
  node_indices.push_back(leaf_buffer_);

  // Usual thing

  using sts::util::beagle_check;

  beagle_check(beagleUpdateTransitionMatrices(beagle_instance_,
                                              0,
                                              node_indices.data(),
                                              NULL,
                                              NULL,
                                              branch_lengths.data(),
                                              node_indices.size()));
  beagle_check(beagleUpdatePartials(beagle_instance_, operations.data(), operations.size(), scratch2_));

  std::vector<int> scale_indices(operations.size());
  for (size_t i = 0; i < operations.size(); ++i)
    scale_indices[i] = operations[i].destinationPartials;

  beagle_check(beagleAccumulateScaleFactors(beagle_instance_, scale_indices.data(), scale_indices.size(),
                                            scratch2_));
  const int categoryWeightIdx = 0;
  const int stateFreqIdx = 0;
  double logLike;
  beagle_check(beagleCalculateRootLogLikelihoods(beagle_instance_,
                                                 &scratch2_,
                                                 &categoryWeightIdx,
                                                 &stateFreqIdx,
                                                 &scratch2_,
                                                 1,
                                                 &logLike));

  logLike += ctl_.sumAdditionalLogLikes();

  return logLike;
}

}} // namespace sts::online
