#ifndef STS_ONLINE_ATTACHMENT_LIKELIHOOD_H
#define STS_ONLINE_ATTACHMENT_LIKELIHOOD_H

#include <memory>
#include <string>

#include "composite_tree_likelihood.h"
#include "beagle_tree_likelihood.h"

namespace sts { namespace online {

class AttachmentLikelihood {
 private:
  CompositeTreeLikelihood ctl_;
  bpp::Node* edge_;
  double distal_length_;
  bool dirty_;

  std::shared_ptr<BeagleTreeLikelihood> btl_;
  BeagleBuffer b1_;
  BeagleBuffer b2_;
  int beagle_instance_;
  int distal_buffer_;
  int proximal_buffer_;
  int leaf_buffer_;
  int scratch1_;
  int scratch2_;

 public:
  AttachmentLikelihood(const CompositeTreeLikelihood& ctl,
                       bpp::Node* edge,
                       const std::string& new_leaf_name,
                       double distal_length);

  void setDistalLength(double distal_length);
  double operator()(double pendant_length);
};

}} // namespace sts::online

#endif // STS_ONLINE_ATTACHMENT_LIKELIHOOD_H
