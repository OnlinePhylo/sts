#ifndef STS_ONLINE_ATTACHMENT_LIKELIHOOD_H
#define STS_ONLINE_ATTACHMENT_LIKELIHOOD_H

#include <memory>
#include <string>

#include "composite_tree_likelihood.h"
#include "beagle_tree_likelihood.h"

namespace sts { namespace online {

class AttachmentLikelihood {
 private:
  const CompositeTreeLikelihood& ctl_;
  const bpp::Node* edge_;
  double distal_length_;
  bool dirty_;

  std::shared_ptr<BeagleTreeLikelihood> btl_;
  std::unique_ptr<BeagleBuffer> b1_;
  std::unique_ptr<BeagleBuffer> b2_;
    
  std::unique_ptr<BeagleBuffer> b3_;
  std::unique_ptr<BeagleBuffer> b4_;
  int beagle_instance_;
  int distal_buffer_;
  int proximal_buffer_;
  int leaf_buffer_;
  int scratch1_;
  int scratch2_;
  int scratch3_;
  int scratch4_;

 public:
  AttachmentLikelihood(const CompositeTreeLikelihood& ctl);
    
    void initialize(const bpp::Node* edge, const std::string& new_leaf_name, double distal_length);
    
    void finalize();

  void setDistalLength(double distal_length);
  double operator()(double pendant_length);
    
    double derivatives_pendant(double pendant_length, double &d1, double &d2);
    
    double derivatives_distal(double pendant_length, double &d1, double &d2);
};

}} // namespace sts::online

#endif // STS_ONLINE_ATTACHMENT_LIKELIHOOD_H
