#ifndef STS_ONLINE_ATTACHMENT_LOCATION_H
#define STS_ONLINE_ATTACHMENT_LOCATION_H

#include <Bpp/Phyl/Node.h>

namespace sts
{
namespace online
{
// An attachment location on a tree
struct AttachmentLocation
{
    AttachmentLocation(bpp::Node* node, double distal) :
        node(node), distal(distal) {};

    bpp::Node* node;
    double distal;
};
}
}


#endif
