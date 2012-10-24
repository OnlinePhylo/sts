#ifndef STS_PARTICLE_DETAIL_EDGE_FWD_HPP
#define STS_PARTICLE_DETAIL_EDGE_FWD_HPP

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <memory>
#include <string>
#include <unordered_map>
#include <Bpp/Phyl/Tree.h>


namespace std {

template<typename T>
struct hash< boost::shared_ptr< T > >
{
  inline std::size_t operator()(const boost::shared_ptr<T>& p) const
  {
    return reinterpret_cast<size_t>( p.get() );
  }
};

}

namespace sts
{
namespace particle
{
class phylo_node;

/// \class edge
/// An edge
class edge
{
public:
    /// Initialize with a node and distance.
    edge(boost::shared_ptr<phylo_node>, double);

    double length;
    boost::shared_ptr<phylo_node> node;

    /// Make an edge from a bpp Tree and node number
    static boost::shared_ptr<edge> of_tree(boost::shared_ptr<sts::likelihood::online_calculator>, bpp::TreeTemplate<bpp::Node> &, int, std::unordered_map<boost::shared_ptr<phylo_node>, std::string>&);
};
}
}
#endif
