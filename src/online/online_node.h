#ifndef STS_ONLINE_ONLINE_NODE_H
#define STS_ONLINE_ONLINE_NODE_H

#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

namespace sts { namespace online {

class Online_node : public bpp::Node
{
friend class bpp::TreeTemplateTools;

private:
    bool dirty;

public:
    Online_node() : bpp::Node(), dirty(true) {}
    Online_node(int id) : bpp::Node(id), dirty(true) {}
    Online_node(const std::string& name) : bpp::Node(name), dirty(true) {}
    Online_node(int id, const std::string& name) : bpp::Node(id, name), dirty(true) {}

protected:
    Online_node(const bpp::Node& node) : bpp::Node(node), dirty(true) {}
    Online_node(const Online_node& node):
        bpp::Node(node), dirty(node.dirty)
    {}

    Online_node& operator=(const Online_node& node)
    {
        bpp::Node::operator=(node);
        dirty = node.dirty;
        return *this;
    }

    Online_node* clone() const { return new Online_node(*this); }
public:
    virtual ~Online_node() {}

    const Online_node* getFather() const { return dynamic_cast<const Online_node *>(father_); }
    Online_node* getFather() { return dynamic_cast<Online_node *>(father_); }
    Online_node* removeFather() {
        make_dirty();
        Online_node* f = dynamic_cast<Online_node *>(father_);
        father_ = 0;
        return f;
    }
    const Online_node* getSon(size_t i) const throw (bpp::IndexOutOfBoundsException)
    {
        return dynamic_cast<Online_node *>(sons_[i]);
    }
    Online_node* getSon(size_t i) throw (bpp::IndexOutOfBoundsException)
    {
        return dynamic_cast<Online_node *>(sons_[i]);
    }
    std::vector<const Online_node*> getNeighbors() const
    {
        std::vector<const Node*> neighbors = bpp::Node::getNeighbors();
        std::vector<const Online_node*> neighbors2(neighbors.size());
        for (size_t i = 0; i < neighbors.size(); i++)
            neighbors2[i] = dynamic_cast<const Online_node*>(neighbors[i]);
        return neighbors2;
    }
    std::vector<Online_node*> getNeighbors()
    {
        std::vector<Node*> neighbors = bpp::Node::getNeighbors();
        std::vector<Online_node*> neighbors2(neighbors.size());
        for (size_t i = 0; i < neighbors.size(); i++)
            neighbors2[i] = dynamic_cast<Online_node*>(neighbors[i]);
        return neighbors2;
    }

    Online_node* operator[](int i)
    {
        return dynamic_cast<Online_node*>((i < 0) ? father_ : sons_[i]);
    }

    const Online_node* operator[](int i) const
    {
        return dynamic_cast<const Online_node *>((i < 0) ? father_ : sons_[i]);
    }

    inline void make_dirty() { dirty = true; }
    inline void make_clean() { dirty = false; }
    inline bool is_dirty() const { return dirty; }

    void setDistanceToFather(double distance)
    {
        make_dirty();
        bpp::Node::setDistanceToFather(distance);
    }
    void deleteDistanceToFather()
    {
        make_dirty();
        bpp::Node::deleteDistanceToFather();
    }

};

}}

#endif  //STS_ONLINE_ONLINE_NODE_H

