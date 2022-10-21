#ifndef __PANREAD_H_INCLUDED__ // if panread.h hasn't been included yet...
#define __PANREAD_H_INCLUDED__

#include <string>
#include <cstdint>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <minihit.h>
#include "minihits.h"

class pangenome::Read {
private:
    // TODO: derive this from MinimizerHits?
    // TODO: or maybe keep it here but without the read id, since it is duplicated?
    std::vector<MinimizerHit*> hits; // store all Minimizer Hits mapping to this read
    std::vector<WeakNodePtr> nodes;

public:
    const uint32_t id; // read id
    std::vector<bool> node_orientations;

    // constructor/destructors
    Read(const uint32_t);
    virtual ~Read();

    // TODO: this can be a source of time inneficiency at the cost of using less memory
    // TODO: check if we should fallback to representing hits as
    // std::unordered_map<uint32_t, std::vector<MinimizerHitPtr>> directly
    std::unordered_map<uint32_t, std::vector<MinimizerHitPtr>>
    get_hits_as_unordered_map() const;

    const std::vector<WeakNodePtr>& get_nodes() const { return nodes; }
    // TODO: this getter allows the caller to the private attribute nodes, use with
    // care...
    // TODO: replace/remove this?
    std::vector<WeakNodePtr>& get_nodes() { return nodes; }

    // TODO: just used in tests, move to private and friend the test class
    void set_nodes(const std::vector<WeakNodePtr>& nodes) { this->nodes = nodes; }

    std::vector<WeakNodePtr>::iterator find_node_by_id(uint32_t node_id);

    // modifiers
    void add_node(const NodePtr& nodePtr)
    {
        nodes.push_back(nodePtr);
        nodes.shrink_to_fit();
    }
    void add_orientation(bool orientation)
    {
        node_orientations.push_back(orientation);
        node_orientations.shrink_to_fit();
    }

    void add_hits(
        const NodePtr& node_ptr, const MinimizerHits& cluster);

    std::pair<uint32_t, uint32_t> find_position(const std::vector<uint_least32_t>&,
        const std::vector<bool>&, const uint16_t min_overlap = 1);

    // remove nodes
    void remove_all_nodes_with_this_id(uint32_t node_id);
    std::vector<WeakNodePtr>::iterator remove_node_with_iterator(
        std::vector<WeakNodePtr>::iterator nit);

    // replace nodes
    void replace_node_with_iterator(
        std::vector<WeakNodePtr>::iterator n_original, NodePtr n);

    bool operator==(const Read& y) const;

    bool operator!=(const Read& y) const;

    bool operator<(const Read& y) const;

    friend std::ostream& operator<<(std::ostream& out, const Read& r);

    friend class pangenome::Graph;

    friend class pangenome::Node;
};

#endif
