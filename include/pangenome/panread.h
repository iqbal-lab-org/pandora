#ifndef __PANREAD_H_INCLUDED__   // if panread.h hasn't been included yet...
#define __PANREAD_H_INCLUDED__

#include <string>
#include <cstdint>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <minihit.h>
#include "minihits.h"
#include "pangenome/ns.cpp"


class pangenome::Read {
private:
    //derive this from MinimizerHits?
    //or maybe keep it here but without the read id, since it is duplicated?
    //todo: there are things that can be done here
    //todo: maybe use sdsl?
    std::vector<MinimizerHit*> hits; //store all Minimizer Hits mapping to this read
    std::vector<WeakNodePtr> nodes;

public:
    const uint32_t id; // corresponding the the read id
    std::vector<bool> node_orientations;

    //constructor/destructors
    Read(const uint32_t);
    virtual ~Read();


    //getters
    std::unordered_map<uint32_t, std::vector<MinimizerHitPtr>> getHits() const {
        std::unordered_map<uint32_t, std::vector<MinimizerHitPtr>> hitsMap; //this will map node_ids from the pangenome::Graph to their minimizer hits
        for (const MinimizerHit * const minihit : hits) {
            //gets the nodeId
            uint32_t nodeId = minihit->get_prg_id(); //prg_id == node_id in pangenome::Graph

            //checks if we have an entry for this nodeId in hitsMap
            if (hitsMap.find(nodeId) == hitsMap.end())
                //no, add it
                hitsMap[nodeId] = std::vector<MinimizerHitPtr>();

            //add this minihit to hitsMap
            hitsMap[nodeId].push_back(std::make_shared<MinimizerHit>(*minihit)); //TODO: I think here we don't really need to create a shared pointer - a raw pointer is fine
        }

        return hitsMap;
    }
    const std::vector<WeakNodePtr>& get_nodes() const {
        return nodes;
    }

    //TODO: this can modify nodes, use with care...
    //TODO: replace this?
    std::vector<WeakNodePtr>& get_nodes() {
        return nodes;
    }

    std::vector<WeakNodePtr>::iterator find_node_by_id (uint32_t node_id);

    //modifiers
    void add_node(const NodePtr &nodePtr) {
        nodes.push_back(nodePtr);
        nodes.shrink_to_fit();
    }
    void add_orientation(bool orientation) {
        node_orientations.push_back(orientation);
        node_orientations.shrink_to_fit();
    }


    void add_hits(const uint32_t, std::set<MinimizerHitPtr, pComp> &);

    std::pair<uint32_t, uint32_t>
    find_position(const std::vector<uint_least32_t> &, const std::vector<bool> &, const uint16_t min_overlap = 1);

    // remove nodes
    void remove_all_nodes_with_this_id(uint32_t node_id);
    std::vector<WeakNodePtr>::iterator remove_node_with_iterator(std::vector<WeakNodePtr>::iterator nit);

    // replace nodes
    void replace_node_with_iterator(std::vector<WeakNodePtr>::iterator n_original, NodePtr n);


    bool operator==(const Read &y) const;

    bool operator!=(const Read &y) const;

    bool operator<(const Read &y) const;

    friend std::ostream &operator<<(std::ostream &out, const Read &r);

    friend class pangenome::Graph;

    friend class pangenome::Node;
};

#endif
