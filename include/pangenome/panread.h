#ifndef __PANREAD_H_INCLUDED__   // if panread.h hasn't been included yet...
#define __PANREAD_H_INCLUDED__

#include <string>
#include <cstdint>
#include <vector>
#include <iostream>
#include <unordered_map>
#include "minihits.h"
#include "pangenome/ns.cpp"


class pangenome::Read {
private:
    //derive this from MinimizerHits?
    //or maybe keep it here but without the read id, since it is duplicated?
    //todo: there are things that can be done here
    //todo: maybe use sdsl?
    std::vector<std::pair<uint32_t, std::vector<MinimizerHit*>*>> hits; // from node id to cluster of hits against that node in this read

public:
    const uint32_t id; // corresponding the the read id
    std::vector<bool> node_orientations;
    std::vector<NodePtr> nodes;

    std::unordered_map<uint32_t, std::vector<MinimizerHitPtr>> getHits() const {
        std::unordered_map<uint32_t, std::vector<MinimizerHitPtr>> hitsMap;
        for (const auto &pairNodeIdAndClusterHits : hits) {
            std::vector<MinimizerHitPtr> minimizerHits;
            for (MinimizerHit* rawPointer : *(pairNodeIdAndClusterHits.second))
                minimizerHits.push_back(std::make_shared<MinimizerHit>(*rawPointer));
            hitsMap[pairNodeIdAndClusterHits.first] = minimizerHits;
        }

        return hitsMap;
    }




    Read(const uint32_t);
    virtual ~Read();

    void add_hits(const uint32_t, std::set<MinimizerHitPtr, pComp> &);

    std::pair<uint32_t, uint32_t>
    find_position(const std::vector<uint_least32_t> &, const std::vector<bool> &, const uint16_t min_overlap = 1);

    void remove_node(NodePtr);

    std::vector<NodePtr>::iterator remove_node(std::vector<NodePtr>::iterator);

    // remove nodes
    void replace_node(std::vector<NodePtr>::iterator, NodePtr);
    // replace nodes

    bool operator==(const Read &y) const;

    bool operator!=(const Read &y) const;

    bool operator<(const Read &y) const;

    friend std::ostream &operator<<(std::ostream &out, const Read &r);

    friend class pangenome::Graph;

    friend class pangenome::Node;
};

#endif
