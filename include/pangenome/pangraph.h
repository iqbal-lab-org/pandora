#ifndef __PANGRAPH_H_INCLUDED__   // if pangraph.h hasn't been included yet...
#define __PANGRAPH_H_INCLUDED__

struct MinimizerHit;
class LocalPRG;
class KmerNode;

#include <cstring>
#include <unordered_map>
#include <ostream>
#include <functional>
#include <memory>
#include "minihits.h"
#include "pangenome/ns.cpp"


typedef std::shared_ptr<KmerNode> KmerNodePtr;

class pangenome::Graph {
protected:
    std::unordered_map<std::string, SamplePtr> samples;
public:
    std::unordered_map<uint32_t, ReadPtr> reads;
    std::unordered_map<uint32_t, NodePtr> nodes;

    Graph();
    ~Graph();
    void clear();

    // graph additions/removals
    void add_node(const uint32_t, const std::string, uint32_t,
                  const std::set<MinimizerHitPtr, pComp> &); // used by pandora map
    void add_node(const uint32_t, const std::string &, const std::string &, const std::vector<KmerNodePtr> &,
                  const LocalPRG *); // used by pandora compare

    std::unordered_map<uint32_t, NodePtr>::iterator remove_node(NodePtr);
    void remove_read(const uint32_t);

    void remove_low_covg_nodes(const uint &);

    void add_hits_to_kmergraphs(const std::vector<LocalPRG *> &);

    // graph comparison
    bool operator==(const Graph &y) const;
    bool operator!=(const Graph &y) const;

    // graph read/write
    void save_matrix(const std::string &);

    friend std::ostream &operator<<(std::ostream &out, const Graph &m);

};

#endif

