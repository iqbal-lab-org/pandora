#ifndef __PANGRAPH_H_INCLUDED__   // if pangraph.h hasn't been included yet...
#define __PANGRAPH_H_INCLUDED__

class PanNode;
class PanEdge;
class PanRead;
class PanSample;
struct MinimizerHit;
class LocalPRG;
class KmerNode;

#include <cstring>
#include <map>
#include <ostream>
#include <functional>
#include <minihits.h>

class PanGraph {
    uint next_id;
    std::vector<PanEdge *> edges;
    std::map<uint32_t, PanRead *> reads;
    std::map<std::string, PanSample *> samples;
public:

    std::map<uint32_t, PanNode *> nodes;

    PanGraph();
    ~PanGraph();
    void clear();

    // graph additions/removals
    void add_node(const uint32_t, const std::string, uint32_t,
                  const std::set<MinimizerHit *, pComp> &); // used by pandora map
    void add_node(const uint32_t, const std::string &, const std::string &, const std::vector<KmerNode *> &,
                  const LocalPRG *); // used by pandora compare
    PanEdge *add_edge(const uint32_t &, const uint32_t &, const uint &);
    void add_edge(const uint32_t &, const uint32_t &, const uint &, const uint &);
    std::vector<PanEdge *>::iterator remove_edge(PanEdge *);
    std::map<uint32_t, PanNode *>::iterator remove_node(PanNode *);
    std::vector<PanEdge *>::iterator add_shortcut_edge(const std::vector<PanEdge *>::iterator, PanRead *);

    // graph manipulation
    std::vector<PanEdge *>::iterator split_node_by_edges(PanNode *, PanEdge *, PanEdge *);
    void split_nodes_by_reads(const uint &, const uint &);

    // graph cleaning
    void read_clean(const uint &);
    void remove_low_covg_nodes(const uint &);
    void remove_low_covg_edges(const uint &);
    void clean(const uint32_t &);

    // graph magic
    void add_hits_to_kmergraphs(const std::vector<LocalPRG *> &);

    // graph comparison
    bool operator==(const PanGraph &y) const;
    bool operator!=(const PanGraph &y) const;

    // graph read/write
    void write_gfa(const std::string &);
    void save_matrix(const std::string &);

    friend std::ostream &operator<<(std::ostream &out, const PanGraph &m);

    friend class PanGraphTest_add_node_Test;
    friend class PanGraphTest_add_node_sample_Test;
    friend class PanGraphTest_add_edge_Test;
    friend class PanGraphTest_add_shortcut_edge_Test;
    friend class PanGraphTest_read_clean_Test;
    friend class PanGraphTest_remove_edge_Test;
    friend class PanGraphTest_split_node_by_edges_Test;
    friend class PanReadTest_replace_edge1_Test;
    friend class PanReadTest_replace_edge2_Test;
    friend class PanReadTest_remove_edge1_Test;
    friend class PanReadTest_remove_edge2_Test;
    friend class PanReadTest_replace_node_Test;
    friend class PanReadTest_remove_node_Test;
};

#endif

