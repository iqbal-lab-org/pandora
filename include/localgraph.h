#ifndef __LOCALGRAPH_H_INCLUDED__ // if localgraph.h hasn't been included yet...
#define __LOCALGRAPH_H_INCLUDED__

#include "IITree.h"
#include "interval.h"
#include "localnode.h"
#include "prg/path.h"
#include <cstdint>
#include <cstring>
#include <iostream>
#include <map>
#include <vector>

class LocalGraph {
public:
    std::map<uint32_t, LocalNodePtr> nodes; // representing nodes in graph
    IITree<uint32_t, LocalNodePtr> intervalTree; // TODO: move to private
    std::map<uint32_t, LocalNodePtr> startIndexOfZeroLengthIntervals;
    std::map<uint32_t, LocalNodePtr> startIndexOfAllIntervals;

    LocalGraph();

    ~LocalGraph();

    void add_node(const uint32_t& id, const std::string& seq, const Interval& pos);

    void add_edge(const uint32_t&, const uint32_t&);

    void write_gfa(const std::string&) const;

    void read_gfa(const std::string&);

    std::vector<PathPtr> walk(const uint32_t&, const uint32_t&, const uint32_t&) const;

    std::vector<PathPtr> walk_back(
        const uint32_t&, const uint32_t&, const uint32_t&) const;

    LocalNodePtr get_previous_node(const LocalNodePtr) const;

    std::vector<LocalNodePtr> nodes_along_string(
        const std::string&, bool end_to_end = false) const;

    std::vector<LocalNodePtr> top_path() const;

    std::vector<LocalNodePtr> bottom_path() const;

    bool operator==(const LocalGraph& y) const;

    bool operator!=(const LocalGraph& y) const;

    friend std::ostream& operator<<(std::ostream& out, LocalGraph const& data);
};

#endif
