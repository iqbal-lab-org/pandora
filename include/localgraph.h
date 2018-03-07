#ifndef __LOCALGRAPH_H_INCLUDED__   // if localgraph.h hasn't been included yet...
#define __LOCALGRAPH_H_INCLUDED__

class LocalNode;

#include <cstring>
#include <map>
#include <vector>
#include <iostream>
#include "interval.h"
#include "path.h"
#include "localnode.h"

class LocalGraph {
public:
    std::map<uint32_t, LocalNodePtr> nodes; // representing nodes in graph

    LocalGraph();

    ~LocalGraph();

    void add_node(const uint32_t &id, const std::string &seq, const Interval &pos);

    void add_edge(const uint32_t &, const uint32_t &);

    void write_gfa(const std::string &) const;

    void read_gfa(const std::string &);

    std::vector<Path> walk(const uint32_t &, const uint32_t &, const uint32_t &) const;

    std::vector<Path> walk_back(const uint32_t &, const uint32_t &, const uint32_t &) const;

    LocalNodePtr get_previous_node(const LocalNodePtr) const;

    std::vector<LocalNodePtr> nodes_along_string(const std::string &) const;

    std::vector<LocalNodePtr> top_path() const;

    std::vector<LocalNodePtr> bottom_path() const;

    bool operator==(const LocalGraph &y) const;

    bool operator!=(const LocalGraph &y) const;

    friend std::ostream &operator<<(std::ostream &out, LocalGraph const &data);
};

#endif
