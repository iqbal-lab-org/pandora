#ifndef __PANSAMPLE_H_INCLUDED__   // if pansample.h hasn't been included yet...
#define __PANSAMPLE_H_INCLUDED__

#include <string>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include "pangenome/ns.cpp"

class KmerNode;

typedef std::shared_ptr<KmerNode> KmerNodePtr;

class pangenome::Sample {
public:
    const std::string name; // first column in index of read files
    vector<NodePtr> nodes;
    vector<bool> node_orientations;
    std::unordered_map<uint32_t, std::vector<std::vector<KmerNodePtr>>> paths; // from prg id (or unique id) to kmernnode path(s) through each node

    Sample(const std::string &);

    void add_path(const uint32_t, const std::vector<KmerNodePtr> &);

    bool operator==(const Sample &y) const;

    bool operator!=(const Sample &y) const;

    bool operator<(const Sample &y) const;

    friend std::ostream &operator<<(std::ostream &out, const Sample &r);
};

#endif
