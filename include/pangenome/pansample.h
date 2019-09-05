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
    const uint32_t sample_id;
    std::vector<WeakNodePtr> nodes;
    std::vector<bool> node_orientations;
    std::unordered_map<uint32_t, std::vector<std::vector<KmerNodePtr>>> paths; // from prg id (or unique id) to kmernnode path(s) through each node

    Sample(const std::string &, const uint32_t &id);

    void add_path(const uint32_t, const std::vector<KmerNodePtr> &);

    bool operator==(const Sample &y) const;

    bool operator!=(const Sample &y) const;

    bool operator<(const Sample &y) const;

    friend std::ostream &operator<<(std::ostream &out, const Sample &s);
};

struct pangenome::SamplePtrSorterBySampleId {
    bool operator()(const pangenome::SamplePtr &lhs, const pangenome::SamplePtr &rhs) const;
};

#endif
