#ifndef __EXTRACT_READS_H_INCLUDED__
#define __EXTRACT_READS_H_INCLUDED__

#include <set>
#include <algorithm>
#include <memory>
#include <cassert>
#include <utility>
#include <vector>
#include "localnode.h"
#include "minihits.h"
#include "denovo_discovery/local_assembly.h"
#include "interval.h"
#include "minihit.h"
#include "pangenome/ns.cpp"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "localPRG.h"
#include "fastaq_handler.h"
#include "fastaq.h"
#include "utils.h"
#include "prg/path.h"


using PanNodePtr = std::shared_ptr<pangenome::Node>;


struct ReadCoordinate {
    uint32_t id;
    uint32_t start;
    uint32_t end;
    bool is_forward;

    ReadCoordinate(uint32_t id, uint32_t start, uint32_t end, bool is_forward);

    bool operator<(const ReadCoordinate &y) const;

    bool operator==(const ReadCoordinate &y) const;

    bool operator!=(const ReadCoordinate &y) const;

    friend std::ostream &operator<<(std::ostream &, ReadCoordinate const &);
};


template<>
struct std::hash<ReadCoordinate> {
    size_t operator()(const ReadCoordinate &coordinate) const;
};


struct PathComponents {
    prg::Path flank_left;
    prg::Path slice;
    prg::Path flank_right;

    PathComponents();

    PathComponents(prg::Path flank_left, prg::Path slice, prg::Path flank_right);

    bool operator==(const PathComponents &other) const;

    bool operator!=(const PathComponents &other) const;
};


PathComponents find_interval_and_flanks_in_localpath(const Interval &interval,
                                                     const std::vector<LocalNodePtr> &local_node_max_likelihood_path);

std::vector<MinimizerHitPtr>
find_hits_inside_path(const std::vector<MinimizerHitPtr> &read_hits, const prg::Path &local_path);

#endif
