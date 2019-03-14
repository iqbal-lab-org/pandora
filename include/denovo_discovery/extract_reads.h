#ifndef __EXTRACT_READS_H_INCLUDED__   // if extract_reads.h hasn't been included yet...
#define __EXTRACT_READS_H_INCLUDED__

#include <set>
#include <algorithm>
#include <memory>
#include "interval.h"
#include "localnode.h"
#include "minihits.h"
#include "pangenome/ns.cpp"

#include "denovo_discovery/local_assembly.h"
#include "gene_interval_info.h"

using PanNodePtr = std::shared_ptr<pangenome::Node>;

struct ReadCoordinate {
  uint32_t id;
  uint32_t start;
  uint32_t end;
  bool strand;

  bool operator<(const ReadCoordinate &y) const;

  bool operator==(const ReadCoordinate &y) const;

  bool operator!=(const ReadCoordinate &y) const;

  friend std::ostream &operator<<(std::ostream &, ReadCoordinate const &);
};

struct PathComponents {
  prg::Path flank_left;
  prg::Path slice;
  prg::Path flank_right;
};

std::vector<Interval>
identify_regions(const std::vector<uint32_t> &covgs, const uint32_t &threshold, const uint32_t &min_length,
                 const uint_least16_t &padding = 0);

PathComponents
find_interval_in_localpath(const Interval &interval, const std::vector<LocalNodePtr> &lmp);

std::set<MinimizerHitPtr, pComp_path>
hits_inside_path(const std::set<MinimizerHitPtr, pComp_path> &read_hits, const prg::Path &local_path);

std::set<ReadCoordinate>
get_read_overlap_coordinates(const PanNodePtr &, const prg::Path &local_path, const uint32_t &min_number_hits = 2);

using ReadPileup = std::vector<std::string>;

namespace denovo_discovery {
void add_pnode_coordinate_pairs(std::vector<std::shared_ptr<LocalPRG>> &prgs,
                                std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
                                const PanNodePtr &pnode,
                                const std::vector<LocalNodePtr> &local_node_path,
                                const std::vector<KmerNodePtr> &kmer_node_path,
                                const uint_least16_t &padding_size = 0,
                                const uint32_t &low_coverage_threshold = 2,
                                const uint32_t &interval_min_length = 5,
                                const uint32_t &min_number_hits = 2);

std::map<GeneIntervalInfo, ReadPileup>
collect_read_pileups(const std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &, const boost::filesystem::path &);
}

#endif
