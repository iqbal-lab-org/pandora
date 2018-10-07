#ifndef __EXTRACT_READS_H_INCLUDED__   // if extract_reads.h hasn't been included yet...
#define __EXTRACT_READS_H_INCLUDED__

#include <set>
#include <algorithm>
#include <memory>
#include "interval.h"
#include "localnode.h"
#include "minihits.h"
#include "pangenome/ns.cpp"
#include "local_assembly.h"

typedef std::shared_ptr<pangenome::Node> PanNodePtr;

std::vector<Interval>
identify_regions(const std::vector<uint32_t> &, const uint32_t &threshold = 0, const uint32_t &min_length = 0);

vector<LocalNodePtr>
find_interval_in_localpath(const Interval &, const vector<LocalNodePtr> &, const unsigned int &buff);

std::set<MinimizerHitPtr, pComp_path> hits_along_path(const std::set<MinimizerHitPtr, pComp_path> &,
                                                      const std::vector<LocalNodePtr> &);

void get_read_overlap_coordinates(PanNodePtr, std::set<std::vector<uint32_t>> &, std::vector<LocalNodePtr> &);

void save_read_strings_to_denovo_assemble(const std::string &,
                                          const std::string &,
                                          const PanNodePtr,
                                          const std::vector<LocalNodePtr> &,
                                          const std::vector<KmerNodePtr> &,
                                          const unsigned int &buff = 0,
                                          const uint32_t &threshold = 2,
                                          const uint32_t &min_length = 5
);


Interval apply_buffer_to_interval(const Interval &interval, const int32_t &buff);

#endif
