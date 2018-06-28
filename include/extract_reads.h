#ifndef __EXTRACT_READS_H_INCLUDED__   // if extract_reads.h hasn't been included yet...
#define __EXTRACT_READS_H_INCLUDED__

#include "interval.h"

std::vector<Interval> identify_regions(const std::vector<uint32_t>&, const uint32_t& threshold=0, const uint32_t& min_length=0);

std::vector<LocalNodePtr> find_interval_in_localpath(const Interval&, const std::vector<LocalNodePtr>&);

std::set<MinimizerHitPtr, pComp_path>& hits_along_path(const std::set<MinimizerHitPtr, pComp_path>&,
                                                       const std::vector<LocalNodePtr>&);

void get_read_overlap_coordinates(PanNodePtr, std::vector<std::vector<uint32_t>>&, std::vector<LocalNodePtr>&);


#endif
