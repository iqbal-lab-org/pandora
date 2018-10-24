#ifndef __GENE_INTERVAL_INFO_H_INCLUDED__   // if gene_interval_info.h hasn't been included yet...
#define __GENE_INTERVAL_INFO_H_INCLUDED__

#include <ostream>
#include <memory>
#include "pangenome/pannode.h"
#include "interval.h"


typedef std::shared_ptr<pangenome::Node> PanNodePtr;

struct GeneIntervalInfo {
    PanNodePtr pnode;
    Interval interval;
    std::string seq;

    bool operator<(const GeneIntervalInfo &y) const;

    bool operator==(const GeneIntervalInfo &y) const;

    bool operator!=(const GeneIntervalInfo &y) const;

    friend std::ostream &operator<<(std::ostream &out, const GeneIntervalInfo &y);
};

#endif
