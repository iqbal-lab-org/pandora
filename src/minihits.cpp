#include <cassert>
#include <functional>
#include <iostream>
#include <cstring>
#include "minihits.h"
#include "minihit.h"
#include "minirecord.h"
#include "minimizer.h"

using namespace std;

MinimizerHits::MinimizerHits() {
    uhits.reserve(20000);
}

MinimizerHits::~MinimizerHits()
{
    for (auto c: hits)
    {
        delete c;
    }
}

void MinimizerHits::add_hit(const uint32_t i, const Minimizer* m, const MiniRecord* r)
{
    MinimizerHit *mh;
    mh = new MinimizerHit(i, m, r);
    //set<MinimizerHit*, pComp>::iterator it=hits.find(mh);
    unordered_set<MinimizerHit*>::iterator it=uhits.find(mh);
    if(it==uhits.end())
    {
        uhits.insert(mh);
    } else {
        delete mh;
    }
}

void MinimizerHits::sort()
{
    hits.insert(uhits.begin(), uhits.end());
    uhits.clear();
    return;
}

/*std::ostream& operator<< (std::ostream & out, MinimizerHits const& m) {
    out << "(" << m.read_id << ", " << m.read_interval << ", " << m.prg_id << ", " << m.prg_path << ", " << strand << ")";
    return out ;
}*/

bool pComp::operator()(MinimizerHit* lhs, MinimizerHit* rhs) {
        return (*lhs)<(*rhs);
}

bool pComp_path::operator()(MinimizerHit* lhs, MinimizerHit* rhs) {
    //want those that match against the same prg_path together
    if (lhs->prg_path<rhs->prg_path) { return true;}
    if (rhs->prg_path<lhs->prg_path) { return false;}
    //separated into two categories, corresponding to a forward, and a rev-complement hit, note fwd come first
    if (lhs->strand>rhs->strand) {return true;}
    if (rhs->strand>lhs->strand) {return false;}
    // finally, make sure that hits from separate reads aren't removed from the set as "=="
    if (lhs->read_id<rhs->read_id) { return true;}
    if (rhs->read_id<lhs->read_id) { return false;}
    if (lhs->read_interval<rhs->read_interval) { return true;}
    if (rhs->read_interval<lhs->read_interval) { return false;}
    return false;
}

bool clusterComp::operator()(set<MinimizerHit*, pComp> lhs, set<MinimizerHit*, pComp> rhs) {
    if ((*lhs.begin())->read_interval.start < (*rhs.begin())->read_interval.start) { return true; }
    if ((*rhs.begin())->read_interval.start < (*lhs.begin())->read_interval.start) { return false; }
    return false;
}
