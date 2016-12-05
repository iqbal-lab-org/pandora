#include <cassert>
#include <functional>
#include <iostream>
#include <cstring>
#include "minihits.h"
#include "minihit.h"
#include "minirecord.h"
#include "minimizer.h"

using namespace std;

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
    set<MinimizerHit*, pComp>::iterator it=hits.find(mh);
    if(it==hits.end())
    {
        hits.insert(mh);
        //cout << "Added hit " << *mh << endl;
        //cout << *mh << endl;
    } else {
        delete mh;
    }
}

/*std::ostream& operator<< (std::ostream & out, MinimizerHits const& m) {
    out << "(" << m.read_id << ", " << m.read_interval << ", " << m.prg_id << ", " << m.prg_path << ", " << strand << ")";
    return out ;
}*/

bool pComp::operator()(MinimizerHit* lhs, MinimizerHit* rhs) {
        return (*lhs)<(*rhs);
}

bool pComp_path::operator()(MinimizerHit* lhs, MinimizerHit* rhs) {
	if (!(lhs->prg_path==rhs->prg_path))
	{
            return lhs->prg_path<rhs->prg_path;
	} else {
	    return lhs->read_id<rhs->read_id;
	}
}

bool clusterComp::operator()(set<MinimizerHit*, pComp> lhs, set<MinimizerHit*, pComp> rhs) {
    if ((*lhs.begin())->read_interval.start < (*rhs.begin())->read_interval.start) { return true; }
    if ((*rhs.begin())->read_interval.start < (*lhs.begin())->read_interval.start) { return false; }
    return false;
}
