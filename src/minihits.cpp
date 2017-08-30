#include <cassert>
#include <functional>
#include <iostream>
#include <sstream>
#include <cstring>
#include "minihits.h"
#include "minihit.h"
#include "minirecord.h"
#include "minimizer.h"
#include "utils.h" // for pointer_values_equal

using namespace std;

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

MinimizerHits::MinimizerHits(const uint& num_hits) {
    uhits.reserve(num_hits);
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
    //unordered_set<MinimizerHit*, Hash, pEq>::iterator it=uhits.find(mh);
    //if(it==uhits.end())
    //pointer_values_equal<MinimizerHit> eq = { mh };
    //if (find_if(uhits.begin(), uhits.end(), eq) == uhits.end())
    //{
    uhits.insert(mh);
	//cout << "added hit " << *mh << endl;
    //} else {
    //    delete mh;
    //}
}

void MinimizerHits::sort()
{
    if (hits.max_size()>uhits.size())
    {
        hits.insert(uhits.begin(), uhits.end());
        /*for (auto h=uhits.begin(); h!=uhits.end(); ++h)
	{
	    auto f = hits.find(*h);
	    if (f == hits.end())
	    {
		hits.insert(*h);
	    } else {
		cout << "Minihit " << **h << " was already in the set: " << **f << endl;
	    }
	}*/
        assert(hits.size() == uhits.size() || assert_msg("Expected uhits.size()=" << uhits.size() << " to equal hits.size()=" << hits.size()));
        uhits.clear();
    } else {
	cerr << "Could not create a set big enough for " << uhits.size() << " elements. The max size is " << hits.max_size() << endl;
	exit (EXIT_FAILURE);
    }
    return;
}

/*std::ostream& operator<< (std::ostream & out, MinimizerHits const& m) {
    out << "(" << m.read_id << ", " << m.read_interval << ", " << m.prg_id << ", " << m.prg_path << ", " << strand << ")";
    return out ;
}*/

bool pComp::operator () (MinimizerHit* lhs, MinimizerHit* rhs) {
        return (*lhs)<(*rhs);
}

bool pEq::operator() (const MinimizerHit* lhs, const MinimizerHit* rhs) const {
        return *lhs == *rhs;
}

size_t Hash::operator()(const MinimizerHit* mh) const
{
    stringstream ss;
    string temp;
    ss << *mh;
    ss >> temp;
    return hash<string>()(temp);
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
    if ((*lhs.begin())->read_id < (*rhs.begin())->read_id) { return true; }
    if ((*rhs.begin())->read_id < (*lhs.begin())->read_id) { return false; }
    if ((*lhs.begin())->read_interval.start < (*rhs.begin())->read_interval.start) { return true; }
    if ((*rhs.begin())->read_interval.start < (*lhs.begin())->read_interval.start) { return false; }
    return false;
}
