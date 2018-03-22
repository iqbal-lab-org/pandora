#include <cassert>
#include <functional>
#include <iostream>
#include <sstream>
#include <memory>
#include "minihits.h"
#include "minihit.h"
#include "minirecord.h"
#include "minimizer.h"
#include "utils.h" // for pointer_values_equal

using namespace std;

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

MinimizerHits::MinimizerHits(const uint &num_hits) {
    uhits.reserve(num_hits);
}

void MinimizerHits::clear() {
    /*for (auto c: hits)
    {
        delete c;
    }*/
    hits.clear();

    /*for (auto c: uhits)
    {   
        delete c;
    }*/
    uhits.clear();
}

MinimizerHits::~MinimizerHits() {
    clear();
}

void MinimizerHits::add_hit(const uint32_t i, const Minimizer& m, const MiniRecord *r) {
    MinimizerHitPtr mh(make_shared<MinimizerHit>(i, m, r));
    uhits.insert(mh);
}

void MinimizerHits::sort() {
    if (hits.max_size() > uhits.size()) {
        hits.insert(uhits.begin(), uhits.end());
        assert(hits.size() == uhits.size() ||
               assert_msg("Expected uhits.size()=" << uhits.size() << " to equal hits.size()=" << hits.size()));
        uhits.clear();
    } else {
        cerr << "Could not create a set big enough for " << uhits.size() << " elements. The max size is "
             << hits.max_size() << endl;
        exit(EXIT_FAILURE);
    }
}

/*std::ostream& operator<< (std::ostream & out, MinimizerHits const& m) {
    out << "(" << m.read_id << ", " << m.read_start_position << ", " << m.prg_id << ", " << m.prg_path << ", " << strand << ")";
    return out ;
}*/

bool pComp::operator()(const MinimizerHitPtr &lhs, const MinimizerHitPtr &rhs) {
    return (*lhs) < (*rhs);
}

bool pEq::operator()(const MinimizerHitPtr &lhs, const MinimizerHitPtr &rhs) const {
    return *lhs == *rhs;
}

/*size_t Hash::operator()(const MinimizerHit* mh) const
{
    stringstream ss;
    string temp;
    ss << *mh;
    ss >> temp;
    return hash<string>()(temp);
}*/

bool pComp_path::operator()(const MinimizerHitPtr &lhs, const MinimizerHitPtr &rhs) {
    // should be same id
    if (lhs->prg_id < rhs->prg_id) { return true; }
    if (rhs->prg_id < lhs->prg_id) { return false; }
    //want those that match against the same prg_path together
    if (lhs->prg_path < rhs->prg_path) { return true; }
    if (rhs->prg_path < lhs->prg_path) { return false; }
    //separated into two categories, corresponding to a forward, and a rev-complement hit, note fwd come first
    if (lhs->strand > rhs->strand) { return true; }
    if (rhs->strand > lhs->strand) { return false; }
    // finally, make sure that hits from separate reads aren't removed from the set as "=="
    if (lhs->read_id < rhs->read_id) { return true; }
    if (rhs->read_id < lhs->read_id) { return false; }
    if (lhs->read_start_position < rhs->read_start_position) { return true; }
    if (rhs->read_start_position < lhs->read_start_position) { return false; }
    return false;
}

bool clusterComp::operator()(set<MinimizerHitPtr, pComp> lhs, set<MinimizerHitPtr, pComp> rhs) {
    if ((*lhs.begin())->read_id < (*rhs.begin())->read_id) { return true; }
    if ((*rhs.begin())->read_id < (*lhs.begin())->read_id) { return false; }
    if ((*lhs.begin())->read_start_position < (*rhs.begin())->read_start_position) { return true; }
    if ((*rhs.begin())->read_start_position < (*lhs.begin())->read_start_position) { return false; }
    if (lhs.size() > rhs.size()) { return true; } // want bigger first!
    if (rhs.size() > lhs.size()) { return false; }
    if ((*lhs.begin())->prg_id < (*rhs.begin())->prg_id) { return true; }
    if ((*rhs.begin())->prg_id < (*lhs.begin())->prg_id) { return false; }
    if ((*lhs.begin())->prg_path < (*rhs.begin())->prg_path) { return true; }
    if ((*rhs.begin())->prg_path < (*lhs.begin())->prg_path) { return false; }
    if ((*lhs.begin())->strand < (*rhs.begin())->strand) { return true; }
    if ((*rhs.begin())->strand < (*lhs.begin())->strand) { return false; }
    return false;
}

bool clusterComp_size::operator()(set<MinimizerHitPtr, pComp> lhs, set<MinimizerHitPtr, pComp> rhs) {
    if ((*lhs.begin())->read_id < (*rhs.begin())->read_id) { return true; }
    if ((*rhs.begin())->read_id < (*lhs.begin())->read_id) { return false; }
    if (lhs.size() > rhs.size()) { return true; }
    if (rhs.size() > lhs.size()) { return false; }
    if ((*lhs.begin())->read_start_position < (*rhs.begin())->read_start_position) { return true; }
    if ((*rhs.begin())->read_start_position < (*lhs.begin())->read_start_position) { return false; }
    if ((*lhs.begin())->prg_id < (*rhs.begin())->prg_id) { return true; }
    if ((*rhs.begin())->prg_id < (*lhs.begin())->prg_id) { return false; }
    if ((*lhs.begin())->prg_path < (*rhs.begin())->prg_path) { return true; }
    if ((*rhs.begin())->prg_path < (*lhs.begin())->prg_path) { return false; }
    if ((*lhs.begin())->strand < (*rhs.begin())->strand) { return true; }
    if ((*rhs.begin())->strand < (*lhs.begin())->strand) { return false; }
    return false;
}
