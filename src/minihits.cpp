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


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

MinimizerHits::MinimizerHits(const uint32_t &num_hits) {
    uhits.reserve(num_hits);
}

void MinimizerHits::clear() {
    /*for (const auto &c: hits)
    {
        delete c;
    }*/
    hits.clear();

    /*for (const auto &c: uhits)
    {   
        delete c;
    }*/
    uhits.clear();
}

MinimizerHits::~MinimizerHits() {
    clear();
}

void MinimizerHits::add_hit(const uint32_t i, const Minimizer &minimizerFromRead, const MiniRecord &minimizerFromPRG) {
    MinimizerHitPtr mh(std::make_shared<MinimizerHit>(i, minimizerFromRead, minimizerFromPRG));
    uhits.insert(mh);
}

void MinimizerHits::sort() {
    if (hits.max_size() > uhits.size()) {
        hits.insert(uhits.begin(), uhits.end());
        assert(hits.size() == uhits.size() ||
               assert_msg("Expected uhits.size()=" << uhits.size() << " to equal hits.size()=" << hits.size()));
        uhits.clear();
    } else {
        std::cerr << "Could not create a set big enough for " << uhits.size() << " elements. The max size is "
                  << hits.max_size() << std::endl;
        std::exit(EXIT_FAILURE);
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
    if (lhs->get_prg_id() < rhs->get_prg_id()) { return true; }
    if (rhs->get_prg_id() < lhs->get_prg_id()) { return false; }
    //want those that match against the same prg_path together
    if (lhs->get_prg_path() < rhs->get_prg_path()) { return true; }
    if (rhs->get_prg_path() < lhs->get_prg_path()) { return false; }
    //separated into two categories, corresponding to a forward, and a rev-complement hit, note fwd come first
    if (lhs->is_forward() > rhs->is_forward()) { return true; }
    if (rhs->is_forward() > lhs->is_forward()) { return false; }
    // finally, make sure that hits from separate reads aren't removed from the set as "=="
    if (lhs->get_read_id() < rhs->get_read_id()) { return true; }
    if (rhs->get_read_id() < lhs->get_read_id()) { return false; }
    if (lhs->get_read_start_position() < rhs->get_read_start_position()) { return true; }
    if (rhs->get_read_start_position() < lhs->get_read_start_position()) { return false; }
    return false;
}

bool clusterComp::operator()(std::set<MinimizerHitPtr, pComp> lhs, std::set<MinimizerHitPtr, pComp> rhs) {
    if ((*lhs.begin())->get_read_id() < (*rhs.begin())->get_read_id()) { return true; }
    if ((*rhs.begin())->get_read_id() < (*lhs.begin())->get_read_id()) { return false; }
    if ((*lhs.begin())->get_read_start_position() < (*rhs.begin())->get_read_start_position()) { return true; }
    if ((*rhs.begin())->get_read_start_position() < (*lhs.begin())->get_read_start_position()) { return false; }
    if (lhs.size() > rhs.size()) { return true; } // want bigger first!
    if (rhs.size() > lhs.size()) { return false; }
    if ((*lhs.begin())->get_prg_id() < (*rhs.begin())->get_prg_id()) { return true; }
    if ((*rhs.begin())->get_prg_id() < (*lhs.begin())->get_prg_id()) { return false; }
    if ((*lhs.begin())->get_prg_path() < (*rhs.begin())->get_prg_path()) { return true; }
    if ((*rhs.begin())->get_prg_path() < (*lhs.begin())->get_prg_path()) { return false; }
    if ((*lhs.begin())->is_forward() < (*rhs.begin())->is_forward()) { return true; }
    if ((*rhs.begin())->is_forward() < (*lhs.begin())->is_forward()) { return false; }
    return false;
}

bool clusterComp_size::operator()(std::set<MinimizerHitPtr, pComp> lhs, std::set<MinimizerHitPtr, pComp> rhs) {
    if ((*lhs.begin())->get_read_id() < (*rhs.begin())->get_read_id()) { return true; }
    if ((*rhs.begin())->get_read_id() < (*lhs.begin())->get_read_id()) { return false; }
    if (lhs.size() > rhs.size()) { return true; }
    if (rhs.size() > lhs.size()) { return false; }
    if ((*lhs.begin())->get_read_start_position() < (*rhs.begin())->get_read_start_position()) { return true; }
    if ((*rhs.begin())->get_read_start_position() < (*lhs.begin())->get_read_start_position()) { return false; }
    if ((*lhs.begin())->get_prg_id() < (*rhs.begin())->get_prg_id()) { return true; }
    if ((*rhs.begin())->get_prg_id() < (*lhs.begin())->get_prg_id()) { return false; }
    if ((*lhs.begin())->get_prg_path() < (*rhs.begin())->get_prg_path()) { return true; }
    if ((*rhs.begin())->get_prg_path() < (*lhs.begin())->get_prg_path()) { return false; }
    if ((*lhs.begin())->is_forward() < (*rhs.begin())->is_forward()) { return true; }
    if ((*rhs.begin())->is_forward() < (*lhs.begin())->is_forward()) { return false; }
    return false;
}
