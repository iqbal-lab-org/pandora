#include <functional>
#include "minihits.h"
#include "minihit.h"

void MinimizerHits::insert(const uint32_t i, const Minimizer& minimizer_from_read,
    const MiniRecord& minimizer_from_PRG) {
    MinimizerHitPtr mh(
        std::make_shared<MinimizerHit>(i, minimizer_from_read, minimizer_from_PRG));
    this->insert(mh);
}

bool pComp::operator()(const MinimizerHitPtr& lhs, const MinimizerHitPtr& rhs) const
{
    return (*lhs) < (*rhs);
}

bool pEq::operator()(const MinimizerHitPtr& lhs, const MinimizerHitPtr& rhs) const
{
    return *lhs == *rhs;
}

bool pCompReadPositionFirst::operator()(const MinimizerHitPtr& lhs, const MinimizerHitPtr& rhs) const {
    if (lhs->get_read_id() < rhs->get_read_id()) {
        return true;
    }
    if (rhs->get_read_id() < lhs->get_read_id()) {
        return false;
    }
    if (lhs->get_read_start_position() < rhs->get_read_start_position()) {
        return true;
    }
    if (rhs->get_read_start_position() < lhs->get_read_start_position()) {
        return false;
    }
    if (lhs->same_strands() > rhs->same_strands()) {
        return true;
    }
    if (rhs->same_strands() > lhs->same_strands()) {
        return false;
    }
    if (lhs->get_prg_id() < rhs->get_prg_id()) {
        return true;
    }
    if (rhs->get_prg_id() < lhs->get_prg_id()) {
        return false;
    }
    if (lhs->get_prg_path() < rhs->get_prg_path()) {
        return true;
    }
    if (rhs->get_prg_path() < lhs->get_prg_path()) {
        return false;
    }
    return false;
}

bool MinimizerHits::operator<(const MinimizerHits &rhs) const
{
    auto first_hit_from_left = *this->begin();
    auto first_hit_from_right = *rhs.begin();

    if (first_hit_from_left->get_read_id() < first_hit_from_right->get_read_id()) {
        return true;
    }
    if (first_hit_from_right->get_read_id() < first_hit_from_left->get_read_id()) {
        return false;
    }
    if (first_hit_from_left->get_read_start_position()
        < first_hit_from_right->get_read_start_position()) {
        return true;
    }
    if (first_hit_from_right->get_read_start_position()
        < first_hit_from_left->get_read_start_position()) {
        return false;
    }
    if (this->size() > rhs.size()) {
        return true;
    } // want bigger first!
    if (rhs.size() > this->size()) {
        return false;
    }
    if (first_hit_from_left->get_prg_id() < first_hit_from_right->get_prg_id()) {
        return true;
    }
    if (first_hit_from_right->get_prg_id() < first_hit_from_left->get_prg_id()) {
        return false;
    }
    if (first_hit_from_left->get_prg_path() < first_hit_from_right->get_prg_path()) {
        return true;
    }
    if (first_hit_from_right->get_prg_path() < first_hit_from_left->get_prg_path()) {
        return false;
    }
    if (first_hit_from_left->same_strands() < first_hit_from_right->same_strands()) {
        return true;
    }
    if (first_hit_from_right->same_strands() < first_hit_from_left->same_strands()) {
        return false;
    }
    return false;
}
