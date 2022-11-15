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
    if (lhs->get_kmer_node_id() < rhs->get_kmer_node_id()) {
        return true;
    }
    if (rhs->get_kmer_node_id() < lhs->get_kmer_node_id()) {
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
    if (first_hit_from_left->get_kmer_node_id() < first_hit_from_right->get_kmer_node_id()) {
        return true;
    }
    if (first_hit_from_right->get_kmer_node_id() < first_hit_from_left->get_kmer_node_id()) {
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

std::pair<uint32_t, uint32_t> MinimizerHits::get_strand_counts() const {
    std::vector<char> strands;
    strands.reserve(hits.size());
    for (const MinimizerHitPtr &hit : hits) {
        strands.push_back("-+"[hit->same_strands()]);
    }
    uint32_t plus_strand_count = std::count(strands.begin(), strands.end(), '+');
    uint32_t minus_strand_count = strands.size() - plus_strand_count;
    return std::make_pair(plus_strand_count, minus_strand_count);
}
