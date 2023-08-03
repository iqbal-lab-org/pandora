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

    /*
     * We don't take the PRG into account anymore when sorting clusters because we now
     * randomly multimap equally best clusters to random PRGs
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
    */

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

uint32_t MinimizerHits::read_span_size() const
{
    return back()->get_read_start_position() - front()->get_read_start_position();
}

double MinimizerHits::overlap_amount(const MinimizerHits& cluster) const {
    // Calculate the overlap
    const uint32_t start = std::max(this->front()->get_read_start_position(),
                                    cluster.front()->get_read_start_position());
    const uint32_t end = std::min(this->back()->get_read_start_position(),
                                  cluster.back()->get_read_start_position());

    const bool no_overlap = start > end;
    if (no_overlap) {
        return false;
    }

    const uint32_t overlap = end - start;
    const uint32_t shortest_read_span = std::min(this->read_span_size(), cluster.read_span_size());

    return (double)overlap / shortest_read_span;
}

double MinimizerHits::target_coverage() const {
    return (double)get_number_of_unique_mini_in_cluster() / (double)prg_max_path_lengths->at(front()->get_prg_id());
}

bool MinimizerHits::is_preferred_to(const MinimizerHits& cluster, double minimisers_tolerance) const {
    double margin = minimisers_tolerance *
        std::max(this->get_number_of_unique_mini_in_cluster(), cluster.get_number_of_unique_mini_in_cluster());
    double difference = std::abs((double)(this->get_number_of_unique_mini_in_cluster() - cluster.get_number_of_unique_mini_in_cluster()));

    if (difference > margin) {
        return this->get_number_of_unique_mini_in_cluster() > cluster.get_number_of_unique_mini_in_cluster();
    }

    return this->target_coverage() > cluster.target_coverage();
}
