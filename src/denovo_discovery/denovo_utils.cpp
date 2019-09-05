#include "denovo_discovery/denovo_utils.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


PathComponents find_interval_and_flanks_in_localpath(const Interval &interval,
                                                     const std::vector<LocalNodePtr> &local_node_max_likelihood_path) {
    if (interval.empty()) {
        return PathComponents();
    }

    uint32_t start { 0 };
    uint32_t end { 0 };
    bool found_start { false };
    bool found_end { false };
    uint32_t total_bases_traversed { 0 };
    std::vector<Interval> intervals_found, flank_left_intervals, flank_right_intervals;
    const auto interval_end { interval.get_end() };

    for (const auto &current_node : local_node_max_likelihood_path) {
        total_bases_traversed += current_node->pos.length;
        start = current_node->pos.start;
        const auto current_node_end { current_node->pos.get_end() };

        const auto havent_reached_interval_start { interval.start >= total_bases_traversed };
        if (havent_reached_interval_start) {
            flank_left_intervals.emplace_back(current_node->pos);
            continue;
        }

        const auto found_node_interval_starts_in { (not found_start and interval.start < total_bases_traversed) };
        if (found_node_interval_starts_in) {
            start = current_node_end - (total_bases_traversed - interval.start);
            found_start = true;

            const auto interval_starts_after_start_of_current_node {
                    interval.start > (total_bases_traversed - current_node->pos.length) };
            if (interval_starts_after_start_of_current_node) {
                flank_left_intervals.emplace_back(Interval(current_node->pos.start, start));
            }

            const auto interval_doesnt_end_in_this_node { interval_end > total_bases_traversed };
            if (interval_doesnt_end_in_this_node) {
                intervals_found.emplace_back(Interval(start, current_node_end));
                continue;
            }
        }

        const auto found_node_interval_ends_in { (not found_end and interval_end <= total_bases_traversed) };
        if (found_node_interval_ends_in) {
            end = current_node_end - (total_bases_traversed - interval_end);

            const auto interval_ends_before_end_of_current_node { interval_end < total_bases_traversed };
            if (interval_ends_before_end_of_current_node) {
                flank_right_intervals.emplace_back(Interval(end, current_node_end));
            }
            intervals_found.emplace_back(Interval(start, end));
            found_end = true;
            continue;
        }

        const auto node_is_completely_contained_in_interval {
                interval.start < total_bases_traversed and interval_end > total_bases_traversed };
        if (node_is_completely_contained_in_interval) {
            intervals_found.emplace_back(Interval(start, start + current_node->pos.length));
            continue;
        }

        const auto past_interval_end { found_end and interval_end < total_bases_traversed };
        if (past_interval_end) {
            flank_right_intervals.emplace_back(current_node->pos);
        }
    }

    prg::Path interval_as_path;
    interval_as_path.initialize(intervals_found);
    prg::Path left_flank_path;
    left_flank_path.initialize(flank_left_intervals);
    prg::Path right_flank_path;
    right_flank_path.initialize(flank_right_intervals);

    PathComponents found_components { left_flank_path, interval_as_path, right_flank_path };

    return found_components;
}


std::vector<MinimizerHitPtr>
find_hits_inside_path(const std::vector<MinimizerHitPtr> &read_hits, const prg::Path &local_path) {
    std::vector<MinimizerHitPtr> hits_inside_local_path;

    if (local_path.empty()) {
        return hits_inside_local_path;
    }

    for (const auto &current_read_hit : read_hits) {
        for (const auto &interval : local_path) {
            const auto &prg_path_of_current_read_hit = current_read_hit->get_prg_path();
            const auto hit_is_to_left_of_path_start { interval.start > prg_path_of_current_read_hit.get_end() };
            const auto hit_is_to_right_of_current_interval {
                    interval.get_end() < prg_path_of_current_read_hit.get_start() };

            if (hit_is_to_left_of_path_start) {
                break;
            } else if (hit_is_to_right_of_current_interval) {
                continue;
            } else if (prg_path_of_current_read_hit.is_subpath(local_path)) {
                hits_inside_local_path.push_back(current_read_hit);
                break;
            }
        }
    }

    return hits_inside_local_path;
}


ReadCoordinate::ReadCoordinate(uint32_t id, uint32_t start, uint32_t end, bool is_forward)
        : id(id), start(start), end(end), is_forward(is_forward) {}


bool ReadCoordinate::operator<(const ReadCoordinate &y) const {
    if (this->id < y.id) {
        return true;
    }
    if (this->id > y.id) {
        return false;
    }

    if (this->start < y.start) {
        return true;
    }
    if (this->start > y.start) {
        return false;
    }

    if (this->end < y.end) {
        return true;
    }
    if (this->end > y.end) {
        return false;
    }

    if (this->is_forward and not y.is_forward) {
        return true;
    }
    if (y.is_forward and not this->is_forward) {
        return false;
    }

    return false;
}


bool ReadCoordinate::operator==(const ReadCoordinate &y) const {
    return ((this->id == y.id) and (this->start == y.start) and (this->end == y.end) and
            (this->is_forward == y.is_forward));
}


bool ReadCoordinate::operator!=(const ReadCoordinate &y) const {
    return (!(*this == y));
}


size_t std::hash<ReadCoordinate>::operator()(const ReadCoordinate &coordinate) const {
    return hash<int>()(coordinate.start) ^ hash<int>()(coordinate.id) ^ hash<int>()(coordinate.end) ^
           hash<bool>()(coordinate.is_forward);
}


std::ostream &operator<<(std::ostream &out, ReadCoordinate const &y) {
    out << "{" << y.id << ", " << y.start << ", " << y.end << ", " << y.is_forward << "}";
    return out;
}


PathComponents::PathComponents() = default;


PathComponents::PathComponents(prg::Path flank_left, prg::Path slice, prg::Path flank_right)
        : flank_left(std::move(flank_left)), slice(std::move(slice)), flank_right(std::move(flank_right)) {}


bool PathComponents::operator==(const PathComponents &other) const {
    return (flank_left == other.flank_left and flank_right == other.flank_right and slice == other.slice);

}


bool PathComponents::operator!=(const PathComponents &other) const {
    return (!(*this == other));
}