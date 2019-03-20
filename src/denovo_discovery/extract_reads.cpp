#include <cassert>
#include <vector>
#include <set>
#include <utility>
#include "denovo_discovery/extract_reads.h"
#include "interval.h"
#include "minihit.h"
#include "pangenome/ns.cpp"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "localPRG.h"
#include "fastaq_handler.h"
#include "fastaq.h"
#include "utils.h"
#include "prg/path.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


using std::make_pair;

std::vector<Interval>
identify_regions(const std::vector<uint32_t> &covg_at_each_position, const uint32_t &min_required_covg,
                 const uint32_t &min_length, const uint_least16_t &padding) {
    uint32_t start_position { 0 };
    uint32_t end_position { 0 };
    std::vector<Interval> identified_regions;
    bool found_start { false };

    for (uint32_t current_pos = 0; current_pos < covg_at_each_position.size(); ++current_pos) {
        if (covg_at_each_position[current_pos] <= min_required_covg and start_position == 0) {
            start_position = current_pos;
            end_position = 0;
            found_start = true;
        } else if (found_start and covg_at_each_position[current_pos] > min_required_covg and end_position == 0) {
            end_position = current_pos;
            if (end_position - start_position >= min_length) {
                const auto interval_start { (start_position <= padding) ? 0 : start_position - padding };
                identified_regions.emplace_back(Interval(interval_start, end_position + padding));
            }
            start_position = 0;
        }
    }
    if (found_start and end_position == 0 and covg_at_each_position.size() - start_position >= min_length) {
        end_position = covg_at_each_position.size();
        identified_regions.emplace_back(Interval(start_position, end_position));
    }
    return identified_regions;
}


PathComponents
find_interval_in_localpath(const Interval &interval, const std::vector<LocalNodePtr> &local_node_max_likelihood_path) {
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

        const auto havent_reached_interval_start { interval.start > total_bases_traversed };
        if (havent_reached_interval_start) {
            flank_left_intervals.emplace_back(current_node->pos);
            continue;
        }

        const auto found_node_interval_starts_in { (not found_start and interval.start <= total_bases_traversed) };
        if (found_node_interval_starts_in) {
            start = current_node_end - (total_bases_traversed - interval.start);
            found_start = true;
            flank_left_intervals.emplace_back(Interval(current_node->pos.start, start));

            const auto interval_doesnt_end_in_this_node { interval_end > total_bases_traversed };
            if (interval_doesnt_end_in_this_node) {
                intervals_found.emplace_back(Interval(start, current_node_end));
                continue;
            }
        }

        const auto found_node_interval_ends_in { (not found_end and interval_end <= total_bases_traversed) };
        if (found_node_interval_ends_in) {
            end = current_node_end - (total_bases_traversed - interval_end);

            if (interval_end != total_bases_traversed) {
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

std::set<MinimizerHitPtr, pComp_path>
find_hits_inside_path(const std::set<MinimizerHitPtr, pComp_path> &read_hits, const prg::Path &local_path) {
    std::set<MinimizerHitPtr, pComp_path> hits_inside_local_path;

    if (local_path.path.empty()) {
        return hits_inside_local_path;
    }

    for (const auto &current_read_hit : read_hits) {
        for (const auto &interval : local_path.path) {
            const auto hit_is_to_left_of_path_start { interval.start > current_read_hit->prg_path.get_end() };
            const auto hit_is_to_right_of_current_interval {
                    interval.get_end() < current_read_hit->prg_path.get_start() };

            if (hit_is_to_left_of_path_start) {
                break;
            } else if (hit_is_to_right_of_current_interval) {
                continue;
            } else if (current_read_hit->prg_path.is_subpath(local_path)) {
                hits_inside_local_path.insert(current_read_hit);
                break;
            }
        }
    }

    return hits_inside_local_path;
}


std::set<ReadCoordinate>
get_read_overlap_coordinates(const PanNodePtr &pan_node, const prg::Path &local_path, const uint32_t &min_number_hits) {
    std::set<ReadCoordinate> read_overlap_coordinates;

    for (const auto &current_read: pan_node->reads) {
        const auto read_hits_inside_path { find_hits_inside_path(current_read->hits.at(pan_node->prg_id), local_path) };

        if (read_hits_inside_path.size() < min_number_hits) {
            continue;
        }

        // todo create function to do the following
        const auto read_hits_iter { read_hits_inside_path.cbegin() };
        uint32_t start { (*read_hits_iter)->read_start_position };
        uint32_t end { 0 };

        for (const auto &read_hit : read_hits_inside_path) {
            start = std::min(start, read_hit->read_start_position);
            end = std::max(end, read_hit->read_start_position + read_hit->prg_path.length());
        }

        assert(end > start);

        read_overlap_coordinates.emplace(current_read->id, start, end, (*read_hits_iter)->strand);
    }
    return read_overlap_coordinates;
}

void denovo_discovery::add_pnode_coordinate_pairs(const std::vector<std::shared_ptr<LocalPRG>> &all_prgs,
                                                  std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
                                                  const PanNodePtr &pan_node,
                                                  const std::vector<LocalNodePtr> &local_node_path,
                                                  const std::vector<KmerNodePtr> &kmer_node_path,
                                                  const uint_least16_t &interval_padding,
                                                  const uint32_t &interval_min_required_covg,
                                                  const uint32_t &interval_min_length) {
    const uint32_t sample_id { 0 };
    const auto covgs_along_localnode_path {
            get_covgs_along_localnode_path(pan_node, local_node_path, kmer_node_path, sample_id) };
    const auto candidate_intervals {
            identify_regions(covgs_along_localnode_path, interval_min_required_covg, interval_min_length,
                             interval_padding) };

    if (candidate_intervals.empty()) {
        return;
    }

    BOOST_LOG_TRIVIAL(debug) << "there are " << candidate_intervals.size() << " intervals";

    for (const auto &current_interval: candidate_intervals) {
        BOOST_LOG_TRIVIAL(debug) << "Looking at interval: " << current_interval;
        BOOST_LOG_TRIVIAL(debug) << "For gene: " << pan_node->get_name();

        const auto interval_path_components { find_interval_in_localpath(current_interval, local_node_path) };
        const auto read_overlap_coordinates { get_read_overlap_coordinates(pan_node, interval_path_components.slice) };

        BOOST_LOG_TRIVIAL(debug) << "there are " << read_overlap_coordinates.size() << " read_overlap_coordinates";

        const auto sequence_of_interval {
                all_prgs[pan_node->prg_id]->string_along_path(interval_path_components.slice) };
        const auto sequence_of_left_flank {
                all_prgs[pan_node->prg_id]->string_along_path(interval_path_components.flank_left) };
        const auto sequence_of_right_flank {
                all_prgs[pan_node->prg_id]->string_along_path(interval_path_components.flank_right) };
        const GeneIntervalInfo interval_info { pan_node, current_interval, sequence_of_interval, sequence_of_left_flank,
                                               sequence_of_right_flank, };

        for (const auto &read_coordinate: read_overlap_coordinates) {
            pangraph_coordinate_pairs.emplace(std::make_pair(read_coordinate, interval_info));
        }
    }
}


ReadCoordinate::ReadCoordinate(uint32_t id, uint32_t start, uint32_t end, bool strand) : id(id), start(start), end(end),
                                                                                         strand(strand) {}

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

    if (this->strand and not y.strand) {
        return true;
    }
    if (y.strand and not this->strand) {
        return false;
    }

    return false;
}

bool ReadCoordinate::operator==(const ReadCoordinate &y) const {
    return ((this->id == y.id) and (this->start == y.start) and (this->end == y.end) and (this->strand == y.strand));
}

bool ReadCoordinate::operator!=(const ReadCoordinate &y) const {
    return (!(*this == y));
}

std::ostream &operator<<(std::ostream &out, ReadCoordinate const &y) {
    out << "{" << y.id << ", " << y.start << ", " << y.end << ", " << y.strand << "}";
    return out;
}

std::map<GeneIntervalInfo, ReadPileup> denovo_discovery::collect_read_pileups(
        const std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
        const boost::filesystem::path &readfilepath) {
    FastaqHandler readfile(readfilepath.string());
    std::map<GeneIntervalInfo, ReadPileup> pileup;
    uint32_t last_id { 0 };

    for (const auto &pangraph_coordinate: pangraph_coordinate_pairs) {
        const auto &read_coordinate { pangraph_coordinate.first };
        const auto &interval_info { pangraph_coordinate.second };

        assert (last_id <= read_coordinate.id);
        readfile.get_id(read_coordinate.id);
        const auto end { std::min(read_coordinate.end, (uint32_t) readfile.read.length()) };
        auto sequence = readfile.read.substr(read_coordinate.start, end - read_coordinate.start);

        if (!read_coordinate.strand) {
            sequence = rev_complement(sequence);
        }

        pileup[interval_info].emplace_back(sequence);
        last_id = read_coordinate.id;
    }
    return pileup;
}

PathComponents::PathComponents() {}

PathComponents::PathComponents(prg::Path flank_left, prg::Path slice, prg::Path flank_right) : flank_left(flank_left),
                                                                                               slice(slice),
                                                                                               flank_right(
                                                                                                       flank_right) {}

bool PathComponents::operator==(const PathComponents &other) const {
    return (flank_left == other.flank_left and flank_right == other.flank_right and slice == other.slice);

}

bool PathComponents::operator!=(const PathComponents &other) const {
    return (!(*this == other));
}