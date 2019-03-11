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
typedef prg::Path Path;

std::vector<Interval>
identify_regions(const std::vector<uint32_t> &covgs, const uint32_t &threshold, const uint32_t &min_length) {
    // threshold is to be less than or equal to in intervals [ , )
    const auto padding = 21;
    uint32_t start = 0, end = 0;
    std::vector<Interval> regions;
    bool found_start = false;

    for (uint32_t i = 0; i < covgs.size(); ++i) {
        if (covgs[i] <= threshold and start == 0) {
            start = i;
            end = 0;
            found_start = true;
        } else if (found_start and covgs[i] > threshold and end == 0) {
            end = i;
            if (end - start >= min_length) {
                const auto interval_start = (start <= padding) ? 0 : start - padding;
                regions.emplace_back(Interval(interval_start, end + padding));
            }
            start = 0;
        }
    }
    if (found_start and end == 0 and covgs.size() - start >= min_length) {
        end = covgs.size();
        regions.emplace_back(Interval(start, end));
    }
    return regions;
}

prg::Path
find_interval_in_localpath(const Interval &interval,
                           const std::vector<LocalNodePtr> &lmp,
                           const uint32_t &buff) {
    uint32_t start = 0, end = 0;
    bool found_start = false;
    uint32_t total = 0;
    std::vector<Interval> sub_localpath_intervals;

    for (uint_least16_t i = 0; i < lmp.size(); ++i) {
        const auto len = lmp[i]->pos.length;
        total += len;
        start = lmp[i]->pos.start;
        const auto lmp_end = lmp[i]->pos.get_end();
        const auto interval_end = interval.get_end();

        if (interval.start > total) {  // still haven't reached our interval start
            continue;
        }
        if (not found_start and interval.start <= total) {  // we have found the local node in which interval starts
            start = lmp_end - (total - interval.start);
            found_start = true;
            if (interval_end > total) {
                sub_localpath_intervals.emplace_back(Interval(start, lmp_end - start));
                continue;
            }
        }
        if (interval_end <= total) {  // we have found the local node in which our interval ends
            end = lmp_end - (total - interval_end);
            sub_localpath_intervals.emplace_back(Interval(start, end));
            break;
        }

        if (interval.start < total and interval_end > total) {  // we have found a local node which falls in the middle of our interval range
            sub_localpath_intervals.emplace_back(Interval(start, len));
        }
    }

    prg::Path sub_localpath;
    sub_localpath.initialize(sub_localpath_intervals);

    BOOST_LOG_TRIVIAL(debug) << "For interval " << interval.start << ", " << interval.get_end() << " got sub_localpath: ";
    for (auto &x : sub_localpath.path) {
        BOOST_LOG_TRIVIAL(debug) << x.start << ", " << x.get_end();
    }

    return sub_localpath;
}

std::set<MinimizerHitPtr, pComp_path>
hits_inside_path(const std::set<MinimizerHitPtr, pComp_path> &read_hits, const prg::Path &local_path) {
    std::set<MinimizerHitPtr, pComp_path> subset;

    if (local_path.path.empty()) {
        return subset;
    }

    for (const auto &hit_ptr : read_hits) {
        for (const auto &interval : local_path.path) {
            if (interval.start > hit_ptr->prg_path.get_end()) {
                break;
            } else if (interval.get_end() < hit_ptr->prg_path.get_start()) {
                continue;
            } else if (hit_ptr->prg_path.is_subpath(local_path)) {
                subset.insert(hit_ptr);
                break;
            }
        }
    }

    return subset;
}


std::set<ReadCoordinate>
get_read_overlap_coordinates(const PanNodePtr &pnode, const prg::Path &local_path, const uint32_t &min_number_hits) {
    std::set<ReadCoordinate> read_overlap_coordinates = {};
    BOOST_LOG_TRIVIAL(debug) << "read ids covering pannode to consider ";
    for (const auto &r : pnode->reads)
        BOOST_LOG_TRIVIAL(debug) << r->id;
    for (const auto &read_ptr: pnode->reads) {
        auto read_hits_inside_path = hits_inside_path(read_ptr->hits.at(pnode->prg_id), local_path);
        BOOST_LOG_TRIVIAL(debug) << "read_hits_inside_path size " << read_hits_inside_path.size();
        if (read_hits_inside_path.size() < min_number_hits)
            continue;

        auto hit_ptr_iter = read_hits_inside_path.begin();
        uint32_t start = (*hit_ptr_iter)->read_start_position;
        uint32_t end = 0;
        for (const auto &hit_ptr : read_hits_inside_path) {
            start = std::min(start, hit_ptr->read_start_position);
            end = std::max(end, hit_ptr->read_start_position + hit_ptr->prg_path.length());
        }
        assert(end > start);

        ReadCoordinate coordinate = {read_ptr->id, start, end, (*hit_ptr_iter)->strand};
        read_overlap_coordinates.insert(coordinate);
    }
    return read_overlap_coordinates;
}

void denovo_discovery::add_pnode_coordinate_pairs(std::vector<std::shared_ptr<LocalPRG>> &prgs,
        std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
        const PanNodePtr &pnode,
        const std::vector<LocalNodePtr> &local_node_path,
        const std::vector<KmerNodePtr> &kmer_node_path,
        const uint32_t &padding_size,
        const uint32_t &low_coverage_threshold,
        const uint32_t &interval_min_length,
        const uint32_t &min_number_hits) {
    uint32_t sample_id = 0;
    auto covgs = get_covgs_along_localnode_path(pnode, local_node_path, kmer_node_path, sample_id);
    auto intervals = identify_regions(covgs, low_coverage_threshold, interval_min_length);
    if (intervals.empty())
        return;

    BOOST_LOG_TRIVIAL(debug) << "there are " << intervals.size() << " intervals";
    for (const auto &interval: intervals) {
        BOOST_LOG_TRIVIAL(debug) << "Looking at interval: " << interval;
        BOOST_LOG_TRIVIAL(debug) << "For gene: " << pnode->get_name();

        auto local_path = find_interval_in_localpath(interval, local_node_path, padding_size);

        auto read_overlap_coordinates = get_read_overlap_coordinates(pnode, local_path,
                                                                     min_number_hits);
        BOOST_LOG_TRIVIAL(debug) << "there are " << read_overlap_coordinates.size() << " read_overlap_coordinates";

        GeneIntervalInfo interval_info{
                pnode,
                interval,
                prgs[pnode->prg_id]->string_along_path(local_path)
        };
        for (const auto &read_coordinate: read_overlap_coordinates)
            pangraph_coordinate_pairs.insert(std::make_pair(read_coordinate, interval_info));
    }
}

bool ReadCoordinate::operator<(const ReadCoordinate &y) const {
    if (this->id < y.id)
        return true;
    if (this->id > y.id)
        return false;

    if (this->start < y.start)
        return true;
    if (this->start > y.start)
        return false;

    if (this->end < y.end)
        return true;
    if (this->end > y.end)
        return false;

    if (this->strand and not y.strand)
        return true;
    if (y.strand and not this->strand)
        return false;

    return false;
}

bool ReadCoordinate::operator==(const ReadCoordinate &y) const {
    return (
            (this->id == y.id)
            and (this->start == y.start)
            and (this->end == y.end)
            and (this->strand == y.strand)
    );
}

bool ReadCoordinate::operator!=(const ReadCoordinate &y) const {
    return (!(*this == y));
}

std::ostream &operator<<(std::ostream &out, ReadCoordinate const &y) {
    out << "{" << y.id << ", " << y.start << ", " << y.end << ", " << y.strand << "}";
    return out;
}

std::map<GeneIntervalInfo, ReadPileup>
denovo_discovery::collect_read_pileups(
        const std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
        const boost::filesystem::path &readfilepath,
        const uint32_t &padding_size) {
    FastaqHandler readfile(readfilepath.string());

    std::map<GeneIntervalInfo, ReadPileup> pileup = {};

    uint32_t last_id = 0;
    for (const auto &pangraph_coordinate: pangraph_coordinate_pairs) {
        const auto &read_coordinate = pangraph_coordinate.first;
        const auto &interval_info = pangraph_coordinate.second;

        assert (last_id <= read_coordinate.id);
        readfile.get_id(read_coordinate.id);
        uint32_t start = padding_size > read_coordinate.start ? 0 : read_coordinate.start - padding_size;
        uint32_t end = std::min(read_coordinate.end + padding_size,
                                (uint32_t) readfile.read.length());

        auto sequence = readfile.read.substr(start, end - start);
        if (!read_coordinate.strand)
            sequence = rev_complement(sequence);
        pileup[interval_info].emplace_back(sequence);
        last_id = read_coordinate.id;
    }
    return pileup;
}