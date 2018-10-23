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

using namespace std;
using std::make_pair;
typedef prg::Path Path;

vector<Interval>
identify_regions(const vector<uint32_t> &covgs, const uint32_t &threshold, const uint32_t &min_length) {
    // threshold is to be less than or equal to in intervals [ , )
    //cout << "identify regions " << threshold << " " << min_length << endl;
    uint32_t start = 0, end = 0;
    vector<Interval> regions;
    bool found_start = false;
    for (uint32_t i = 0; i < covgs.size(); ++i) {
        if (covgs[i] <= threshold and start == 0) {
            //cout << "covgs[" << i << "] " << covgs[i] << "<=" << threshold << " threshold" << endl;
            start = i;
            end = 0;
            found_start = true;
        } else if (found_start and covgs[i] > threshold and end == 0) {
            end = i;
            if (end - start >= min_length) {
                regions.emplace_back(Interval(start, end));
                //cout << "end - start " << end << "-" << start << " = " << end - start << " >= " << min_length
                // << " min_length" << endl;
            }
            start = 0;
        }
    }
    if (found_start and end == 0 and covgs.size() - start >= min_length) {
        end = covgs.size();
        regions.emplace_back(Interval(start, end));
        //cout << "end - start " << end << "-" << start << " = " << end - start << " >= " << min_length
        // << " min_length" << endl;
    }
    return regions;
}

vector<LocalNodePtr>
find_interval_in_localpath(const Interval &interval, const vector<LocalNodePtr> &lmp,
                           const unsigned int &buff) {
    uint32_t added = 0;
    uint16_t start = 0, end = 0;
    bool found_end = false;
    for (uint_least16_t i = 0; i < lmp.size(); ++i) {
        if (added + buff <= interval.start and added + buff+ lmp[i]->pos.length >= interval.start) {
            start = i;
        }
        if (added <= interval.get_end() + buff and added + lmp[i]->pos.length >= interval.get_end() + buff) {
            found_end = true;
        }
        end = i;
        if (found_end) {
            break;
        }
        added += lmp[i]->pos.length;
    }

    vector<LocalNodePtr> sub_localpath(lmp.begin() + start, lmp.begin() + end + 1);
    BOOST_LOG_TRIVIAL(debug) << "Interval" << interval << " along path ";
    for (const auto &l : lmp)
        BOOST_LOG_TRIVIAL(debug) << *l << " ";
    BOOST_LOG_TRIVIAL(debug) << endl;
    BOOST_LOG_TRIVIAL(debug) << "(plus buffer " << +buff << ") corresponds to sub path ";
    for (const auto &l : sub_localpath)
        BOOST_LOG_TRIVIAL(debug) << *l << " ";
    BOOST_LOG_TRIVIAL(debug) << endl;
    return sub_localpath;
}

set<MinimizerHitPtr, pComp_path> hits_inside_path(const set<MinimizerHitPtr, pComp_path> &read_hits,
                                                 const vector<LocalNodePtr> &lmp) {
    set<MinimizerHitPtr, pComp_path> subset;

    if (lmp.empty()) {
        return subset;
    }

    deque<Interval> d;
    for (const auto &n : lmp) {
        d.push_back(n->pos);
    }

    prg::Path local_path;
    local_path.initialize(d);

    for (const auto &hit_ptr : read_hits) {
        for (const auto &interval : local_path.path) {
            if (interval.start > hit_ptr->prg_path.get_end()) {
                break;
            } else if (interval.get_end() < hit_ptr->prg_path.get_start()) {
                continue;
            } else if (hit_ptr->prg_path.is_subpath(local_path)) {
                //and not n.second->path.is_branching(local_path))
                subset.insert(hit_ptr);
                break;
            }
        }
    }

    return subset;
}


std::set<ReadCoordinate> get_read_overlap_coordinates(const PanNodePtr &pnode,
                                                      const std::vector<LocalNodePtr> &local_node_path,
                                                      const uint32_t& min_number_hits) {
    std::set<ReadCoordinate> read_overlap_coordinates = {};
    cout << "read ids " << endl;
    for (const auto &r : pnode->reads)
        cout << r->id << endl;
    for (const auto &read_ptr: pnode->reads) {
        auto read_hits_inside_path = hits_inside_path(read_ptr->hits.at(pnode->prg_id), local_node_path);
        cout << "read_hits_inside_path size " << read_hits_inside_path.size() << endl;
        if (read_hits_inside_path.size() < min_number_hits)
            continue;

        auto hit_ptr_iter = read_hits_inside_path.begin();
        uint32_t start = (*hit_ptr_iter)->read_start_position;
        uint32_t end = 0;
        for (const auto &hit_ptr : read_hits_inside_path) {
            start = min(start, hit_ptr->read_start_position);
            end = max(end, hit_ptr->read_start_position + hit_ptr->prg_path.length());
        }
        assert(end > start);

        ReadCoordinate coordinate = {read_ptr->id, start, end, (*hit_ptr_iter)->strand};
        read_overlap_coordinates.insert(coordinate);
    }
    return read_overlap_coordinates;
}

void denovo_discovery::add_pnode_coordinate_pairs(
        std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
        const PanNodePtr &pnode,
        const std::vector<LocalNodePtr> &local_node_path,
        const std::vector<KmerNodePtr> &kmer_node_path,
        const uint32_t &padding_size,
        const uint32_t &low_coverage_threshold,
        const uint32_t &interval_min_length,
        const uint32_t& min_number_hits) {
    auto covgs = get_covgs_along_localnode_path(pnode, local_node_path, kmer_node_path);
    auto intervals = identify_regions(covgs, low_coverage_threshold, interval_min_length);
    if (intervals.empty())
        return;

    cout << "there are " << intervals.size() << " intervals" << endl;
    for (const auto &interval: intervals) {
        BOOST_LOG_TRIVIAL(debug) << "Looking at interval: " << interval;
        BOOST_LOG_TRIVIAL(debug) << "For gene: " << pnode->get_name();

        auto sub_local_node_path = find_interval_in_localpath(interval, local_node_path, padding_size);
        cout << "sub local node path ";
        for (const auto &l : sub_local_node_path)
            cout << *l;
        cout << endl;
        auto read_overlap_coordinates = get_read_overlap_coordinates(pnode, sub_local_node_path, min_number_hits);
        cout << "there are " << read_overlap_coordinates.size() << " read_overlap_coordinates" << endl;

        GeneIntervalInfo interval_info{
                pnode,
                interval,
                LocalPRG::string_along_path(sub_local_node_path)
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

std::ostream &operator<<(std::ostream &out, ReadCoordinate const &y) {
    out << "{" << y.id << ", " << y.start << ", " << y.end << ", " << y.strand << "}";
    return out;
}

std::map<GeneIntervalInfo, ReadPileup>
denovo_discovery::collect_read_pileups(const std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> &pangraph_coordinate_pairs,
                     const boost::filesystem::path &readfilepath,
                     const uint32_t &padding_size) {
    FastaqHandler readfile(readfilepath.string());

    std::map<GeneIntervalInfo, ReadPileup> pileup = {};

    for (const auto &pangraph_coordinate: pangraph_coordinate_pairs) {
        const auto &read_coordinate = pangraph_coordinate.first;
        const auto &interval_info = pangraph_coordinate.second;

        readfile.get_id(read_coordinate.id);
        uint32_t start = padding_size > read_coordinate.start ? 0 : read_coordinate.start - padding_size;
        uint32_t end = std::min(read_coordinate.end + padding_size,
                                (uint32_t) readfile.read.length());

        auto sequence = readfile.read.substr(start, end - start);
        pileup[interval_info].emplace_back(sequence);
    }
    return pileup;
}