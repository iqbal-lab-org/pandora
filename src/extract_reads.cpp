#include <cassert>
#include <vector>
#include <set>
#include "extract_reads.h"
#include "interval.h"
#include "path.h"
#include "minihit.h"
#include "pangenome/ns.cpp"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

vector<Interval> identify_regions(const vector<uint32_t>& covgs, const uint32_t& threshold, const uint32_t& min_length) {
    // threshold is to be less than or equal to in intervals [ , )
    //cout << "identify regions " << threshold << " " << min_length << endl;
    uint32_t start = 0, end = 0;
    vector<Interval> regions;
    bool found_start = false;
    for (uint32_t i=0; i<covgs.size(); ++i){
        if (covgs[i] <= threshold and start == 0){
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

vector<LocalNodePtr> find_interval_in_localpath(const Interval& interval, const vector<LocalNodePtr>& lmp) {
    uint32_t added = 0;
    uint8_t level=0, start_level=0, lowest_level = 0;
    uint16_t start=0, end=0;
    bool found_end = false;
    for (uint_least16_t i=0; i<lmp.size(); ++i){
        if (added <= interval.start and added + lmp[i]->pos.length >= interval.start) {
            start = i;
            start_level = level;
            lowest_level = level;
        }
        if (start_level > 0 and level < start_level) {
            lowest_level = level;
        }
        if (added <= interval.get_end() and added + lmp[i]->pos.length >= interval.get_end()) {
            found_end = true;
        }
        end = i;
        if (found_end and level == lowest_level){
            break;
        }
        added += lmp[i]->pos.length;

        if (lmp[i]->outNodes.size() > 1) {
            level += 1;
        } else {
            level -=1;
        }
    }
    //cout << "found start node " << +start << " level " << +start_level << ", and end node " << +end
    //     << ", with lowest level " << +lowest_level << endl;

    if (start_level > lowest_level) {
        // Now extend the start of the interval found so starts at lowest_level
        for (uint_least16_t i=start; i>0; --i) {
            if (lmp[i-1]->outNodes.size() > 1) {
                start_level -= 1;
            } else {
                start_level += 1;
            }

            if (start_level == lowest_level) {
                start = i-1;
                //cout << "new start " << +start << "at level " << +start_level << endl;
                break;
            }
        }
    }

    vector<LocalNodePtr> sub_localpath(lmp.begin()+start, lmp.begin()+end+1);
    return sub_localpath;
}

set<MinimizerHitPtr, pComp_path> hits_along_path(const set<MinimizerHitPtr, pComp_path>& read_hits,
                                                  const vector<LocalNodePtr>& lmp){
    set<MinimizerHitPtr, pComp_path> subset;

    if (lmp.empty())
        return subset;

    deque<Interval> d;
    for (auto n : lmp) {
        d.push_back(n->pos);
    }

    Path local_path;
    local_path.initialize(d);

    for (auto hit_ptr : read_hits) {
        for (auto interval : local_path.path) {
            if (interval.start > hit_ptr->prg_path.get_end())
                break;
            else if (interval.get_end() < hit_ptr->prg_path.get_start())
                continue;
            else if (hit_ptr->prg_path.is_subpath(local_path)) {
                //and not n.second->path.is_branching(local_path))
                subset.insert(hit_ptr);
                break;
            }
        }
    }

    return subset;
}

void get_read_overlap_coordinates(PanNodePtr pnode, vector<vector<uint32_t>>& read_overlap_coordinates,
                                  vector<LocalNodePtr>& lmp)
{
    read_overlap_coordinates.reserve(pnode->reads.size());
    vector<uint32_t> coordinate;

    auto read_count = 0;
    for (const auto read_ptr : pnode->reads)
    {
        read_count++;
        auto read_hits_along_path = hits_along_path(read_ptr->hits.at(pnode->prg_id), lmp);
        if (read_hits_along_path.size() < 2)
            continue;

        auto hit_ptr_iter = read_hits_along_path.begin();
        uint32_t start = (*hit_ptr_iter)->read_start_position;
        uint32_t end = 0;
        for (const auto hit_ptr : read_hits_along_path)
        {
            start = min(start, hit_ptr->read_start_position);
            end = max(end, hit_ptr->read_start_position + hit_ptr->prg_path.length());
        }

        assert(end > start or assert_msg("Error finding the read overlap coordinates for node " << pnode->name << " and read "
                                          << read_ptr->id << " (the " << read_count << "th on this node)" << endl
                                          << "Found end " << end << " after found start " << start));
        coordinate = {read_ptr->id, start, end, (*hit_ptr_iter)->strand};
        read_overlap_coordinates.push_back(coordinate);
    }

    if (read_overlap_coordinates.size() > 0) {
        sort(read_overlap_coordinates.begin(), read_overlap_coordinates.end(),
             [](const vector<uint32_t>& a, const vector<uint32_t>& b) {
                 for (uint32_t i=0; i<a.size(); ++i) {
                     if (a[i] != b[i]) { return a[i] < b[i]; }}
                 return false;});
    }

}

