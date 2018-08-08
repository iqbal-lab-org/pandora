#include <cassert>
#include <vector>
#include <set>
#include "extract_reads.h"
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
    uint8_t level = 0, start_level = 0, lowest_level = 0;
    uint16_t start = 0, end = 0;
    bool found_end = false;
    for (uint_least16_t i = 0; i < lmp.size(); ++i) {
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
        if (found_end and level == lowest_level) {
            break;
        }
        added += lmp[i]->pos.length;

        if (lmp[i]->outNodes.size() > 1) {
            level += 1;
        } else {
            level -= 1;
        }
    }
    //cout << "found start node " << +start << " level " << +start_level << ", and end node " << +end
    //     << ", with lowest level " << +lowest_level << endl;
    if (buff > 0) {
        auto start_buffer_added{lmp[start]->pos.length};
        auto end_buffer_added{lmp[end]->pos.length};
        while (start > 0 and start_buffer_added < buff) {
            start -= 1;
            start_buffer_added += lmp[start]->pos.length;
        }
        while (end < lmp.size() and end_buffer_added < buff) {
            end += 1;
            end_buffer_added += lmp[end]->pos.length;
        }
    }

    if (start_level > lowest_level) {
        // Now extend the start of the interval found so starts at lowest_level
        for (uint_least16_t i = start; i > 0; --i) {
            if (lmp[i - 1]->outNodes.size() > 1) {
                start_level -= 1;
            } else {
                start_level += 1;
            }

            if (start_level == lowest_level) {
                start = i - 1;
                break;
            }
        }
    }

    vector<LocalNodePtr> sub_localpath(lmp.begin() + start, lmp.begin() + end + 1);
    return sub_localpath;
}

set<MinimizerHitPtr, pComp_path> hits_along_path(const set<MinimizerHitPtr, pComp_path> &read_hits,
                                                 const vector<LocalNodePtr> &lmp) {
    set<MinimizerHitPtr, pComp_path> subset;

    if (lmp.empty()) {
        return subset;
    }

    deque<Interval> d;
    for (auto n : lmp) {
        d.push_back(n->pos);
    }

    prg::Path local_path;
    local_path.initialize(d);

    for (auto hit_ptr : read_hits) {
        for (auto interval : local_path.path) {
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

void get_read_overlap_coordinates(PanNodePtr pnode, vector<vector<uint32_t>> &read_overlap_coordinates,
                                  vector<LocalNodePtr> &lmp) {
    read_overlap_coordinates.clear();
    read_overlap_coordinates.reserve(pnode->reads.size());
    vector<uint32_t> coordinate;

    auto read_count = 0;
    for (const auto read_ptr : pnode->reads) {
        read_count++;
        auto read_hits_along_path = hits_along_path(read_ptr->hits.at(pnode->prg_id), lmp);
        if (read_hits_along_path.size() < 2) {
            continue;
        }

        auto hit_ptr_iter = read_hits_along_path.begin();
        uint32_t start = (*hit_ptr_iter)->read_start_position;
        uint32_t end = 0;
        for (const auto hit_ptr : read_hits_along_path) {
            start = min(start, hit_ptr->read_start_position);
            end = max(end, hit_ptr->read_start_position + hit_ptr->prg_path.length());
        }

        assert(end > start or
               assert_msg("Error finding the read overlap coordinates for node " << pnode->name << " and read "
                                                                                 << read_ptr->id << " (the "
                                                                                 << read_count << "th on this node)"
                                                                                 << endl
                                                                                 << "Found end " << end
                                                                                 << " after found start " << start));
        coordinate = {read_ptr->id, start, end, (*hit_ptr_iter)->strand};
        read_overlap_coordinates.push_back(coordinate);
    }

    if (read_overlap_coordinates.size() > 0) {
        sort(read_overlap_coordinates.begin(), read_overlap_coordinates.end(),
             [](const vector<uint32_t> &a, const vector<uint32_t> &b) {
                 for (uint32_t i = 0; i < a.size(); ++i) {
                     if (a[i] != b[i]) { return a[i] < b[i]; }
                 }
                 return false;
             });
    }

}

void save_read_strings_to_denovo_assemble(const string &readfilepath,
                                          const string &outdir,
                                          const PanNodePtr pnode,
                                          const vector<LocalNodePtr> &lmp,
                                          const vector<KmerNodePtr> &kmp,
                                          const unsigned int &buff,
                                          const uint32_t &threshold,
                                          const uint32_t &min_length
) {

    // level for boost logging
    logging::core::get()->set_filter(logging::trivial::severity >= g_log_level);

    vector<uint32_t> covgs = get_covgs_along_localnode_path(pnode, lmp, kmp);
    vector<Interval> intervals = identify_regions(covgs, threshold, min_length);

    if (intervals.empty()) {
        return;
    }

    cout << now() << "Save mapped read strings and coordinates" << endl;
    make_dir(outdir);

    FastaqHandler readfile(readfilepath);
    Fastaq fa;
    uint32_t start, end;
    vector<vector<uint32_t>> read_overlap_coordinates;
    vector<LocalNodePtr> sub_lmp;

    for (auto interval : intervals) {

        BOOST_LOG_TRIVIAL(debug) << "Looking at interval: " << interval;

        sub_lmp = find_interval_in_localpath(interval, lmp, buff);
        get_read_overlap_coordinates(pnode, read_overlap_coordinates, sub_lmp);

        uint16_t j = 0;
        for (auto coord : read_overlap_coordinates) {
            cout << "\nLooking at coordinate j = " << +j << " {" << coord[0] << "," << coord[1] << "," << coord[2]
                 << ","
                 << coord[3] << "}" << endl;
            j++;
            readfile.get_id(coord[0]);
            start = (uint32_t) std::max((int32_t) coord[1] - (int32_t) buff, 0);
            end = min(coord[2] + (uint32_t) buff, (uint32_t) readfile.read.length());

            assert(coord[1] < coord[2] or assert_msg("For read #" << coord[0] << " " << readfile.name << " pannode "
                                                                  << pnode->get_name() << " interval ["
                                                                  << coord[1] << ", " << coord[2]
                                                                  << ") found start " << start
                                                                  << " end " << end
                                                                  << " compared to readlength "
                                                                  << readfile.read.length()));
            assert(start <= coord[1] or assert_msg("For read #" << coord[0] << " " << readfile.name << " pannode "
                                                                << pnode->get_name() << " interval ["
                                                                << coord[1] << ", " << coord[2]
                                                                << ") found start " << start
                                                                << " end " << end
                                                                << " compared to readlength "
                                                                << readfile.read.length()));
            assert(start <= readfile.read.length() or
                   assert_msg("For read #" << coord[0] << " " << readfile.name << " pannode "
                                           << pnode->get_name() << " interval ["
                                           << coord[1] << ", " << coord[2]
                                           << ") found start " << start
                                           << " end " << end
                                           << " compared to readlength "
                                           << readfile.read.length()));
            assert(coord[2] <= readfile.read.length() or
                   assert_msg("For read #" << coord[0] << " " << readfile.name << " pannode "
                                           << pnode->get_name() << " interval ["
                                           << coord[1] << ", " << coord[2]
                                           << ") found start " << start
                                           << " end " << end
                                           << " compared to readlength "
                                           << readfile.read.length()));
            assert(end >= coord[2] or assert_msg("For read #" << coord[0] << " " << readfile.name << " pannode "
                                                              << pnode->get_name() << " interval ["
                                                              << coord[1] << ", " << coord[2]
                                                              << ") found start " << start
                                                              << " end " << end
                                                              << " compared to readlength "
                                                              << readfile.read.length()));
            assert(start < end or assert_msg("For read #" << coord[0] << " " << readfile.name << " pannode "
                                                          << pnode->get_name() << " interval ["
                                                          << coord[1] << ", " << coord[2]
                                                          << ") found start " << start
                                                          << " end " << end
                                                          << " compared to readlength "
                                                          << readfile.read.length()));

            string header = "pandora: " + to_string(coord[0]) + " " + to_string(start) + ":" + to_string(end);
            if (coord[3] == true) {
                header += " + ";
            } else {
                header += " - ";
            }
            string sequence = readfile.read.substr(start, end - start);

            fa.add_entry(readfile.name, sequence, header);
        }
        const auto filepath = outdir + "/" + pnode->get_name() + "." + to_string(interval.start) + "-" +
                              to_string(interval.get_end()) + ".fa";
        fa.save(filepath);
        BOOST_LOG_TRIVIAL(debug) << "Graph slice for local assembly saved as " << filepath;

        // get sub_lmp path as string
        const auto sub_lmp_as_string{LocalPRG::string_along_path(sub_lmp)};
        BOOST_LOG_TRIVIAL(debug) << "sub_lmp for interval is " << sub_lmp_as_string;

        const unsigned long max_path_length{
                sub_lmp_as_string.length() + (interval.length * 5)};  // arbitrary at the moment


        BOOST_LOG_TRIVIAL(debug) << "Max path length is calculated as " << std::to_string(sub_lmp_as_string.length())
                                 << " + (" << std::to_string(interval.length) << " * 5) = "
                                 << std::to_string(max_path_length) << "\n";

        if (g_local_assembly_kmer_size > sub_lmp_as_string.length()) {
            BOOST_LOG_TRIVIAL(warning) << "Local assembly kmer size " << std::to_string(g_local_assembly_kmer_size)
                                       << " is greater than the length of the interval string "
                                       << std::to_string(sub_lmp_as_string.length())
                                       << ". Skipping local assembly for "
                                       << filepath << "\n";
        } else {
            // get start and end kmer from sub_lmp path
            auto start_kmer{sub_lmp_as_string.substr(0, g_local_assembly_kmer_size)};
            auto end_kmer{sub_lmp_as_string.substr(sub_lmp_as_string.size() - g_local_assembly_kmer_size)};
            // create outpath for local assembly file
            const auto out_path = filepath.substr(0, filepath.rfind('.')) +
                                  "_local_assembly_K" +
                                  std::to_string(g_local_assembly_kmer_size) + ".fa";

            const auto len_threshold{150};
            if (max_path_length < len_threshold) {
                // run local assembly
                BOOST_LOG_TRIVIAL(info) << now() << " Running local assembly for "
                                        << pnode->get_name() + "." + to_string(interval.start) + "-" +
                                           to_string(interval.get_end())
                                        << "\n";

                local_assembly(filepath, start_kmer, end_kmer, out_path, g_local_assembly_kmer_size, max_path_length);

                BOOST_LOG_TRIVIAL(info) << now() << " Finished local assembly for "
                                        << pnode->get_name() + "." + to_string(interval.start) + "-" +
                                           to_string(interval.get_end())
                                        << "\n";
            } else {
                BOOST_LOG_TRIVIAL(debug) << "Max path length " << std::to_string(max_path_length) << " is greater than "
                                         << std::to_string(len_threshold) << ". Skipping local assembly.";
            }
        }
        read_overlap_coordinates.clear();
        fa.clear();
    }

    readfile.close();
}


Interval apply_buffer_to_interval(const Interval &interval, const int32_t &buff) {
    uint32_t start;
    if (buff < interval.start) {
        start = interval.start - buff;
    } else {  // start_diff <= 0:
        start = 0;
    }
    const uint32_t end{interval.get_end() + buff};
    Interval padded_interval{start, end};
    return padded_interval;
}