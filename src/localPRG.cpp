#include <fstream>
#include <set>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <boost/log/trivial.hpp>

#include "minimizer.h"
#include "localPRG.h"
#include "inthash.h"
#include "utils.h"
#include "fastaq.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


LocalPRG::LocalPRG(uint32_t i, const std::string &n, const std::string &p) : next_id(0), buff(" "), next_site(5), id(i),
                                                                             name(n),
                                                                             seq(p), num_hits(2, 0) {
    std::vector<uint32_t> v;
    // avoid error if a prg contains only empty space as it's sequence
    if (seq.find_first_not_of("\t\n\v\f\r") != std::string::npos) {
        std::vector<uint32_t> b = build_graph(Interval(0, seq.size()), v);
    } else {
        prg.add_node(0, "", Interval(0, 0));
    }
}

bool LocalPRG::isalpha_string(const std::string &s) const {
    // Returns if a string s is entirely alphabetic
    for (char j : s)
        if (isalpha(j) == 0) {
            //std::cout << "Found non-alpha char: " << s[j] << std::endl;
            return false;
        }
    return true;
}

std::string LocalPRG::string_along_path(const Path &p) const {
    //std::cout << p << std::endl;
    assert(p.get_start() <= seq.length());
    assert(p.get_end() <= seq.length());
    std::string s;
    for (const auto &it : p.path) {
        s += seq.substr(it.start, it.length);
        //std::cout << s << std::endl;
    }
    //std::cout << "lengths: " << s.length() << ", " << p.length << std::endl;
    assert(s.length() == p.length() ||
           assert_msg("sequence length " << s.length() << " is not equal to path length " << p.length()));
    return s;
}

std::string LocalPRG::string_along_path(const std::vector<LocalNodePtr> &p) {
    std::string s;
    for (const auto &n: p) {
        s += n->seq;
    }
    return s;
}

std::vector<LocalNodePtr> LocalPRG::nodes_along_path(const Path &p) const {
    std::vector<LocalNodePtr> path_nodes;
    path_nodes.reserve(100);

    for (auto it = p.path.begin(); it != p.path.end(); ++it) {
        const auto &interval = *it;
        uint32_t interval_end = it->get_end();

        // find the appropriate node of the prg
        for (auto n = prg.nodes.begin(); n != prg.nodes.end(); ++n) {
            auto node_ptr = n->second;
            uint32_t node_start = node_ptr->pos.start;
            uint32_t node_end = node_ptr->pos.get_end();

            bool is_overlapping_node = interval_end > node_start and interval.start < node_end;
            bool is_equal_node = interval.start == node_start and interval_end == node_end;
            bool is_null_node = (
                    interval.start == node_start
                    and interval.length == 0
                    and it == --(p.path.end())
                    and n != prg.nodes.begin()
            );

            if (is_overlapping_node or is_equal_node or is_null_node) {
                path_nodes.push_back(n->second);
            } else if (interval_end < node_start) {
                // local nodes are labelled in order of occurance in the linear prg string, no need to search further
                break;
            }
        }
    }
    return path_nodes;
}

std::vector<Interval> LocalPRG::split_by_site(const Interval &i) const {
    // Splits interval by next_site based on substring of seq in the interval
    //std::cout << "splitting by site " << next_site << " in interval " << i << std::endl;

    // Split first by var site
    std::vector<Interval> v;
    v.reserve(4);
    std::string::size_type k = i.start;
    std::string d = buff + std::to_string(next_site) + buff;
    std::string::size_type j = seq.find(d, k);
    while (j != std::string::npos and j + d.size() <= i.get_end()) {
        v.emplace_back(Interval(k, j));
        k = j + d.size();
        j = seq.find(d, k);
    }

    if (j != std::string::npos and j < i.get_end() and j + d.size() > i.get_end()) {
        v.emplace_back(Interval(k, j));
    } else if (j != std::string::npos and j + d.size() == i.get_end()) {
        v.emplace_back(Interval(k, j));
        if (seq.find(buff, j + d.size()) == j + d.size()) {
            v.emplace_back(Interval(j + d.size(), j + d.size()));
        }
    } else {
        v.emplace_back(Interval(k, i.get_end()));
    }

    assert(v[0].start >= i.start);
    for (uint32_t l = 1; l != v.size(); ++l) {
        assert(v[l - 1].get_end() <= v[l].start || assert_msg(
                v[l - 1].get_end() << ">" << v[l].start << " giving overlapping intervals  " << v[l - 1] << " and "
                                   << v[l]));
    }
    assert(v.back().get_end() <= i.get_end());

    // then split by var site + 1
    std::vector<Interval> w;
    w.reserve(20);
    d = buff + std::to_string(next_site + 1) + buff;
    for (uint32_t l = 0; l != v.size(); ++l) {
        k = v[l].start;
        j = seq.find(d, k);
        while (j != std::string::npos and j + d.size() <= v[l].get_end()) {
            w.emplace_back(Interval(k, j));
            k = j + d.size();
            j = seq.find(d, k);
        }
        if (j != std::string::npos and j < v[l].get_end() and j + d.size() > v[l].get_end()) {
            w.emplace_back(Interval(k, j));
        } else if (j != std::string::npos and j + d.size() == v[l].get_end()) {
            w.emplace_back(Interval(k, j));
            if (seq.find(buff, j + d.size()) == j + d.size()) {
                v.emplace_back(Interval(j + d.size(), j + d.size()));
            }
        } else {
            w.emplace_back(Interval(k, v[l].get_end()));
        }
    }
    if (v.size() == w.size() && v.size() == 3) {
        BOOST_LOG_TRIVIAL(warning) << "There was something dodgy with var site " << next_site << ": found no separated "
                    "alternates. I'm going to assume for now that this is as a result of straggly ends of sequences "
                    "which don't align nicely, but you should check this. To handle, add an empty interval alternate.";
        std::vector<Interval> x;
        for (uint32_t l = 0; l != w.size() - 1; ++l) {
            x.push_back(w[l]);
        }
        x.emplace_back(Interval(w[w.size() - 2].get_end(), w[w.size() - 2].get_end()));
        x.push_back(w[w.size() - 1]);
        w = x;
    }

    assert(w[0].start >= i.start);
    for (uint32_t l = 1; l != w.size(); ++l) {
        assert(w[l - 1].get_end() <= w[l].start || assert_msg(
                w[l - 1].get_end() << ">" << w[l].start << " giving overlapping intervals  " << w[l - 1] << " and "
                                   << w[l]
                                   << " when splitting seq :" << seq.substr(i.start, i.length)));
    }
    assert(w.back().get_end() <= i.get_end());
    return w;
}

std::vector<uint32_t> LocalPRG::build_graph(const Interval &i,
                                            const std::vector<uint32_t> &from_ids,
                                            uint32_t current_level) {
    // we will return the ids on the ends of any stretches of graph added
    std::vector<uint32_t> end_ids;
    end_ids.reserve(20);

    // save the start id, so can add 0, and the last id to the index at level 0 at the end
    uint32_t start_id = next_id;

    // add nodes
    std::string s = seq.substr(i.start, i.length); //check length correct with this end...
    if (isalpha_string(s)) // should return true for empty string too
    {
        prg.add_node(next_id, s, i);
        // add edges from previous part of graph to start of this interval
        for (uint32_t j = 0; j != from_ids.size(); j++) {
            prg.add_edge(from_ids[j], next_id);
        }
        end_ids.push_back(next_id);
        next_id++;
    } else {
        // split by next var site
        std::vector<Interval> v = split_by_site(i); // should have length at least 4
        if (v.size() < (uint32_t) 4) {
            BOOST_LOG_TRIVIAL(warning) << "In conversion from linear localPRG string to graph, splitting the string by "
                        "the next var site resulted in the wrong number of intervals. Please check that site numbers "
                        "are flanked by a space on either side. Or perhaps ordering of numbers in GFA is irregular?!"
                        "Size of partition based on site " << next_site << " is " << v.size();
            std::exit(-1);
        }
        next_site += 2;
        // add first interval (should be alpha)
        s = seq.substr(v[0].start, v[0].length);
        if (!(isalpha_string(s))) {
            BOOST_LOG_TRIVIAL(warning) << "In conversion from linear localPRG string to graph, splitting the string by "
                        "the next var site resulted in the first interval being non alphabetic. Please check that site "
                        "numbers are flanked by a space on either side. Or perhaps ordering of numbers in GFA is "
                        "irregular?! After splitting by site " << next_site << " do not have alphabetic sequence before "
                        "var site: " << v[0];
            std::exit(-1);
        }
        prg.add_node(next_id, s, v[0]);
        // add edges from previous part of graph to start of this interval
        for (uint32_t j = 0; j != from_ids.size(); j++) {
            prg.add_edge(from_ids[j], next_id);
        }

        std::vector<uint32_t> mid_ids;
        mid_ids.reserve(20);
        mid_ids.push_back(next_id);
        next_id++;
        // add (recurring as necessary) middle intervals
        for (uint32_t j = 1; j != v.size() - 1; j++) {
            std::vector<uint32_t> w = build_graph(v[j], mid_ids, current_level + 1);
            end_ids.insert(end_ids.end(), w.begin(), w.end());
        }
        // add end interval
        end_ids = build_graph(v.back(), end_ids, current_level);
    }
    if (start_id == 0) {
        assert(end_ids.size() == 1);
    }
    return end_ids;
}

std::vector<Path> LocalPRG::shift(Path p) const {
    // returns all paths of the same length which have been shifted by one position along prg graph
    Path q;
    q = p.subpath(1, p.length() - 1);
    std::vector<LocalNodePtr> n;
    std::vector<Path> return_paths;
    std::deque<Path> short_paths = {q};
    std::vector<Path> k_paths;
    bool non_terminus;
    //uint exp_num_return_seqs = 0;

    // first find extensions of the path
    while (!short_paths.empty()) {
        p = short_paths.front();
        n = nodes_along_path(p);
        short_paths.pop_front();

        // if we can extend within the same localnode, do
        if (p.get_end() < n.back()->pos.get_end()) {
            p.path.back().length += 1;
            k_paths.push_back(p);
        } else if (p.get_end() != (--(prg.nodes.end()))->second->pos.get_end()) {
            for (uint32_t i = 0; i != n.back()->outNodes.size(); ++i) {
                //exp_num_return_seqs += 1;
                short_paths.push_back(p);
                short_paths.back().add_end_interval(
                        Interval(n.back()->outNodes[i]->pos.start, n.back()->outNodes[i]->pos.start));
            }
        }
    }

    // now check if by adding null nodes we reach the end of the prg
    for (uint32_t i = 0; i != k_paths.size(); ++i) {
        short_paths = {k_paths[i]};
        non_terminus = false; // assume there all extensions are terminal i.e. reach end or prg

        while (!short_paths.empty()) {
            p = short_paths.front();
            n = nodes_along_path(p);
            short_paths.pop_front();

            if (n.back()->pos.get_end() == (--(prg.nodes.end()))->second->pos.get_end()) {
                return_paths.push_back(p);
            } else if (n.back()->pos.get_end() == p.get_end()) {
                for (uint32_t j = 0; j != n.back()->outNodes.size(); ++j) {
                    if (n.back()->outNodes[j]->pos.length == 0) {
                        short_paths.push_back(p);
                        short_paths.back().add_end_interval(n.back()->outNodes[j]->pos);
                    } else {
                        non_terminus = true;
                    }
                }
            } else {
                non_terminus = true;
            }
        }
        if (non_terminus) {
            return_paths.push_back(k_paths[i]);
        }
    }

    return return_paths;
}

void LocalPRG::minimizer_sketch(std::shared_ptr<Index> index, const uint32_t w, const uint32_t k) {
    BOOST_LOG_TRIVIAL(info) << "Sketch PRG " << name << " which has " << prg.nodes.size() << " nodes";

    // clean up after any previous runs
    // although note we can't clear the index because it is also added to by other LocalPRGs
    kmer_prg.clear();

    // declare variables
    std::vector<Path> walk_paths, shift_paths, v;
    walk_paths.reserve(100);
    shift_paths.reserve(100);
    std::deque<KmerNodePtr> current_leaves, end_leaves;
    std::deque<std::vector<Path>> shifts;
    std::deque<Interval> d;
    Path kmer_path;
    std::string kmer;
    uint64_t smallest;
    std::pair<uint64_t, uint64_t> kh;
    KmerHash hash;
    uint32_t num_kmers_added = 0;
    KmerNodePtr kn, new_kn;
    std::vector<LocalNodePtr> n;
    bool mini_found_in_window;
    size_t num_AT = 0;

    // create a null start node in the kmer graph
    d = {Interval(0, 0)};
    kmer_path.initialize(d);
    kmer_prg.add_node(kmer_path);
    num_kmers_added += 1;

    // if this is a null prg, return the null kmergraph
    if (prg.nodes.size() == 1 and prg.nodes[0]->pos.length < k) {
        return;
    }

    // find first w,k minimizers
    walk_paths = prg.walk(prg.nodes.begin()->second->id, 0, w + k - 1);
    for (uint32_t i = 0; i != walk_paths.size(); ++i) {
        // find minimizer for this path 
        smallest = std::numeric_limits<uint64_t>::max();
        mini_found_in_window = false;
        for (uint32_t j = 0; j != w; j++) {
            kmer_path = walk_paths[i].subpath(j, k);
            if (!kmer_path.path.empty()) {
                kmer = string_along_path(kmer_path);
                kh = hash.kmerhash(kmer, k);
                smallest = std::min(smallest, std::min(kh.first, kh.second));
            }
        }
        for (uint32_t j = 0; j != w; j++) {
            kmer_path = walk_paths[i].subpath(j, k);
            if (!kmer_path.path.empty()) {
                kmer = string_along_path(kmer_path);
                kh = hash.kmerhash(kmer, k);
                n = nodes_along_path(kmer_path);

                if (prg.walk(n.back()->id, n.back()->pos.get_end(), w + k - 1).empty()) {
                    while (kmer_path.get_end() >= n.back()->pos.get_end() and n.back()->outNodes.size() == 1 and
                           n.back()->outNodes[0]->pos.length == 0) {
                        kmer_path.add_end_interval(n.back()->outNodes[0]->pos);
                        n.push_back(n.back()->outNodes[0]);
                    }
                }

                if (kh.first == smallest or kh.second == smallest) {
                    const auto found = find_if(kmer_prg.nodes.begin(), kmer_prg.nodes.end(), condition(kmer_path));
                    if (found == kmer_prg.nodes.end()) {
                        // add to index, kmer_prg
                        num_AT = std::count(kmer.begin(), kmer.end(), 'A') + std::count(kmer.begin(), kmer.end(), 'T');
                        kn = kmer_prg.add_node_with_kh(kmer_path, std::min(kh.first, kh.second), num_AT);
                        index->add_record(std::min(kh.first, kh.second), id, kmer_path, kn->id,
                                          (kh.first <= kh.second));
                        num_kmers_added += 1;
                        if (!mini_found_in_window) {
                            kmer_prg.add_edge(kmer_prg.nodes[0], kn);
                        }
                        mini_found_in_window = true;
                        current_leaves.push_back(kn);
                    }
                }
            }
        }
    }

    // while we have intermediate leaves of the kmergraph, for each in turn, explore the neighbourhood
    // in the prg to find the next minikmers as you walk the prg
    while (!current_leaves.empty()) {
        kn = current_leaves.front();
        current_leaves.pop_front();
        assert(kn->khash < std::numeric_limits<uint64_t>::max());

        // find all paths which are this kmernode shifted by one place along the graph
        shift_paths = shift(kn->path);
        if (shift_paths.empty()) {
            //assert(kn->path.get_start() == 0); not true for a too short test, would be true if all paths long enough to have at least 2 minikmers on...
            end_leaves.push_back(kn);
        }
        for (uint32_t i = 0; i != shift_paths.size(); ++i) {
            v = {shift_paths[i]};
            shifts.push_back(v);
        }
        shift_paths.clear();

        while (!shifts.empty()) {
            v = shifts.front();
            shifts.pop_front();
            assert(v.back().length() == k);
            kmer = string_along_path(v.back());
            kh = hash.kmerhash(kmer, k);
            if (std::min(kh.first, kh.second) <= kn->khash) {
                // found next minimizer
                const auto found = find_if(kmer_prg.nodes.begin(), kmer_prg.nodes.end(), condition(v.back()));
                if (found == kmer_prg.nodes.end()) {
                    num_AT = std::count(kmer.begin(), kmer.end(), 'A') + std::count(kmer.begin(), kmer.end(), 'T');
                    new_kn = kmer_prg.add_node_with_kh(v.back(), std::min(kh.first, kh.second), num_AT);
                    index->add_record(std::min(kh.first, kh.second), id, v.back(), new_kn->id, (kh.first <= kh.second));
                    kmer_prg.add_edge(kn, new_kn);
                    if (v.back().get_end() == (--(prg.nodes.end()))->second->pos.get_end()) {
                        end_leaves.push_back(new_kn);
                    } else if (find(current_leaves.begin(), current_leaves.end(), new_kn) == current_leaves.end()) {
                        current_leaves.push_back(new_kn);
                    }
                    num_kmers_added += 1;
                } else {
                    kmer_prg.add_edge(kn, *found);
                    if (v.back().get_end() == (--(prg.nodes.end()))->second->pos.get_end()) {
                        end_leaves.push_back(*found);
                    } else if (find(current_leaves.begin(), current_leaves.end(), *found) ==
                               current_leaves.end()) {
                        current_leaves.push_back(*found);
                    }
                }
            } else if (v.size() == w) {
                // the old minimizer has dropped out the window, minimizer the w new kmers
                smallest = std::numeric_limits<uint64_t>::max();
                mini_found_in_window = false;
                for (uint32_t j = 0; j != w; j++) {
                    kmer = string_along_path(v[j]);
                    kh = hash.kmerhash(kmer, k);
                    smallest = std::min(smallest, std::min(kh.first, kh.second));
                    //std::cout << min(kh.first, kh.second) << " ";
                }
                for (uint32_t j = 0; j != w; j++) {
                    kmer = string_along_path(v[j]);
                    kh = hash.kmerhash(kmer, k);
                    if (kh.first == smallest or kh.second == smallest) {
                        const auto found = find_if(kmer_prg.nodes.begin(), kmer_prg.nodes.end(), condition(v[j]));
                        if (found == kmer_prg.nodes.end()) {
                            num_AT = std::count(kmer.begin(), kmer.end(), 'A') +
                                     std::count(kmer.begin(), kmer.end(), 'T');
                            new_kn = kmer_prg.add_node_with_kh(v[j], std::min(kh.first, kh.second), num_AT);
                            index->add_record(std::min(kh.first, kh.second), id, v[j], new_kn->id,
                                              (kh.first <= kh.second));

                            // if there is more than one mini in the window, edge should go to the first, and from the first to the second
                            if (!mini_found_in_window) {
                                kmer_prg.add_edge(kn, new_kn);
                            }
                            mini_found_in_window = true;

                            if (v.back().get_end() == (--(prg.nodes.end()))->second->pos.get_end()) {
                                end_leaves.push_back(new_kn);
                            } else if (find(current_leaves.begin(), current_leaves.end(), new_kn) ==
                                       current_leaves.end()) {
                                current_leaves.push_back(new_kn);
                            }
                            num_kmers_added += 1;
                        } else {
                            if (!mini_found_in_window) {
                                kmer_prg.add_edge(kn, *found);
                            }
                            mini_found_in_window = true;
                            if (v.back().get_end() == (--(prg.nodes.end()))->second->pos.get_end()) {
                                end_leaves.push_back(*found);
                            } else if (find(current_leaves.begin(), current_leaves.end(), *found) ==
                                       current_leaves.end()) {
                                current_leaves.push_back(*found);
                            }
                        }
                    }
                }
            } else if (v.back().get_end() == (--(prg.nodes.end()))->second->pos.get_end()) {
                end_leaves.push_back(kn);
            } else {
                shift_paths = shift(v.back());
                for (uint32_t i = 0; i != shift_paths.size(); ++i) {
                    shifts.push_back(v);
                    shifts.back().push_back(shift_paths[i]);
                }
                shift_paths.clear();
            }
        }
    }

    // create a null end node, and for each end leaf add an edge to this terminus
    assert(!end_leaves.empty());
    d = {Interval((--(prg.nodes.end()))->second->pos.get_end(), (--(prg.nodes.end()))->second->pos.get_end())};
    kmer_path.initialize(d);
    kn = kmer_prg.add_node(kmer_path);
    num_kmers_added += 1;
    for (uint32_t i = 0; i != end_leaves.size(); ++i) {
        kmer_prg.add_edge(end_leaves[i], kn);
    }

    // print, check and return
    assert(num_kmers_added == 0 or kmer_prg.nodes.size() == num_kmers_added ||
           assert_msg("nodes.size(): " << kmer_prg.nodes.size() << " and num minikmers: " << num_kmers_added));
    kmer_prg.remove_shortcut_edges();
    kmer_prg.check();
}

bool intervals_overlap(const Interval &first, const Interval &second) {
    return ((first == second)
            or (second.length == 0 and (first.start == second.start or first.get_end() == second.get_end()))
            or (first.start < second.get_end() and first.get_end() > second.start));
}

std::vector<KmerNodePtr>
LocalPRG::kmernode_path_from_localnode_path(const std::vector<LocalNodePtr> &localnode_path) const {
    std::vector<KmerNodePtr> kmernode_path;

    if (localnode_path.empty())
        return kmernode_path;

    std::deque<Interval> d;
    for (const auto &n: localnode_path) {
        d.push_back(n->pos);
    }

    Path local_path;
    local_path.initialize(d);

    for (const auto &n: kmer_prg.sorted_nodes) {
        for (const auto &interval : local_path.path) {
            if (interval.start > n->path.get_end())
                break;
            else if (interval.get_end() < n->path.get_start())
                continue;
            else if ((intervals_overlap(interval, n->path.path[0]) or intervals_overlap(interval, n->path.path.back()))
                     and not local_path.is_branching(n->path)) {
                //and not n.second->path.is_branching(local_path))
                kmernode_path.push_back(n);
                break;
            }
        }
    }

    assert(!kmernode_path.empty());
    return kmernode_path;
}


std::vector<LocalNodePtr>
LocalPRG::localnode_path_from_kmernode_path(const std::vector<KmerNodePtr> &kmernode_path, const uint32_t w) const {
    BOOST_LOG_TRIVIAL(debug) << "Convert kmernode path to localnode path";
    std::vector<LocalNodePtr> localnode_path, kmernode, walk_path;
    std::vector<Path> walk_paths;
    for (uint32_t i = 0; i != kmernode_path.size(); ++i) {
        if (i != 0 and kmernode_path[i]->path.length() == 0) // only have null paths at beginning and end
        {
            break;
        }
        kmernode = nodes_along_path(kmernode_path[i]->path);

        // if the start of the new localnode path is after the end of the previous, join up WLOG with top path
        while (!localnode_path.empty() and !localnode_path.back()->outNodes.empty() and
               kmernode[0]->id > localnode_path.back()->outNodes[0]->id) {
            localnode_path.push_back(localnode_path.back()->outNodes[0]);
        }
        // if new the localnodes in kmernode overlap old ones, truncate localnode_path
        while (!localnode_path.empty() and kmernode[0]->id <= localnode_path.back()->id) {
            localnode_path.pop_back();
        }
        localnode_path.insert(localnode_path.end(), kmernode.begin(), kmernode.end());
    }

    // extend to beginning of graph if possible
    bool overlap;
    if (localnode_path[0]->id != 0) {
        walk_paths = prg.walk(0, 0, w);
        for (uint32_t i = 0; i != walk_paths.size(); ++i) {
            walk_path = nodes_along_path(walk_paths[i]);
            // does it overlap
            uint32_t n = 0, m = 0;
            overlap = false;
            for (uint32_t j = 0; j != walk_path.size(); ++j) {
                if (walk_path[j] == localnode_path[n]) {
                    if (!overlap) {
                        m = j;
                    }
                    overlap = true;
                    if (n + 1 >= localnode_path.size()) {
                        break;
                    } else {
                        ++n;
                    }
                } else if (overlap) {
                    overlap = false;
                    break;
                }
            }
            if (overlap) {
                localnode_path.insert(localnode_path.begin(), walk_path.begin(), walk_path.begin() + m);
                break;
            }
        }
        if (localnode_path[0]->id != 0) {
            //std::cout << "localnode path still does not start with 0" << std::endl;
            // add the first path to start
            LocalNodePtr next = nullptr;
            while (localnode_path[0]->id != 0 and next != localnode_path[0]) {
                //std::cout << "look for a previous node to " << *localnode_path[0] << std::endl;
                next = prg.get_previous_node(localnode_path[0]);
                if (next != nullptr) {
                    //std::cout << "add " << *next << std::endl;
                    localnode_path.insert(localnode_path.begin(), next);
                }
                //std::cout << "new start node is " << *localnode_path[0] << std::endl;
            }
        }
    }

    // extend to end of graph if possible
    if (localnode_path.back()->id != prg.nodes.size() - 1) {
        walk_paths = prg.walk_back(prg.nodes.size() - 1, seq.length(), w);
        for (uint32_t i = 0; i != walk_paths.size(); ++i) {
            walk_path = nodes_along_path(walk_paths[i]);

            // does it overlap
            uint32_t n = localnode_path.size();
            uint32_t m = 0;
            overlap = false;
            for (uint32_t j = walk_path.size(); j != 0; --j) {
                if (walk_path[j - 1] == localnode_path[n - 1]) {
                    if (!overlap) {
                        m = j;
                    }
                    overlap = true;
                    if (n - 1 == 0) {
                        break;
                    } else {
                        --n;
                    }
                } else if (overlap) {
                    overlap = false;
                    break;
                }
            }
            if (overlap) {
                localnode_path.insert(localnode_path.end(), walk_path.begin() + m, walk_path.end());
                break;
            }
        }
        if (localnode_path.back()->id != prg.nodes.size() - 1) {
            //std::cout << "localnode path still does not end with " << prg.nodes.size()-1 << std::endl;
            // add the first path to end
            while (localnode_path.back()->id != prg.nodes.size() - 1 and !localnode_path.back()->outNodes.empty()) {
                //std::cout << "extend " << *localnode_path.back();
                localnode_path.push_back(localnode_path.back()->outNodes[0]);
                //std::cout << " with " << *localnode_path.back() << std::endl;
            }
        }
    }
    return localnode_path;
}

std::vector<uint32_t> get_covgs_along_localnode_path(const PanNodePtr pnode,
                                                     const std::vector<LocalNodePtr> &localnode_path,
                                                     const std::vector<KmerNodePtr> &kmernode_path,
                                                     const uint32_t &sample_id) {
    // defines estimated per base coverage for the bases of localnode_path based on the coverages from the
    // kmernode_path kmers
    // NB the kmernode_path may be a vector of pointers to nodes in the localPRG copy of the kmergraph,
    // and it is the copy of the kmergraph in the pangraph pnode which stores coverages, which is why we
    // need both

    //define 0 coverage for each base in localnode path
    std::vector<std::vector<uint32_t>> coverages;
    for (const auto &n: localnode_path) {
        std::vector<uint32_t> v(n->pos.length, 0);
        coverages.push_back(v);
    }

    // collect covgs
    uint32_t j = 0, k = 0, start = 0, end = 0;
    for (const auto &kmernode_ptr: kmernode_path) {
        //std::cout << kmernode_ptr->path << std::endl;

        if (kmernode_ptr->path.length() == 0)
            continue;

        while (j < localnode_path.size() and localnode_path[j]->pos.get_end() < kmernode_ptr->path.get_start()) {
            //std::cout << localnode_path[j]->pos.get_end() << " < " << kmernode_ptr->path.get_start() << std::endl;
            j++;
            //std::cout << "j = " << j << std::endl;
        }
        k = j;
        for (const auto &interval : kmernode_ptr->path.path) {
            //std::cout << "k = " << k << " and " << localnode_path[k]->pos << " ~ " << interval << std::endl;
            assert(localnode_path[k]->pos.start <= interval.start
                   and localnode_path[k]->pos.get_end() >= interval.get_end());
            start = interval.start - localnode_path[k]->pos.start;
            end = std::min(start + interval.length, localnode_path[k]->pos.get_end());
            //std::cout << "add from " << start << " to " << end << std::endl;
            for (uint32_t l = start; l < end; ++l) {
                assert(kmernode_ptr->id < pnode->kmer_prg.nodes.size() and
                       pnode->kmer_prg.nodes[kmernode_ptr->id] != nullptr);
                coverages[k][l] = std::max(coverages[k][l],
                                           pnode->kmer_prg.nodes[kmernode_ptr->id]->get_covg(0, sample_id) +
                                           pnode->kmer_prg.nodes[kmernode_ptr->id]->get_covg(1, sample_id));
            }
            k++;
        }
    }

    std::vector<uint32_t> return_coverages;
    for (const auto &v : coverages) {
        return_coverages.insert(return_coverages.end(), v.begin(), v.end());
    }

    return return_coverages;
}

void LocalPRG::write_covgs_to_file(const boost::filesystem::path &filepath, const std::vector<uint32_t> &covgs) const {
    std::ofstream handle;
    handle.open(filepath.string());
    assert (!handle.fail() or assert_msg("Could not open file " << filepath.string()));

    handle << ">" << name << std::endl;
    for (const auto &i : covgs) {
        handle << i << " ";
    }
    handle << std::endl;

    handle.close();
}

void LocalPRG::write_path_to_fasta(const boost::filesystem::path &filepath,
                                   const std::vector<LocalNodePtr> &lmp,
                                   const float &ppath) const {
    std::ofstream handle;
    handle.open(filepath.string());
    assert (!handle.fail() or assert_msg("Could not open file " << filepath.string()));

    handle << ">" << name << "\tlog P(data|sequence)=" << ppath << std::endl;
    for (uint32_t j = 0; j != lmp.size(); ++j) {
        handle << lmp[j]->seq;
    }
    handle << std::endl;

    handle.close();
}

void LocalPRG::append_path_to_fasta(const boost::filesystem::path &filepath,
                                    const std::vector<LocalNodePtr> &lmp,
                                    const float &ppath) const {
    std::ofstream handle;
    handle.open(filepath.string(), std::ios::app);
    assert (!handle.fail() or assert_msg("Could not open file " << filepath.string()));

    handle << ">" << name << "\tlog P(data|sequence)=" << ppath << std::endl;
    for (uint32_t j = 0; j != lmp.size(); ++j) {
        handle << lmp[j]->seq;
    }
    handle << std::endl;

    handle.close();
}

void LocalPRG::write_aligned_path_to_fasta(const boost::filesystem::path &filepath,
                                           const std::vector<LocalNodePtr> &lmp,
                                           const float &ppath) const {
    std::ofstream handle;
    handle.open(filepath.string());
    assert (!handle.fail() or assert_msg("Could not open file " << filepath.string()));

    handle << ">" << name << "\tlog P(data|sequence)=" << ppath << std::endl;

    uint32_t i = 0;
    for (const auto &c : prg.nodes) {
        if (c.second == lmp[i]) {
            handle << lmp[i]->seq;
            i++;
        } else {
            std::string s(c.second->seq.length(), '-'); // s == "------"
            handle << s;
        }
    }
    handle << std::endl;

    handle.close();
}

void LocalPRG::build_vcf(VCF &vcf, const std::vector<LocalNodePtr> &ref) const {
    BOOST_LOG_TRIVIAL(info) << "Build VCF for prg " << name;
    assert(!prg.nodes.empty()); //otherwise empty nodes -> segfault

    std::vector<LocalNodePtr> varpath;
    varpath.reserve(100);
    std::vector<LocalNodePtr> bottompath;
    bottompath.reserve(100);
    uint32_t ref_i = 0;
    auto ref_length = string_along_path(ref).length();

    std::deque<std::vector<LocalNodePtr>> paths;

    std::vector<std::vector<LocalNodePtr>> alts;
    alts.reserve(100);

    std::vector<uint32_t> level_start;

    int level = 0, max_level = 0;
    uint32_t pos;
    std::string vartype = "GRAPHTYPE=SIMPLE";
    std::string ref_seq;
    std::string alt_seq;

    // simple case
    if (ref.size() == 1) {
        return;
    }

    while (ref_i < ref.size() - 1) {
        // first update the level we are at
        if (ref[ref_i]->outNodes.size() > 1) {
            level += 1;
            max_level = std::max(level, max_level);
            level_start.push_back(ref_i);

            if (level > 1) {
                vartype = "GRAPHTYPE=NESTED";
            }
        } else {
            // we have come down a level, add the alts compared to this region
            level -= 1;
            assert(level >= 0);
            assert(!level_start.empty());
            // define ref and pos
            pos = 0;
            ref_seq = "";
            for (uint32_t j = 0; j <= level_start.back(); ++j) {
                pos += ref[j]->seq.length();
            }
            for (uint32_t j = level_start.back() + 1; j <= ref_i; ++j) {
                ref_seq += ref[j]->seq;
            }

            // initialise alt paths
            for (uint32_t n = 0; n < ref[level_start.back()]->outNodes.size(); ++n) {
                if (ref[level_start.back()]->outNodes[n] != ref[level_start.back() + 1]) {
                    varpath = {ref[level_start.back()]->outNodes[n]};
                    paths.push_back(varpath);
                }
            }

            // extend alt paths to end of site
            while (!paths.empty()) {
                varpath = paths[0];
                paths.pop_front();
                if (varpath.back()->outNodes[0]->id == ref[ref_i]->outNodes[0]->id) {
                    alts.push_back(varpath);
                } else {
                    for (uint32_t j = 0; j != varpath.back()->outNodes.size(); ++j) {
                        paths.push_back(varpath);
                        paths.back().push_back(varpath.back()->outNodes[j]);
                    }
                }

                // if have too many alts, just give bottom path and top path
                if (paths.size() > 1000) {
                    paths.clear();
                    alts.clear();
                    bottompath.push_back(ref[level_start.back()]->outNodes.back());
                    while (!bottompath.back()->outNodes.empty() and
                           bottompath.back()->outNodes[0]->id != ref[ref_i]->outNodes[0]->id) {
                        bottompath.push_back(bottompath.back()->outNodes.back());
                    }
                    alts.push_back(bottompath);

                    bottompath.clear();
                    bottompath.push_back(ref[level_start.back()]->outNodes[0]);
                    while (!bottompath.back()->outNodes.empty() and
                           bottompath.back()->outNodes[0]->id != ref[ref_i]->outNodes[0]->id) {
                        bottompath.push_back(bottompath.back()->outNodes[0]);
                    }
                    alts.push_back(bottompath);
                    bottompath.clear();

                    vartype = "GRAPHTYPE=TOO_MANY_ALTS";
                    break;
                }
            }

            // add sites to vcf
            assert(pos + ref_seq.length() <= ref_length);
            for (auto &alt : alts) {
                for (auto &j : alt) {
                    alt_seq += j->seq;
                }
                if (ref_seq != alt_seq) {
                    vcf.add_record(name, pos, ref_seq, alt_seq, ".", vartype);
                }
                alt_seq = "";
            }
            alts.clear();

            level_start.pop_back();
            if (level == 0) {
                assert(level_start.empty());
                vartype = "GRAPHTYPE=SIMPLE";
            }
        }
        ref_i++;
    }
}

void
LocalPRG::add_sample_gt_to_vcf(VCF &vcf,
                               const std::vector<LocalNodePtr> &rpath,
                               const std::vector<LocalNodePtr> &sample_path,
                               const std::string &sample_name) const {
    BOOST_LOG_TRIVIAL(info) << "Update VCF with sample path";
    /*for (uint i=0; i!=sample_path.size(); ++i)
    {
	std::cout << *sample_path[i] << " ";
    }
    std::cout << std::endl;
    std::cout << now() << "Using ref path" << std::endl;
    for (uint i=0; i!=rpath.size(); ++i)
    {
        std::cout << *rpath[i] << " ";
    }
    std::cout << std::endl;*/
    assert(!prg.nodes.empty()); //otherwise empty nodes -> segfault

    // if prg has only one node, simple case
    if (prg.nodes.size() == 1) {
        vcf.samples.push_back(sample_name);
    }

    std::vector<LocalNodePtr> refpath, samplepath;
    refpath.reserve(100);
    assert(not rpath.empty());
    refpath.push_back(rpath[0]);
    samplepath.reserve(100);
    assert(not sample_path.empty());
    samplepath.push_back(sample_path[0]);
    uint32_t ref_i = 1, sample_id = 1, pos = 0, pos_to = 0;
    std::vector<uint32_t> sample_covg(6, 0);
    std::string ref;
    std::string alt;
    bool found_new_site = false;

    while (!refpath.back()->outNodes.empty() or refpath.size() > 1) {
        if (refpath.back()->id < samplepath.back()->id) {
            assert(rpath.size() > ref_i);
            refpath.push_back(rpath[ref_i]);
            found_new_site = true;
            ref_i++;
        } else if (samplepath.back()->id < refpath.back()->id) {
            assert(sample_path.size() > sample_id);
            samplepath.push_back(sample_path[sample_id]);
            found_new_site = true;
            sample_id++;
        } else if (found_new_site) {
            // refpath back == samplepath back
            // add ref allele from previous site to this one
            //std::cout << "update with ref alleles from " << pos << " to " << pos_to << std::endl;
            vcf.add_sample_ref_alleles(sample_name, name, pos, pos_to);
            pos = pos_to;

            // add new site to vcf
            //std::cout << "find ref seq" << std::endl;
            for (uint32_t j = 1; j < refpath.size() - 1; ++j) {
                ref += refpath[j]->seq;
                //std::cout << ref << std::endl;
            }
            //std::cout << "find alt seq" << std::endl;
            for (uint32_t j = 1; j < samplepath.size() - 1; ++j) {
                alt += samplepath[j]->seq;
                //std::cout << alt << std::endl;
            }

            //std::cout << "add sample gt" << std::endl;
            vcf.add_sample_gt(sample_name, name, pos, ref, alt);
            found_new_site = false;

            // prepare for next iteration
            for (uint32_t j = 1; j < refpath.size() - 1; ++j) {
                pos += refpath[j]->pos.length;
            }
            refpath.erase(refpath.begin(), refpath.end() - 1);
            if (refpath.back()->id != prg.nodes.size() - 1) {
                ref = "";
                alt = "";
                assert(not refpath.empty());
                pos += refpath.back()->pos.length;
                assert(rpath.size() > ref_i);
                refpath.push_back(rpath[ref_i]);
                ref_i++;
                samplepath.erase(samplepath.begin(), samplepath.end() - 1);
                assert(sample_path.size() > sample_id);
                samplepath.push_back(sample_path[sample_id]);
                sample_id++;
            }
            pos_to = pos;
        } else {
            //std::cout << pos_to << std::endl;
            refpath.erase(refpath.begin(), refpath.end() - 1);
            if (refpath.back()->id != prg.nodes.size() - 1) {
                ref = "";
                alt = "";
                assert(not refpath.empty());
                pos_to += refpath.back()->pos.length;
                assert(rpath.size() > ref_i);
                refpath.push_back(rpath[ref_i]);
                ref_i++;
                samplepath.erase(samplepath.begin(), samplepath.end() - 1);
                assert(sample_path.size() > sample_id);
                samplepath.push_back(sample_path[sample_id]);
                sample_id++;
            }
        }
    }
    //std::cout << "add last ref alleles" << std::endl;
    vcf.add_sample_ref_alleles(sample_name, name, pos, pos_to);
    //std::cout << "done" << std::endl;
}

// Find the path through the PRG which deviates at pos from the ref path with alt sequence
std::vector<LocalNodePtr> LocalPRG::find_alt_path(const std::vector<LocalNodePtr> &ref_path,
                                                  const uint32_t pos,
                                                  const std::string &ref,
                                                  const std::string &alt) const {
    //std::cout << now() << "Find alt path for PRG " << name << " variant " << pos << " " << ref << " " << alt << std::endl;
    std::vector<LocalNodePtr> alt_path, considered_path;
    std::deque<std::vector<LocalNodePtr>> paths_in_progress;
    uint32_t ref_added = 0, pos_along_ref_path = 0;

    std::string working_alt = alt;
    if (alt == ".")
        working_alt = "";
    std::string working_ref = ref;
    if (ref == ".")
        working_ref = "";

    for (const auto &n: ref_path) {
        if (ref_added < pos) {
            ref_added += n->pos.length;
            alt_path.push_back(n);
            pos_along_ref_path++;
        } else {
            break;
        }
    }
    //std::cout << "pos " << pos << " pos_along_ref_path " << pos_along_ref_path << " ref_path.size() " << ref_path.size() << std::endl;

    // find the localnodeptr we want to make our way back to
    while (pos_along_ref_path < ref_path.size() - 1
           and (ref_added < pos + working_ref.length() or ref_path[pos_along_ref_path]->pos.length == 0)) {
        ref_added += ref_path[pos_along_ref_path]->pos.length;
        pos_along_ref_path++;
    }
    assert(pos_along_ref_path < ref_path.size());
    auto ref_node_to_find = ref_path[pos_along_ref_path];
    //std::cout << "trying to find " << *ref_node_to_find;

    // find an alt path with the required sequence
    if (alt_path.empty() and not ref_path.empty() and ref_path[0]->pos.length == 0)
        alt_path.push_back(ref_path[0]);
    assert(!alt_path.empty());
    for (const auto &m: alt_path.back()->outNodes) {
        paths_in_progress.push_back({m});
    }

    while (!paths_in_progress.empty()) {
        considered_path = paths_in_progress.front();
        paths_in_progress.pop_front();

        /*std::cout << "consider path ";
        for (const auto &t : considered_path){
            std::cout << t->pos << " ";
        }
        std::cout << std::endl;*/

        auto considered_seq = string_along_path(considered_path);
        //std::cout << "considered_seq " << considered_seq << std::endl;

        if (considered_seq == working_alt) {
            // check if merge with ref path
            if (find(considered_path.back()->outNodes.begin(), considered_path.back()->outNodes.end(), ref_node_to_find)
                != considered_path.back()->outNodes.end()) {
                alt_path.insert(alt_path.end(), considered_path.begin(), considered_path.end());
                alt_path.insert(alt_path.end(), ref_path.begin() + pos_along_ref_path, ref_path.end());
                /*std::cout << "found alt path ";
                for (const auto &t : considered_path){
                    std::cout << t->pos << " ";
                }
                std::cout << std::endl;*/
                return alt_path;
            } else {
                for (const auto &m: considered_path.back()->outNodes) {
                    paths_in_progress.push_back(considered_path);
                    paths_in_progress.back().push_back(m);
                }
            }
        } else if (considered_seq.length() <= working_alt.length()
                   and considered_seq == working_alt.substr(0, considered_seq.length())) {
            for (const auto &m: considered_path.back()->outNodes) {
                paths_in_progress.push_back(considered_path);
                paths_in_progress.back().push_back(m);
            }
        }
    }
    assert (true or assert_msg("Should have found an alt path!!"));
    return alt_path; // this never happens
}

void LocalPRG::append_kmer_covgs_in_range(const KmerGraph &kg,
                                          const std::vector<KmerNodePtr> &kmer_path,
                                          const std::vector<LocalNodePtr> &local_path,
                                          const uint32_t &pos_from,
                                          const uint32_t &pos_to,
                                          std::vector<uint32_t> &fwd_covgs,
                                          std::vector<uint32_t> &rev_covgs,
                                          const uint32_t &sample_id) const {
    //std::cout << "START" << std::endl;
    assert(fwd_covgs.empty());
    assert(rev_covgs.empty());
    assert(kmer_path.size() > 1);
    //std::cout << "add kmer coverages from " << pos_from  << " to " << pos_to << std::endl;

    // start by adding coverage before we get to first kmer
    uint32_t added = 0, k = 0;
    //std::cout << "first kmer" << kmer_path[1]->path << std::endl;
    for (const auto &n: local_path) {
        //std::cout << n->pos << " ";
        if (n->pos.length == 0) {
            continue;
        } else if (n->pos.get_end() < kmer_path[1]->path.get_start()) {
            added += n->pos.length;
        } else if (n->pos.get_end() >= kmer_path[1]->path.get_start()
                   and n->pos.start < kmer_path[1]->path.get_end()) {
            added += kmer_path[1]->path.get_start() - n->pos.start;
            break;
        }
    }

    //std::cout << "added " << added << std::endl;
    KmerNodePtr prev = nullptr;
    for (const auto &n: kmer_path) {
        //std::cout << n->path << " ";
        if (n->path.length() == 0)
            continue;
        else if (prev != nullptr) {
            auto prev_interval_it = prev->path.path.begin();
            while (prev_interval_it->get_end() < n->path.get_start()) {
                added += prev_interval_it->length;
                prev_interval_it++;
            }
            added += n->path.get_start() - prev_interval_it->start;
        } else {
            k = n->path.length();
            //added += k;
        }

        //std::cout << "pos_from:" << pos_from << " < added + k:" << added + k << " and added: " << added << " < pos_to:" << pos_to << std::endl;

        if (pos_from <= added + k and added < pos_to) {
            //std::cout << " add " << n->path << std::endl;
            assert(n->id < kg.nodes.size() and kg.nodes[n->id] != nullptr);
            fwd_covgs.push_back(kg.nodes.at(n->id)->get_covg(0, sample_id));
            rev_covgs.push_back(kg.nodes.at(n->id)->get_covg(1, sample_id));
        } else if (added > pos_to)
            break;

        prev = n;
    }
}

uint32_t sum(const std::vector<uint32_t> &v) {
    return std::accumulate(v.begin(), v.end(), 0);
}

uint32_t mean(const std::vector<uint32_t> &v) {
    /*std::string s;
    for (const auto &i : v)
        s += int_to_string(i) + " ";
    BOOST_LOG_TRIVIAL(debug) << "mean of " << s;*/
    if (v.empty())
        return 0;
    return std::accumulate(v.begin(), v.end(), 0) / v.size();
}

uint32_t median(std::vector<uint32_t> v) {
    /*std::string s;
    for (const auto &i : v)
        s += int_to_string(i) + " ";
    BOOST_LOG_TRIVIAL(debug) << "median of " << s;*/
    if (v.empty())
        return 0;
    std::sort(v.begin(), v.end());
    if (v.size() % 2 == 1) {
        int n = (v.size() + 1) / 2;
        return v[n - 1];
    } else {
        int n1 = (v.size() + 1) / 2;
        int n2 = (v.size() - 1) / 2;
        return (v[n1] + v[n2]) / 2;
    }
}

uint32_t mode(std::vector<uint32_t> v) {
    std::sort(v.begin(), v.end());
    /*std::string s;
    for (const auto &i : v)
        s += int_to_string(i) + " ";
    BOOST_LOG_TRIVIAL(debug) << "mode of " << s;*/
    uint32_t counter = 1;
    uint32_t max_count = 1;
    uint32_t most_common = 0;
    uint32_t last = 0;
    for (const auto &n : v) {
        if (n == last)
            counter++;
        else {
            if (counter > max_count) {
                max_count = counter;
                most_common = last;
            }
            counter = 1;
        }
        last = n;
    }
    return most_common;
}

void LocalPRG::add_sample_covgs_to_vcf(VCF &vcf, const KmerGraph &kg,
                                       const std::vector<LocalNodePtr> &ref_path,
                                       const std::string &sample_name,
                                       const uint32_t &sample_id) const {
    BOOST_LOG_TRIVIAL(info) << "Update VCF with sample covgs";

    assert(!prg.nodes.empty()); //otherwise empty nodes -> segfault
    vcf.sort_records();

    std::vector<LocalNodePtr> alt_path;

    std::vector<KmerNodePtr> ref_kmer_path = kmernode_path_from_localnode_path(ref_path);
    /*std::cout << "ref path: ";
    for (const auto &n : ref_kmer_path){
        std::cout << n->path << " ";
    }
    std::cout << std::endl;*/

    std::vector<KmerNodePtr> alt_kmer_path;

    std::vector<uint32_t> ref_fwd_covgs;
    std::vector<uint32_t> ref_rev_covgs;
    std::vector<uint32_t> alt_fwd_covgs;
    std::vector<uint32_t> alt_rev_covgs;

    for (auto &record : vcf.records) {
        //std::cout << record << std::endl;
        // find corresponding ref kmers
        auto end_pos = record.pos + record.ref.length();
        if (record.ref == ".")
            end_pos = record.pos;

        append_kmer_covgs_in_range(kg, ref_kmer_path, ref_path, record.pos, end_pos,
                                   ref_fwd_covgs, ref_rev_covgs, sample_id);

        // find corresponding alt kmers
        // if sample has alt path, we have the kmer path for this, but otherwise we will need to work it out
        auto sample_it = find(vcf.samples.begin(), vcf.samples.end(), sample_name);
        assert(sample_it != vcf.samples.end());
        auto sample_index = distance(vcf.samples.begin(), sample_it);
        assert((uint) sample_index != vcf.samples.size());
        assert(record.samples.size() > (uint) sample_index);

        record.samples[sample_index]["MEAN_FWD_COVG"] = {};
        record.samples[sample_index]["MEAN_REV_COVG"] = {};
        record.samples[sample_index]["MED_FWD_COVG"] = {};
        record.samples[sample_index]["MED_REV_COVG"] = {};
        record.samples[sample_index]["SUM_FWD_COVG"] = {};
        record.samples[sample_index]["SUM_REV_COVG"] = {};

        record.samples[sample_index]["MEAN_FWD_COVG"].push_back(mean(ref_fwd_covgs));
        record.samples[sample_index]["MEAN_REV_COVG"].push_back(mean(ref_rev_covgs));
        record.samples[sample_index]["MED_FWD_COVG"].push_back(median(ref_fwd_covgs));
        record.samples[sample_index]["MED_REV_COVG"].push_back(median(ref_rev_covgs));
        record.samples[sample_index]["SUM_FWD_COVG"].push_back(sum(ref_fwd_covgs));
        record.samples[sample_index]["SUM_REV_COVG"].push_back(sum(ref_rev_covgs));

        for (const auto &alt_allele : record.alt) {
            alt_path = find_alt_path(ref_path, record.pos, record.ref, alt_allele);
            alt_kmer_path = kmernode_path_from_localnode_path(alt_path);

            // find alt covgs
            end_pos = record.pos + alt_allele.length();
            if (alt_allele == ".")
                end_pos = record.pos;

            append_kmer_covgs_in_range(kg, alt_kmer_path, alt_path, record.pos, end_pos,
                                       alt_fwd_covgs, alt_rev_covgs, sample_id);

            record.samples[sample_index]["MEAN_FWD_COVG"].push_back(mean(alt_fwd_covgs));
            record.samples[sample_index]["MEAN_REV_COVG"].push_back(mean(alt_rev_covgs));
            record.samples[sample_index]["MED_FWD_COVG"].push_back(median(alt_fwd_covgs));
            record.samples[sample_index]["MED_REV_COVG"].push_back(median(alt_rev_covgs));
            record.samples[sample_index]["SUM_FWD_COVG"].push_back(sum(alt_fwd_covgs));
            record.samples[sample_index]["SUM_REV_COVG"].push_back(sum(alt_rev_covgs));
            alt_fwd_covgs.clear();
            alt_rev_covgs.clear();
        }

        record.add_formats({"MEAN_FWD_COVG", "MEAN_REV_COVG",
                            "MED_FWD_COVG", "MED_REV_COVG",
                            "SUM_FWD_COVG", "SUM_REV_COVG"});

        ref_fwd_covgs.clear();
        ref_rev_covgs.clear();
        alt_fwd_covgs.clear();
        alt_rev_covgs.clear();
    }

}

void LocalPRG::add_consensus_path_to_fastaq(Fastaq &output_fq,
                                            PanNodePtr pnode,
                                            std::vector<KmerNodePtr> &kmp,
                                            std::vector<LocalNodePtr> &lmp,
                                            const uint32_t w,
                                            const bool bin,
                                            const uint32_t global_covg,
                                            const uint32_t &sample_id) {
    kmp.clear();
    lmp.clear();
    if (pnode->reads.empty()) {
        BOOST_LOG_TRIVIAL(warning) << "Node " << pnode->get_name() << " has no reads";
        return;
    }
    kmp.reserve(800);

    BOOST_LOG_TRIVIAL(info) << "Find maxpath for " << pnode->get_name();
    float ppath;
    if (bin)
        ppath = pnode->kmer_prg.find_max_path(kmp, sample_id);
    else
        ppath = pnode->kmer_prg.find_nb_max_path(kmp, sample_id);

    lmp.reserve(100);
    lmp = localnode_path_from_kmernode_path(kmp, w);

    std::vector<uint32_t> covgs = get_covgs_along_localnode_path(pnode, lmp, kmp, sample_id);
    auto mode_covg = mode(covgs);
    auto mean_covg = mean(covgs);
    BOOST_LOG_TRIVIAL(info) << "Found global coverage " << global_covg << " and path mode " << mode_covg << " and mean "
                            << mean_covg;
    if (global_covg > 20 and
        (6 * mean(covgs) < global_covg or mean(covgs) > global_covg or (mode(covgs) < 3 and mean(covgs) < 3))) {
        BOOST_LOG_TRIVIAL(info) << "Skip LocalPRG " << name << " as mode and mean along max likelihood path too low";
        kmp.clear();
        return;
    }

    std::string fq_name = pnode->get_name();
    std::string header = " log P(data|sequence)=" + std::to_string(ppath);
    std::string seq = string_along_path(lmp);
    output_fq.add_entry(fq_name, seq, covgs, global_covg, header);
}

std::vector<LocalNodePtr> LocalPRG::get_valid_vcf_reference(const std::string &vcf_reference_sequence) const {
    std::vector<LocalNodePtr> reference_path = {};
    if (vcf_reference_sequence.length() < 30) {
        BOOST_LOG_TRIVIAL(warning) << "Input vcf_ref path was too short to be the ref for PRG " << name;
        return reference_path;
    }

    reference_path = this->prg.nodes_along_string(vcf_reference_sequence);
    if (reference_path.empty()) {
        reference_path = this->prg.nodes_along_string(rev_complement(vcf_reference_sequence));
    }

    if (reference_path.empty())
        return reference_path;

    bool not_starting_at_prg_start = reference_path.front()->pos.start != 0;

    LocalNode last_prg_node = *(*(prg.nodes.rbegin())).second;
    auto final_prg_coordinate = last_prg_node.pos.get_end();
    bool not_ending_at_prg_end = reference_path.back()->pos.get_end() != final_prg_coordinate;

    if (not_starting_at_prg_start or not_ending_at_prg_end) {
        BOOST_LOG_TRIVIAL(warning) << "Input vcf_ref path did not start/end at the beginning/end of PRG "
                                   << name;
        reference_path.clear();
    }
    return reference_path;
}

void LocalPRG::add_variants_to_vcf(VCF &master_vcf,
                                   PanNodePtr pnode,
                                   const std::string &vcf_ref,
                                   const std::vector<KmerNodePtr> &kmp,
                                   const std::vector<LocalNodePtr> &lmp,
                                   const uint32_t &sample_id,
                                   const std::string &sample_name) {
    auto reference_path = get_valid_vcf_reference(vcf_ref);
    if (reference_path.empty()) {
        BOOST_LOG_TRIVIAL(warning) << "Could not find reference sequence for " << name
                                   << "in the PRG so using the consensus path";
        reference_path = lmp;
    }

    VCF vcf;
    build_vcf(vcf, reference_path);
    //std::cout << "add sample gts" << std::endl;
    add_sample_gt_to_vcf(vcf, reference_path, lmp, sample_name);
    //std::cout << "add sample covgs" << std::endl;
    add_sample_covgs_to_vcf(vcf, pnode->kmer_prg, reference_path, sample_name, sample_id);
    //vcf.save("temp.vcf" , true, true, true, true, true, true, true);
    //std::cout << "sort records" << std::endl;
    //vcf.sort_records();
    //std::cout << "append to master vcf" << std::endl;
    vcf.merge_multi_allelic();
    vcf.correct_dot_alleles(string_along_path(reference_path), name);
    master_vcf.append_vcf(vcf);
}


bool operator!=(const std::vector<KmerNodePtr> &lhs, const std::vector<KmerNodePtr> &rhs) {
    if (lhs.size() != rhs.size()) {
        return true;
    }
    for (uint32_t i = 0; i != lhs.size(); ++i) {
        if (lhs[i] != rhs[i]) {
            return true;
        }
    }
    return false;
}

std::ostream &operator<<(std::ostream &out, LocalPRG const &data) {
    out << data.name;
    return out;
}

