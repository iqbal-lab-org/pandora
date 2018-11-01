#include <iostream>
#include <map>
#include <unordered_set>
#include <set>
#include <utility>
#include <vector>
#include <fstream>
#include <cassert>
#include "utils.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"
#include "de_bruijn/graph.h"
#include "minihit.h"


uint_least32_t node_plus_orientation_to_num(const uint_least32_t node_id, const bool orientation) {
    assert(node_id < UINT_LEAST32_MAX / 2);
    uint_least32_t r = 2 * node_id;
    if (orientation) {
        r += 1;
    }
    return r;
}

void num_to_node_plus_orientation(uint_least32_t &node_id, bool &orientation, const uint_least32_t num) {
    if (num % 2 == 1) {
        orientation = true;
        node_id = (num - 1) / 2;
    } else {
        orientation = false;
        node_id = num / 2;
    }
}

uint_least32_t rc_num(const uint_least32_t &num) {
    return num + 1 * (num % 2 == 0) - 1 * (num % 2 == 1);
}

void hashed_node_ids_to_ids_and_orientations(const std::deque<uint_least32_t> &hashed_node_ids,
                                             std::vector<uint_least32_t> &node_ids,
                                             std::vector<bool> &node_orients) {
    node_ids.clear();
    node_orients.clear();

    uint_least32_t node_id;
    bool orientation;
    for (const auto &i : hashed_node_ids) {
        num_to_node_plus_orientation(node_id, orientation, i);
        node_ids.push_back(node_id);
        node_orients.push_back(orientation);
    }
}

bool overlap_forwards(const std::deque<uint_least32_t> &node1, const std::deque<uint_least32_t> &node2) {
    // second deque should extend first by 1
    assert(node1.size() >= node2.size());
    uint32_t i = node1.size() - node2.size() + 1;
    uint32_t j = 0;
    while (i < node1.size() and j < node2.size()) {
        if (node1[i] != node2[j]) {
            return false;
        }
        i++;
        j++;
    }
    return true;
}

bool overlap_backwards(const std::deque<uint_least32_t> &node1, const std::deque<uint_least32_t> &node2) {
    for (uint32_t i = 1; i < std::min(node1.size() + 1, node2.size()); ++i) {
        if (node2[i] != node1[i - 1]) {
            return false;
        }
    }
    return true;
}

std::deque<uint_least32_t> rc_hashed_node_ids(const std::deque<uint_least32_t> &hashed_node_ids) {
    /*for (const auto &n : hashed_node_ids)
    {
        cout << n << " ";
    }*/
    std::deque<uint_least32_t> d;
    for (const auto &i : hashed_node_ids) {
        d.push_front(rc_num(i));
    }
    /*cout << " is now ";
    for (const auto &n : d)
    {
        cout << n << " ";
    }
    cout << endl;*/
    return d;
}

std::deque<uint_least32_t> extend_hashed_pg_node_ids_backwards(const debruijn::Graph &dbg,
                                                               const std::deque<uint32_t> &dbg_node_ids) {
    std::deque<uint_least32_t> hashed_pg_node_ids = dbg.nodes.at(dbg_node_ids.at(0))->hashed_node_ids;
    std::deque<uint_least32_t> rev_node;

    for (uint32_t i = 1; i < dbg_node_ids.size(); ++i) {
        rev_node = rc_hashed_node_ids(dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids);
        if (overlap_backwards(hashed_pg_node_ids, dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids)) {
            hashed_pg_node_ids.push_front(dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids[0]);
        } else if (overlap_backwards(hashed_pg_node_ids, rev_node)) {
            hashed_pg_node_ids.push_front(rc_num(dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids.back()));
        } else {
            hashed_pg_node_ids.clear();
            break;
        }
    }
    return hashed_pg_node_ids;
}

std::deque<uint_least32_t> extend_hashed_pg_node_ids_forwards(const debruijn::Graph &dbg,
                                                              const std::deque<uint32_t> &dbg_node_ids) {
    std::deque<uint_least32_t> hashed_pg_node_ids = dbg.nodes.at(dbg_node_ids.at(0))->hashed_node_ids;
    std::deque<uint_least32_t> rev_node;

    for (uint32_t i = 1; i < dbg_node_ids.size(); ++i) {
        rev_node = rc_hashed_node_ids(dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids);
        if (overlap_forwards(hashed_pg_node_ids, dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids)) {
            hashed_pg_node_ids.push_back(dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids.back());
        } else if (overlap_forwards(hashed_pg_node_ids, rev_node)) {
            hashed_pg_node_ids.push_back(rc_num(dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids[0]));
        } else {
            hashed_pg_node_ids.clear();
            break;
        }
    }
    return hashed_pg_node_ids;
}

void dbg_node_ids_to_ids_and_orientations(const debruijn::Graph &dbg,
                                          const std::deque<uint32_t> &dbg_node_ids,
                                          std::vector<uint_least32_t> &node_ids,
                                          std::vector<bool> &node_orients) {
    node_ids.clear();
    node_orients.clear();

    if (dbg_node_ids.empty()) {
        return;
    }

    std::deque<uint_least32_t> hashed_pg_node_ids = extend_hashed_pg_node_ids_backwards(dbg, dbg_node_ids);
    if (hashed_pg_node_ids.empty())
        hashed_pg_node_ids = extend_hashed_pg_node_ids_forwards(dbg, dbg_node_ids);
    if (hashed_pg_node_ids.empty()) {
        for (const auto &n : dbg_node_ids) {
            std::cout << "(";
            for (const auto &m : dbg.nodes.at(n)->hashed_node_ids) {
                std::cout << m << " ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
    }
    assert(!hashed_pg_node_ids.empty());
    hashed_node_ids_to_ids_and_orientations(hashed_pg_node_ids, node_ids, node_orients);
}

void construct_debruijn_graph(const pangenome::Graph *pg, debruijn::Graph &dbg) {
    dbg.nodes.clear();
    dbg.node_hash.clear();
    //cout << "Removed old nodes" << endl;

    debruijn::OrientedNodePtr prev, current;
    std::deque<uint_least32_t> hashed_ids;

    for (const auto &r : pg->reads) {
        if (r.second->nodes.size() < dbg.size) {
            // can't add anything for this read
            continue;
        }

        prev = std::make_pair(nullptr, false);
        current = std::make_pair(nullptr, false);
        hashed_ids.clear();

        for (uint32_t i = 0; i < r.second->nodes.size(); ++i) {
            hashed_ids.push_back(node_plus_orientation_to_num(r.second->nodes[i]->node_id,
                                                              r.second->node_orientations[i]));

            if (hashed_ids.size() == dbg.size) {
                current = dbg.add_node(hashed_ids, r.first);
                if (prev.first != nullptr and current.first != nullptr) {
                    dbg.add_edge(prev, current);
                }
                prev = current;
                hashed_ids.pop_front();
            }
        }
    }
}

void remove_leaves(pangenome::Graph *pg, debruijn::Graph &dbg, uint_least32_t covg_thresh) {
    std::cout << now() << "Remove leaves of debruijn graph from pangraph" << std::endl;
    std::cout << "Start with " << pg->nodes.size() << " pg.nodes, " << pg->reads.size() << " pg.reads, and "
              << dbg.nodes.size() << " dbg.nodes" << std::endl;
    bool leaves_exist = true;
    std::unordered_set<uint32_t> leaves;
    std::vector<uint_least32_t> node_ids;
    std::vector<bool> node_orients;
    std::pair<uint32_t, uint32_t> pos;
    pangenome::NodePtr node;

    while (leaves_exist) {

        leaves = dbg.get_leaves(covg_thresh);
        std::cout << "there are " << leaves.size() << " leaves" << std::endl;

        if (leaves.empty()) {
            leaves_exist = false;
        }
        std::cout << "leaves exist is " << leaves_exist << std::endl;

        for (const auto &i : leaves) {
            std::cout << std::endl << "looking at leaf " << i << ": ";
            for (const auto &j : dbg.nodes[i]->hashed_node_ids) {
                std::cout << j << " ";
            }
            std::cout << std::endl;

            // look up the node ids and orientations associated with this node
            hashed_node_ids_to_ids_and_orientations(dbg.nodes[i]->hashed_node_ids, node_ids, node_orients);
            std::cout << "looked up node ids" << std::endl;

            // remove the last node from corresponding reads
            assert(not dbg.nodes[i]->read_ids.empty());
            for (const auto &r : dbg.nodes[i]->read_ids) {
                std::cout << "remove from read " << r << ": ";
                for (const auto &n : pg->reads[r]->nodes) {
                    std::cout << n->node_id << " ";
                }
                std::cout << std::endl;
                if (pg->reads[r]->nodes.size() == dbg.size) {
                    std::cout << "remove read";
                    pg->remove_read(r);
                    std::cout << " done" << std::endl;
                } else {
                    std::cout << "remove from read ";
                    pos = pg->reads[r]->find_position(node_ids, node_orients);
                    std::cout << " pos " << pos.first << " " << pos.second;
                    assert(pos.first == 0 or pos.first + node_ids.size() == pg->reads[r]->nodes.size());
                    if (pos.first == 0) {
                        node = pg->reads[r]->nodes[0];
                        pg->reads[r]->remove_node(pg->reads[r]->nodes.begin());
                        node->remove_read(pg->reads[r]);
                    } else if (pos.first + node_ids.size() == pg->reads[r]->nodes.size()) {
                        node = pg->reads[r]->nodes.back();
                        pg->reads[r]->remove_node(--pg->reads[r]->nodes.end());
                        node->remove_read(pg->reads[r]);
                    }
                    std::cout << "read is now " << r << ": ";
                    for (const auto &n : pg->reads[r]->nodes) {
                        std::cout << n->node_id << " ";
                    }
                    std::cout << "done" << std::endl;
                }
            }
            if (node and node->covg == 0) {
                pg->remove_node(node);
            }

            // remove dbg node
            dbg.remove_node(i);

            //cout << "pg is now: " << endl << pg << endl;
        }
    }
    std::cout << "There are now " << pg->nodes.size() << " pg.nodes, " << pg->reads.size() << " pg.reads, and "
              << dbg.nodes.size() << " dbg.nodes" << std::endl;
}

void find_reads_along_tig(const debruijn::Graph &dbg,
                          std::deque<uint32_t> &dbg_node_ids,
                          const pangenome::Graph *pg,
                          std::vector<uint_least32_t> &pg_node_ids,
                          std::vector<bool> &pg_node_orients,
                          std::unordered_set<pangenome::ReadPtr> &reads_along_tig,
                          bool &all_reads_along_tig) {
    // collect the reads covering that tig
    for (const auto &n : dbg_node_ids) {
        for (const auto &r : dbg.nodes.at(n)->read_ids) {
            reads_along_tig.insert(pg->reads.at(r));
        }
    }
    std::cout << "candidate reads ";
    for (const auto &r : reads_along_tig) {
        std::cout << r->id << " ";
    }
    std::cout << std::endl;

    // filter out some which don't really overlap the unitig, keeping those
    // which overlap at least consecutive 2 dbg nodes or only one node
    all_reads_along_tig = true;
    std::cout << "kept reads along tig: ";
    for (auto r = reads_along_tig.begin(); r != reads_along_tig.end();) {
        if ((*r)->nodes.size() > dbg.size and
            (*r)->find_position(pg_node_ids, pg_node_orients, dbg.size + 1).first ==
            std::numeric_limits<uint32_t>::max()) {
            r = reads_along_tig.erase(r);
            all_reads_along_tig = false;
        } else {
            std::cout << (*r)->id << " ";
            ++r;
        }
    }
    //cout << endl;
}

void remove_middle_nodes_of_tig_from_read(pangenome::Graph *pg,
                                          debruijn::Graph &dbg,
                                          pangenome::ReadPtr r,
                                          const std::vector<uint_least32_t> &node_ids,
                                          const std::vector<bool> &node_orients) {
    // suppose tig is subset of read
    // want to remove nodes from position pos (where tig starts in read) + dbg.size()
    // to position start + tig.size()-dbg.size()
    // if pos = 0 (ie read does not include start of tig)

    std::pair<uint32_t, uint32_t> pos = r->find_position(node_ids, node_orients);
    std::cout << "found pos " << pos.first << " " << pos.second << std::endl;
    auto start_shift = pos.first;
    if (pos.first > 0 or pos.second < r->nodes.size() - 1 or node_ids.size() == r->nodes.size())
        start_shift += std::max((int) 0, (int) pos.second - (int) node_ids.size()) + dbg.size;
    else {
        std::vector<uint_least32_t> sub_node_ids(node_ids.begin() + dbg.size, node_ids.end());
        std::vector<bool> sub_node_orients(node_orients.begin() + dbg.size, node_orients.end());
        std::pair<uint32_t, uint32_t> sub_pos = r->find_position(sub_node_ids, sub_node_orients);
        if (sub_pos.first > 0)
            start_shift = sub_pos.first;
    }


    auto end_shift = pos.second;
    if (pos.first > 0 or pos.second < r->nodes.size() - 1 or node_ids.size() == r->nodes.size()) {
        end_shift -= dbg.size - 1;
    } else {
        std::vector<uint_least32_t> sub_node_ids(node_ids.begin(), node_ids.end() - dbg.size);
        std::vector<bool> sub_node_orients(node_orients.begin(), node_orients.end() - dbg.size);
        std::pair<uint32_t, uint32_t> sub_pos = r->find_position(sub_node_ids, sub_node_orients);
        if (sub_pos.second < pos.second)
            end_shift = sub_pos.second + 1;

    }

    std::cout << "remove nodes " << start_shift << " to " << end_shift << std::endl;
    auto it = r->nodes.begin() + start_shift;
    for (auto shift = start_shift; shift < end_shift; ++shift) {
        if (it == r->nodes.end()) {
            break;
        }
        it = pg->remove_node_from_read(it, r);
    }
}

// Remove the internal nodes of low coverage unitigs e.g.
// suppose when dbg kmer size is 3, we have a low covg tig
// 012 -> 126 -> 263 -> 634 -> 345
// then we would remove the 3 internal kmers from the dbg
// and node 6 from the pg->
// If the tig is smaller than k+2 long, currently does nothing
void filter_unitigs(pangenome::Graph *pg, debruijn::Graph &dbg, const uint_least32_t &threshold) {
    std::cout << now() << "Filter unitigs using threshold " << threshold << std::endl;
    std::vector<uint_least32_t> node_ids;
    std::vector<bool> node_orients;
    std::unordered_set<pangenome::ReadPtr> reads_along_tig;
    bool all_reads_tig;

    std::set<std::deque<uint32_t>> unitigs = dbg.get_unitigs();
    std::cout << "have " << unitigs.size() << " tigs" << std::endl;
    for (auto d : unitigs) {
        // look up the node ids and orientations associated with this node
        dbg_node_ids_to_ids_and_orientations(dbg, d, node_ids, node_orients);
        std::cout << "tig: ";
        for (const auto &n : node_ids) {
            std::cout << n << " ";
        }
        std::cout << std::endl;

        // collect the reads covering that tig
        find_reads_along_tig(dbg, d, pg, node_ids, node_orients, reads_along_tig, all_reads_tig);

        // now if the number of reads covering tig falls below threshold, remove the
        // middle nodes of this tig from the reads
        if (reads_along_tig.size() <= threshold) {
            std::cout << "not enough reads, so remove the tig from the reads" << std::endl;
            for (const auto &r : reads_along_tig) {
                std::cout << "read " << r->id << " was ";
                for (const auto &n : r->nodes) {
                    std::cout << n->node_id << " ";
                }
                std::cout << std::endl;
                remove_middle_nodes_of_tig_from_read(pg, dbg, r, node_ids, node_orients);
                std::cout << "read " << r->id << " is now ";
                for (const auto &n : r->nodes) {
                    std::cout << n->node_id << " ";
                }
                std::cout << std::endl;
            }
            // also remove read_ids from each of the corresponding nodes of dbg
            //cout << "now remove read from dbg nodes" << endl;
            for (uint32_t i = 1; i < d.size() - 1; ++i) {
                for (const auto &r :reads_along_tig) {
                    dbg.remove_read_from_node(r->id, d[i]);
                }
            }
        } else {
            std::cout << "tig had enough reads" << std::endl;
        }
        reads_along_tig.clear();
    }
}

void detangle_pangraph_with_debruijn_graph(pangenome::Graph *pg, debruijn::Graph &dbg) {
    std::cout << now() << "Detangle pangraph with debruijn graph" << std::endl;
    std::vector<uint_least32_t> node_ids;
    std::vector<bool> node_orients;
    bool all_reads_tig;
    std::unordered_set<pangenome::ReadPtr> reads_along_tig;

    std::set<std::deque<uint32_t>> unitigs = dbg.get_unitigs();
    for (auto d : unitigs) {
        // look up the node ids and orientations associated with this node
        dbg_node_ids_to_ids_and_orientations(dbg, d, node_ids, node_orients);
        /*cout << "tig: ";
        for (const auto &n : node_ids)
        {
            cout << n << " ";
        }
        cout << endl;*/

        // collect the reads covering that tig
        find_reads_along_tig(dbg, d, pg, node_ids, node_orients, reads_along_tig, all_reads_tig);

        // for each node on tig, for each read covering that node,
        // if we find a read which doesn't lie along whole tig,
        // split that node by reads and create a new node on the tig
        if (!all_reads_tig and !reads_along_tig.empty()) {
            //cout << "not all reads contain tig" << endl;
            for (uint32_t i = 0; i < node_ids.size(); ++i) {
                //cout << "for node " << pg->nodes[node_ids[i]]->node_id << endl;
                for (const auto &r : pg->nodes[node_ids[i]]->reads) {
                    //cout << "for read " << r->id << endl;
                    if (reads_along_tig.find(r) == reads_along_tig.end()) {
                        //cout << "split node" << endl;
                        pg->split_node_by_reads(reads_along_tig, node_ids, node_orients, node_ids[i]);
                        //cout << "done" << endl;
                        break;
                    }
                }
                //cout << pg << endl;
            }
        }
    }
}

void clean_pangraph_with_debruijn_graph(pangenome::Graph *pg, const uint_least32_t size, const uint_least32_t threshold,
                                        const bool illumina) {
    std::cout << now() << "Construct de Bruijn Graph from PanGraph with size " << (uint32_t) size << std::endl;
    debruijn::Graph dbg(size);
    construct_debruijn_graph(pg, dbg);

    if (not illumina)
        remove_leaves(pg, dbg, threshold);
    filter_unitigs(pg, dbg, threshold);
    std::cout << "Finished filtering tigs" << std::endl;

    // update dbg now that have removed leaves and some inner nodes
    std::cout << "Reconstruct dbg" << std::endl;
    construct_debruijn_graph(pg, dbg);
    std::cout << "Now detangle" << std::endl;

    detangle_pangraph_with_debruijn_graph(pg, dbg);
}


enum NodeDirection {
    forward,
    reverse
};


NodeDirection get_pangraph_node_direction(const debruijn::Node &debruijn_node) {
    bool forward_node = debruijn_node.hashed_node_ids[0] % 2 != 0;
    if (forward_node)
        return NodeDirection::forward;
    else
        return NodeDirection::reverse;
}


uint32_t get_pangraph_node_id(const debruijn::Node &debruijn_node) {
    auto direction = get_pangraph_node_direction(debruijn_node);
    if (direction == NodeDirection::forward)
        return (uint32_t) (debruijn_node.hashed_node_ids[0] - 1) / 2;
    else
        return (uint32_t) debruijn_node.hashed_node_ids[0] / 2;
}


namespace gfa {
    namespace dump {
        void header(std::ofstream &gfa_fhandle) {
            gfa_fhandle << "H\tVN:Z:1.0" << std::endl;
        }


        void nodes(const pangenome::Graph *const pangraph, std::ofstream &gfa_fhandle) {
            for (const auto &node : pangraph->nodes)
                gfa_fhandle << "S\t" << node.second->get_name() << "\tN\tFC:i:" << node.second->covg << std::endl;
        }

        void edge(const pangenome::Node &first_node,
                  const NodeDirection &first_node_direction,
                  const pangenome::Node &second_node,
                  const NodeDirection &second_node_direction,
                  std::ofstream &gfa_fhandle) {
            gfa_fhandle << "L\t" << first_node.get_name() << "\t";
            if (first_node_direction == NodeDirection::forward)
                gfa_fhandle << "-";
            else
                gfa_fhandle << "+";

            gfa_fhandle << "\t" << second_node.get_name() << "\t";
            if (second_node_direction == NodeDirection::forward)
                gfa_fhandle << "-";
            else
                gfa_fhandle << "+";

            gfa_fhandle << "\t0M" << std::endl;
        }
    }
}


pangenome::Node convert_node_debruijn_pangraph(const debruijn::Node &debruijn_node,
                                               const pangenome::Graph *const pangraph) {
    auto node_id = get_pangraph_node_id(debruijn_node);
    assert(pangraph->nodes.find(node_id) != pangraph->nodes.end());

    auto node_ptr = pangraph->nodes.at(node_id);
    auto node = *node_ptr;
    return node;
}


void write_pangraph_gfa(const std::string &filepath, const pangenome::Graph *const pangraph) {
    std::ofstream gfa_fhandle;
    gfa_fhandle.open(filepath);
    gfa::dump::header(gfa_fhandle);
    gfa::dump::nodes(pangraph, gfa_fhandle);

    debruijn::Graph debruijn_graph(1);
    construct_debruijn_graph(pangraph, debruijn_graph);

    for (const auto &debruijn_node_id_entry : debruijn_graph.nodes) {
        auto first_debruijn_node = *debruijn_node_id_entry.second;

        auto first_node = convert_node_debruijn_pangraph(first_debruijn_node, pangraph);
        auto first_node_direction = get_pangraph_node_direction(first_debruijn_node);

        for (const auto &second_debruijn_node_id: first_debruijn_node.out_nodes) {
            assert(debruijn_graph.nodes.find(second_debruijn_node_id) != debruijn_graph.nodes.end());
            auto &second_debruijn_node = *debruijn_graph.nodes[second_debruijn_node_id];

            auto second_node = convert_node_debruijn_pangraph(second_debruijn_node, pangraph);
            auto second_node_direction = get_pangraph_node_direction(second_debruijn_node);

            gfa::dump::edge(first_node, first_node_direction,
                            second_node, second_node_direction,
                            gfa_fhandle);

            auto erased = second_debruijn_node.out_nodes.erase(first_debruijn_node.id);
            if (erased)
                continue;
            second_debruijn_node.in_nodes.erase(first_debruijn_node.id);
        }
    }
}
