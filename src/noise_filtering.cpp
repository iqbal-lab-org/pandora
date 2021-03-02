#include <iostream>
#include <map>
#include <unordered_set>
#include <set>
#include <utility>
#include <vector>
#include "utils.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "de_bruijn/graph.h"
#include "minihit.h"

uint_least32_t node_plus_orientation_to_num(
    const uint_least32_t node_id, const bool orientation)
{
    const bool node_id_is_consistent = node_id < UINT_LEAST32_MAX / 2;
    if(!node_id_is_consistent) {
        fatal_error("Error on converting node id and orientation to id only: "
                    "node_id (", node_id, ") should be < than ", UINT_LEAST32_MAX / 2);
    }
    uint_least32_t r = 2 * node_id;
    if (orientation) {
        r += 1;
    }
    return r;
}

void num_to_node_plus_orientation(
    uint_least32_t& node_id, bool& orientation, const uint_least32_t num)
{
    if (num % 2 == 1) {
        orientation = true;
        node_id = (num - 1) / 2;
    } else {
        orientation = false;
        node_id = num / 2;
    }
}

uint_least32_t rc_num(const uint_least32_t& num)
{
    return num + 1 * (num % 2 == 0) - 1 * (num % 2 == 1);
}

void hashed_node_ids_to_ids_and_orientations(
    const std::deque<uint_least32_t>& hashed_node_ids,
    std::vector<uint_least32_t>& node_ids, std::vector<bool>& node_orients)
{
    node_ids.clear();
    node_orients.clear();

    uint_least32_t node_id;
    bool orientation;
    for (const auto& i : hashed_node_ids) {
        num_to_node_plus_orientation(node_id, orientation, i);
        node_ids.push_back(node_id);
        node_orients.push_back(orientation);
    }
}

bool overlap_forwards(
    const std::deque<uint_least32_t>& node1, const std::deque<uint_least32_t>& node2)
{
    // second deque should extend first by 1
    const bool first_node_is_larger_or_same_size = node1.size() >= node2.size();
    if(!first_node_is_larger_or_same_size) {
        fatal_error("Error on checking for overlaps in noise filtering: first node must be larger or have the same size as the second");
    }

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

bool overlap_backwards(
    const std::deque<uint_least32_t>& node1, const std::deque<uint_least32_t>& node2)
{
    for (uint32_t i = 1; i < std::min(node1.size() + 1, node2.size()); ++i) {
        if (node2[i] != node1[i - 1]) {
            return false;
        }
    }
    return true;
}

std::deque<uint_least32_t> rc_hashed_node_ids(
    const std::deque<uint_least32_t>& hashed_node_ids)
{
    std::deque<uint_least32_t> d;
    for (const auto& i : hashed_node_ids) {
        d.push_front(rc_num(i));
    }
    return d;
}

std::deque<uint_least32_t> extend_hashed_pg_node_ids_backwards(
    const debruijn::Graph& dbg, const std::deque<uint32_t>& dbg_node_ids)
{
    std::deque<uint_least32_t> hashed_pg_node_ids
        = dbg.nodes.at(dbg_node_ids.at(0))->hashed_node_ids;
    std::deque<uint_least32_t> rev_node;

    for (uint32_t i = 1; i < dbg_node_ids.size(); ++i) {
        rev_node
            = rc_hashed_node_ids(dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids);
        if (overlap_backwards(hashed_pg_node_ids,
                dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids)) {
            hashed_pg_node_ids.push_front(
                dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids[0]);
        } else if (overlap_backwards(hashed_pg_node_ids, rev_node)) {
            hashed_pg_node_ids.push_front(
                rc_num(dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids.back()));
        } else {
            hashed_pg_node_ids.clear();
            break;
        }
    }
    return hashed_pg_node_ids;
}

std::deque<uint_least32_t> extend_hashed_pg_node_ids_forwards(
    const debruijn::Graph& dbg, const std::deque<uint32_t>& dbg_node_ids)
{
    std::deque<uint_least32_t> hashed_pg_node_ids
        = dbg.nodes.at(dbg_node_ids.at(0))->hashed_node_ids;
    std::deque<uint_least32_t> rev_node;

    for (uint32_t i = 1; i < dbg_node_ids.size(); ++i) {
        rev_node
            = rc_hashed_node_ids(dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids);
        if (overlap_forwards(hashed_pg_node_ids,
                dbg.nodes.at(dbg_node_ids.at(i))->hashed_node_ids)) {
            hashed_pg_node_ids.push_back(
                dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids.back());
        } else if (overlap_forwards(hashed_pg_node_ids, rev_node)) {
            hashed_pg_node_ids.push_back(
                rc_num(dbg.nodes.at(dbg_node_ids[i])->hashed_node_ids[0]));
        } else {
            hashed_pg_node_ids.clear();
            break;
        }
    }
    return hashed_pg_node_ids;
}

void dbg_node_ids_to_ids_and_orientations(const debruijn::Graph& dbg,
    const std::deque<uint32_t>& dbg_node_ids, std::vector<uint_least32_t>& node_ids,
    std::vector<bool>& node_orients)
{
    node_ids.clear();
    node_orients.clear();

    if (dbg_node_ids.empty()) {
        return;
    }

    std::deque<uint_least32_t> hashed_pg_node_ids
        = extend_hashed_pg_node_ids_backwards(dbg, dbg_node_ids);
    if (hashed_pg_node_ids.empty()) {
        hashed_pg_node_ids = extend_hashed_pg_node_ids_forwards(dbg, dbg_node_ids);
    }

    // TODO: give a better name to this bool once we understand what it does
    const bool hashed_pg_node_ids_is_empty = hashed_pg_node_ids.empty();
    if(hashed_pg_node_ids_is_empty) {
        // TODO: improve this message
        fatal_error("Error when noise filtering: hashed_pg_node_ids is empty");
    }
    hashed_node_ids_to_ids_and_orientations(hashed_pg_node_ids, node_ids, node_orients);
}

void construct_debruijn_graph(
    std::shared_ptr<pangenome::Graph> pangraph, debruijn::Graph& dbg)
{
    dbg.nodes.clear();
    dbg.node_hash.clear();

    debruijn::OrientedNodePtr prev, current;
    std::deque<uint_least32_t> hashed_ids;

    for (const auto& r : pangraph->reads) {
        if (r.second->get_nodes().size() < dbg.size) {
            // can't add anything for this read
            continue;
        }

        prev = std::make_pair(nullptr, false);
        current = std::make_pair(nullptr, false);
        hashed_ids.clear();

        for (uint32_t i = 0; i < r.second->get_nodes().size(); ++i) {
            hashed_ids.push_back(
                node_plus_orientation_to_num(r.second->get_nodes()[i].lock()->node_id,
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

void remove_leaves(std::shared_ptr<pangenome::Graph> pangraph, debruijn::Graph& dbg,
    uint_least32_t covg_thresh)
{
    BOOST_LOG_TRIVIAL(debug) << "Remove leaves of debruijn graph from pangraph";
    BOOST_LOG_TRIVIAL(debug) << "Start with " << pangraph->nodes.size() << " pg.nodes, "
                             << pangraph->reads.size() << " pg.reads, and "
                             << dbg.nodes.size() << " dbg.nodes";
    bool leaves_exist = true;
    std::unordered_set<uint32_t> leaves;
    std::vector<uint_least32_t> node_ids;
    std::vector<bool> node_orients;
    std::pair<uint32_t, uint32_t> pos;
    pangenome::WeakNodePtr node;

    while (leaves_exist) {

        leaves = dbg.get_leaves(covg_thresh);
        BOOST_LOG_TRIVIAL(trace) << "there are " << leaves.size() << " leaves";

        if (leaves.empty()) {
            leaves_exist = false;
        }
        BOOST_LOG_TRIVIAL(trace) << "leaves exist is " << leaves_exist;

        for (const auto& i : leaves) {
            // look up the node ids and orientations associated with this node
            hashed_node_ids_to_ids_and_orientations(
                dbg.nodes[i]->hashed_node_ids, node_ids, node_orients);

            const bool dbg_node_has_no_reads = dbg.nodes[i]->read_ids.empty();
            if (dbg_node_has_no_reads) {
                fatal_error("Error when removing leaves from DBG: node has no leaves");
            }

            // remove the last node from corresponding reads
            for (const auto& r : dbg.nodes[i]->read_ids) {
                if (pangraph->reads[r]->get_nodes().size() == dbg.size) {
                    pangraph->remove_read(r);
                } else {
                    pos = pangraph->reads[r]->find_position(node_ids, node_orients);

                    const bool pos_of_nodes_in_read_is_valid = (pos.first == 0) or
                        (pos.first + node_ids.size() == pangraph->reads[r]->get_nodes().size());
                    if (!pos_of_nodes_in_read_is_valid) {
                        fatal_error("Error when removing leaves from DBG: position of "
                                    "DBG nodes in reads are not valid");
                    }

                    if (pos.first == 0) {
                        node = pangraph->reads[r]->get_nodes()[0];
                        pangraph->reads[r]->remove_node_with_iterator(
                            pangraph->reads[r]->get_nodes().begin());
                        node.lock()->remove_read(pangraph->reads[r]);
                    } else if (pos.first + node_ids.size()
                        == pangraph->reads[r]->get_nodes().size()) {
                        node = pangraph->reads[r]->get_nodes().back();
                        pangraph->reads[r]->remove_all_nodes_with_this_id(
                            (--pangraph->reads[r]->get_nodes().end())->lock()->node_id);
                        node.lock()->remove_read(pangraph->reads[r]);
                    }
                }
            }
            auto node_shared_ptr_from_weak_ptr = node.lock();
            if (node_shared_ptr_from_weak_ptr
                and node_shared_ptr_from_weak_ptr->covg == 0) {
                pangraph->remove_node(node_shared_ptr_from_weak_ptr);
            }

            // remove dbg node
            dbg.remove_node(i);
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "There are now " << pangraph->nodes.size()
                             << " pg.nodes, " << pangraph->reads.size()
                             << " pg.reads, and " << dbg.nodes.size() << " dbg.nodes";
}

void find_reads_along_tig(const debruijn::Graph& dbg,
    std::deque<uint32_t>& dbg_node_ids, std::shared_ptr<pangenome::Graph> pangraph,
    std::vector<uint_least32_t>& pg_node_ids, std::vector<bool>& pg_node_orients,
    std::unordered_set<pangenome::ReadPtr>& reads_along_tig, bool& all_reads_along_tig)
{
    // collect the reads covering that tig
    for (const auto& n : dbg_node_ids) {
        for (const auto& r : dbg.nodes.at(n)->read_ids) {
            reads_along_tig.insert(pangraph->reads.at(r));
        }
    }
    // filter out some which don't really overlap the unitig, keeping those
    // which overlap at least consecutive 2 dbg nodes or only one node
    all_reads_along_tig = true;
    for (auto r = reads_along_tig.begin(); r != reads_along_tig.end();) {
        if ((*r)->get_nodes().size() > dbg.size
            and (*r)->find_position(pg_node_ids, pg_node_orients, dbg.size + 1).first
                == std::numeric_limits<uint32_t>::max()) {
            r = reads_along_tig.erase(r);
            all_reads_along_tig = false;
        } else {
            ++r;
        }
    }
}

void remove_middle_nodes_of_tig_from_read(std::shared_ptr<pangenome::Graph> pangenome,
    debruijn::Graph& dbg, pangenome::ReadPtr r,
    const std::vector<uint_least32_t>& node_ids, const std::vector<bool>& node_orients)
{
    // suppose tig is subset of read
    // want to remove nodes from position pos (where tig starts in read) + dbg.size()
    // to position start + tig.size()-dbg.size()
    // if pos = 0 (ie read does not include start of tig)

    std::pair<uint32_t, uint32_t> pos = r->find_position(node_ids, node_orients);
    BOOST_LOG_TRIVIAL(trace) << "found pos " << pos.first << " " << pos.second;
    auto start_shift = pos.first;
    if (pos.first > 0 or pos.second < r->get_nodes().size() - 1
        or node_ids.size() == r->get_nodes().size())
        start_shift
            += std::max((int)0, (int)pos.second - (int)node_ids.size()) + dbg.size;
    else {
        std::vector<uint_least32_t> sub_node_ids(
            node_ids.begin() + dbg.size, node_ids.end());
        std::vector<bool> sub_node_orients(
            node_orients.begin() + dbg.size, node_orients.end());
        std::pair<uint32_t, uint32_t> sub_pos
            = r->find_position(sub_node_ids, sub_node_orients);
        if (sub_pos.first > 0)
            start_shift = sub_pos.first;
    }

    auto end_shift = pos.second;
    if (pos.first > 0 or pos.second < r->get_nodes().size() - 1
        or node_ids.size() == r->get_nodes().size()) {
        end_shift -= dbg.size - 1;
    } else {
        std::vector<uint_least32_t> sub_node_ids(
            node_ids.begin(), node_ids.end() - dbg.size);
        std::vector<bool> sub_node_orients(
            node_orients.begin(), node_orients.end() - dbg.size);
        std::pair<uint32_t, uint32_t> sub_pos
            = r->find_position(sub_node_ids, sub_node_orients);
        if (sub_pos.second < pos.second)
            end_shift = sub_pos.second + 1;
    }

    BOOST_LOG_TRIVIAL(trace) << "remove nodes " << start_shift << " to " << end_shift;
    auto it = r->get_nodes().begin() + start_shift;
    for (auto shift = start_shift; shift < end_shift; ++shift) {
        if (it == r->get_nodes().end()) {
            break;
        }
        it = pangenome->remove_node_from_read(it, r);
    }
}

// Remove the internal nodes of low coverage unitigs e.g.
// suppose when dbg kmer size is 3, we have a low covg tig
// 012 -> 126 -> 263 -> 634 -> 345
// then we would remove the 3 internal kmers from the dbg
// and node 6 from the pg->
// If the tig is smaller than k+2 long, currently does nothing
void filter_unitigs(std::shared_ptr<pangenome::Graph> pangraph, debruijn::Graph& dbg,
    const uint_least32_t& threshold)
{
    BOOST_LOG_TRIVIAL(debug) << "Filter unitigs using threshold " << threshold;
    std::vector<uint_least32_t> node_ids;
    std::vector<bool> node_orients;
    std::unordered_set<pangenome::ReadPtr> reads_along_tig;
    bool all_reads_tig;

    std::set<std::deque<uint32_t>> unitigs = dbg.get_unitigs();
    BOOST_LOG_TRIVIAL(debug) << "have " << unitigs.size() << " tigs";
    for (auto d : unitigs) {
        // look up the node ids and orientations associated with this node
        dbg_node_ids_to_ids_and_orientations(dbg, d, node_ids, node_orients);

        // collect the reads covering that tig
        find_reads_along_tig(
            dbg, d, pangraph, node_ids, node_orients, reads_along_tig, all_reads_tig);

        // now if the number of reads covering tig falls below threshold, remove the
        // middle nodes of this tig from the reads
        if (reads_along_tig.size() <= threshold) {
            BOOST_LOG_TRIVIAL(trace)
                << "not enough reads, so remove the tig from the reads";
            for (const auto& r : reads_along_tig) {
                remove_middle_nodes_of_tig_from_read(
                    pangraph, dbg, r, node_ids, node_orients);
            }
            // also remove read_ids from each of the corresponding nodes of dbg
            for (uint32_t i = 1; i < d.size() - 1; ++i) {
                for (const auto& r : reads_along_tig) {
                    dbg.remove_read_from_node(r->id, d[i]);
                }
            }
        }
        reads_along_tig.clear();
    }
}

void detangle_pangraph_with_debruijn_graph(
    std::shared_ptr<pangenome::Graph> pangraph, debruijn::Graph& dbg)
{
    BOOST_LOG_TRIVIAL(debug) << "Detangle pangraph with debruijn graph";
    std::vector<uint_least32_t> node_ids;
    std::vector<bool> node_orients;
    bool all_reads_tig;
    std::unordered_set<pangenome::ReadPtr> reads_along_tig;

    std::set<std::deque<uint32_t>> unitigs = dbg.get_unitigs();
    for (auto d : unitigs) {
        // look up the node ids and orientations associated with this node
        dbg_node_ids_to_ids_and_orientations(dbg, d, node_ids, node_orients);
        // collect the reads covering that tig
        find_reads_along_tig(
            dbg, d, pangraph, node_ids, node_orients, reads_along_tig, all_reads_tig);

        // for each node on tig, for each read covering that node,
        // if we find a read which doesn't lie along whole tig,
        // split that node by reads and create a new node on the tig
        if (!all_reads_tig and !reads_along_tig.empty()) {
            for (uint32_t i = 0; i < node_ids.size(); ++i) {
                for (const auto& r : pangraph->nodes[node_ids[i]]->reads) {
                    if (reads_along_tig.find(r) == reads_along_tig.end()) {
                        pangraph->split_node_by_reads(
                            reads_along_tig, node_ids, node_orients, node_ids[i]);
                        break;
                    }
                }
            }
        }
    }
}

void clean_pangraph_with_debruijn_graph(std::shared_ptr<pangenome::Graph> pangraph,
    const uint_least32_t size, const uint_least32_t threshold, const bool illumina)
{
    BOOST_LOG_TRIVIAL(debug) << "Construct de Bruijn Graph from PanGraph with size "
                             << (uint32_t)size;
    debruijn::Graph dbg(size);
    construct_debruijn_graph(pangraph, dbg);

    if (not illumina)
        remove_leaves(pangraph, dbg, threshold);
    filter_unitigs(pangraph, dbg, threshold);
    BOOST_LOG_TRIVIAL(debug) << "Finished filtering tigs";

    // update dbg now that have removed leaves and some inner nodes
    BOOST_LOG_TRIVIAL(trace) << "Reconstruct dbg";
    construct_debruijn_graph(pangraph, dbg);

    BOOST_LOG_TRIVIAL(trace) << "Now detangle";
    detangle_pangraph_with_debruijn_graph(pangraph, dbg);
}

enum NodeDirection { forward, reverse };

NodeDirection get_pangraph_node_direction(const debruijn::Node& debruijn_node)
{
    const bool forward_node = debruijn_node.hashed_node_ids[0] % 2 != 0;
    if (forward_node)
        return NodeDirection::forward;
    else
        return NodeDirection::reverse;
}

uint32_t get_pangraph_node_id(const debruijn::Node& debruijn_node)
{
    auto direction = get_pangraph_node_direction(debruijn_node);
    if (direction == NodeDirection::forward)
        return (uint32_t)(debruijn_node.hashed_node_ids[0] - 1) / 2;
    else
        return (uint32_t)debruijn_node.hashed_node_ids[0] / 2;
}

namespace gfa {
namespace dump {
    void header(std::ofstream& gfa_fhandle)
    {
        gfa_fhandle << "H\tVN:Z:1.0" << std::endl;
    }

    void nodes(std::shared_ptr<pangenome::Graph> pangraph, std::ofstream& gfa_fhandle)
    {
        for (const auto& node : pangraph->nodes)
            gfa_fhandle << "S\t" << node.second->get_name()
                        << "\tN\tFC:i:" << node.second->covg << std::endl;
    }

    void edge(const pangenome::Node& first_node,
        const NodeDirection& first_node_direction, const pangenome::Node& second_node,
        const NodeDirection& second_node_direction, std::ofstream& gfa_fhandle)
    {
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

pangenome::Node convert_node_debruijn_pangraph(
    const debruijn::Node& debruijn_node, std::shared_ptr<pangenome::Graph> pangraph)
{
    auto node_id = get_pangraph_node_id(debruijn_node);
    const bool node_exists = pangraph->nodes.find(node_id) != pangraph->nodes.end();
    if (!node_exists) {
        fatal_error("Error converting DBG node to pangraph node: the given DBG node "
                    "does not exist in the pangraph");
    }

    auto node_ptr = pangraph->nodes.at(node_id);
    auto node = *node_ptr;
    return node;
}

void write_pangraph_gfa(
    const fs::path& filepath, std::shared_ptr<pangenome::Graph> pangraph)
{
    fs::ofstream gfa_fhandle(filepath);
    gfa::dump::header(gfa_fhandle);
    gfa::dump::nodes(pangraph, gfa_fhandle);

    debruijn::Graph debruijn_graph(1);
    construct_debruijn_graph(pangraph, debruijn_graph);

    for (const auto& debruijn_node_id_entry : debruijn_graph.nodes) {
        auto first_debruijn_node = *debruijn_node_id_entry.second;

        auto first_node = convert_node_debruijn_pangraph(first_debruijn_node, pangraph);
        auto first_node_direction = get_pangraph_node_direction(first_debruijn_node);

        for (const auto& second_debruijn_node_id : first_debruijn_node.out_nodes) {
            const bool neighbour_node_exists_in_the_graph =
                debruijn_graph.nodes.find(second_debruijn_node_id) != debruijn_graph.nodes.end();
            if (!neighbour_node_exists_in_the_graph) {
                fatal_error("Error writing pangraph to GFA: a neighbour of a node does "
                            "not exist in the graph");
            }

            auto& second_debruijn_node = *debruijn_graph.nodes[second_debruijn_node_id];

            auto second_node
                = convert_node_debruijn_pangraph(second_debruijn_node, pangraph);
            auto second_node_direction
                = get_pangraph_node_direction(second_debruijn_node);

            gfa::dump::edge(first_node, first_node_direction, second_node,
                second_node_direction, gfa_fhandle);

            auto erased = second_debruijn_node.out_nodes.erase(first_debruijn_node.id);
            if (erased)
                continue;
            second_debruijn_node.in_nodes.erase(first_debruijn_node.id);
        }
    }
}
