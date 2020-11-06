#include <sstream>
#include <cassert>
#include <limits>
#include <cstdlib> /* srand, rand */

#include <boost/log/trivial.hpp>

#include "utils.h"
#include "kmernode.h"
#include "kmergraph.h"
#include "localPRG.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace prg;

KmerGraph::KmerGraph()
{
    shortest_path_length = 0;
    k = 0; // nb the kmer size is determined by the first non-null node added
}

// copy constructor
KmerGraph::KmerGraph(const KmerGraph& other)
{
    shortest_path_length = other.shortest_path_length;
    k = other.k;
    KmerNodePtr n;

    // first we need to deallocate for any nodes already got!
    clear();
    nodes.reserve(other.nodes.size());

    // create deep copies of the nodes, minus the edges
    for (const auto& node : other.nodes) {
        n = std::make_shared<KmerNode>(*node);
        assert(nodes.size() == n->id);
        nodes.push_back(n);
        sorted_nodes.insert(n);
    }

    // now need to copy the edges
    for (const auto& node : other.nodes) {
        for (uint32_t j = 0; j < node->out_nodes.size(); ++j) {
            add_edge(nodes.at(node->id), nodes.at(node->out_nodes[j].lock()->id));
        }
    }
}

// Assignment operator
KmerGraph& KmerGraph::operator=(const KmerGraph& other)
{
    // check for self-assignment
    if (this == &other)
        return *this;

    // first we need to deallocate for any nodes already got!
    clear();
    nodes.reserve(other.nodes.size());

    // shallow copy no pointers
    shortest_path_length = other.shortest_path_length;
    k = other.k;
    KmerNodePtr n;

    // create deep copies of the nodes, minus the edges
    for (const auto& node : other.nodes) {
        n = std::make_shared<KmerNode>(*node);
        assert(nodes.size() == n->id);
        nodes.push_back(n);
        sorted_nodes.insert(n);
    }

    // now need to copy the edges
    for (const auto& node : other.nodes) {
        for (uint32_t j = 0; j < node->out_nodes.size(); ++j) {
            add_edge(nodes.at(node->id), nodes.at(node->out_nodes[j].lock()->id));
        }
    }

    return *this;
}

void KmerGraph::clear()
{
    nodes.clear();
    assert(nodes.empty());

    sorted_nodes.clear();
    assert(sorted_nodes.empty());

    shortest_path_length = 0;
    k = 0;
}

KmerNodePtr KmerGraph::add_node(const prg::Path& p)
{ // add this kmer path to this kmer graph
    for (const auto& c : nodes) { // check if this kmer path is already added
        if (c->path == p) { // TODO: overload operator == to receive a prg::Path?
            return c;
        }
    }

    // if we didn't find an existing node, add this kmer path to the graph
    KmerNodePtr n(std::make_shared<KmerNode>(nodes.size(), p)); // create the node
    nodes.push_back(n); // add it to nodes
    sorted_nodes.insert(n);
    assert(k == 0 or p.length() == 0 or p.length() == k);
    if (k == 0 and p.length() > 0) {
        k = p.length();
    }
    assert(nodes.size() < std::numeric_limits<uint32_t>::max()
        || assert_msg("WARNING, reached max kmer graph node size"));
    if (nodes.size() == reserved_size) {
        reserved_size *= 2;
        nodes.reserve(reserved_size);
    }
    return n;
}

KmerNodePtr KmerGraph::add_node_with_kh(
    const prg::Path& p, const uint64_t& kh, const uint8_t& num)
{
    KmerNodePtr n = add_node(p);
    n->khash = kh;
    n->num_AT = num;
    assert(n->khash < std::numeric_limits<uint64_t>::max());
    return n;
}

condition::condition(const prg::Path& p)
    : q(p) {};

bool condition::operator()(const KmerNodePtr kn) const { return kn->path == q; }

void KmerGraph::add_edge(KmerNodePtr from, KmerNodePtr to)
{
    assert(from->id < nodes.size() and nodes[from->id] == from);
    assert(to->id < nodes.size() and nodes[to->id] == to);
    assert(from->path < to->path
        or assert_msg("Cannot add edge from " << from->id << " to " << to->id
                                              << " because " << from->path
                                              << " is not less than " << to->path));

    if (from->find_node_ptr_in_out_nodes(to) == from->out_nodes.end()) {
        from->out_nodes.emplace_back(to);
        to->in_nodes.emplace_back(from);
    }
}

void KmerGraph::remove_shortcut_edges()
{
    BOOST_LOG_TRIVIAL(debug) << "Remove 'bad' edges from kmergraph";
    prg::Path temp_path;
    uint32_t num_removed_edges = 0;
    std::vector<KmerNodePtr> v = {};
    std::deque<std::vector<KmerNodePtr>> d;

    for (const auto& n : nodes) {
        for (const auto& out : n->out_nodes) {
            auto out_node_as_shared_ptr = out.lock();
            for (auto nextOut = out_node_as_shared_ptr->out_nodes.begin();
                 nextOut != out_node_as_shared_ptr->out_nodes.end();) {
                auto nextOutAsSharedPtr = nextOut->lock();
                // if the outnode of an outnode of A is another outnode of A
                if (n->find_node_ptr_in_out_nodes(nextOutAsSharedPtr)
                    != n->out_nodes.end()) {
                    temp_path = get_union(n->path, nextOutAsSharedPtr->path);

                    if (out_node_as_shared_ptr->path.is_subpath(temp_path)) {
                        // remove it from the outnodes
                        BOOST_LOG_TRIVIAL(debug) << "found the union of " << n->path
                                                 << " and " << nextOutAsSharedPtr->path;
                        BOOST_LOG_TRIVIAL(debug)
                            << "result " << temp_path << " contains "
                            << out_node_as_shared_ptr->path;
                        nextOutAsSharedPtr->in_nodes.erase(
                            nextOutAsSharedPtr->find_node_ptr_in_in_nodes(
                                out_node_as_shared_ptr));
                        nextOut = out_node_as_shared_ptr->out_nodes.erase(nextOut);
                        BOOST_LOG_TRIVIAL(debug)
                            << "next out is now " << nextOutAsSharedPtr->path;
                        num_removed_edges += 1;
                        break;
                    } else {
                        nextOut++;
                    }
                } else {
                    nextOut++;
                }
            }
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "Found and removed " << num_removed_edges
                             << " edges from the kmergraph";
}

void KmerGraph::check() const
{
    // should not have any leaves, only nodes with degree 0 are start and end
    for (auto c = sorted_nodes.begin(); c != sorted_nodes.end(); ++c) {
        assert(!(*c)->in_nodes.empty() or (*c) == (*sorted_nodes.begin())
            || assert_msg(
                "node" << **c << " has inNodes size " << (*c)->in_nodes.size()));
        assert(!(*c)->out_nodes.empty() or (*c) == *(sorted_nodes.rbegin())
            || assert_msg("node"
                << **c << " has outNodes size " << (*c)->out_nodes.size()
                << " and isn't equal to back node " << **(sorted_nodes.rbegin())));
        for (const auto& d : (*c)->out_nodes) {
            auto dAsSharedPtr = d.lock();
            assert((*c)->path < dAsSharedPtr->path
                || assert_msg(
                    (*c)->path << " is not less than " << dAsSharedPtr->path));
            assert(find(c, sorted_nodes.end(), dAsSharedPtr) != sorted_nodes.end()
                || assert_msg(dAsSharedPtr->id
                    << " does not occur later in sorted list than " << (*c)->id));
        }
    }
}

void KmerGraph::discover_k()
{
    if (nodes.size() > 0) {
        auto it = nodes.begin();
        it++;
        const auto& knode = **it;
        k = knode.path.length();
    }
}

std::ostream& operator<<(std::ostream& out, KmerGraph const& data)
{
    for (const auto& c : data.nodes) {
        out << *(c) << std::endl;
    }
    return out;
}

// save the KmerGraph as gfa
void KmerGraph::save(const fs::path& filepath, const std::shared_ptr<LocalPRG> localprg)
{
    fs::ofstream handle;
    handle.open(filepath);
    if (handle.is_open()) {
        handle << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << std::endl;
        for (const auto& c : nodes) {
            handle << "S\t" << c->id << "\t";

            if (localprg != nullptr) {
                handle << localprg->string_along_path(c->path);
            } else {
                handle << c->path;
            }

            // TODO: leave as coverage 0 or change this?
            handle << "\tFC:i:" << 0 << "\t"
                   << "\tRC:i:" << 0 << std::endl;

            for (uint32_t j = 0; j < c->out_nodes.size(); ++j) {
                handle << "L\t" << c->id << "\t+\t" << c->out_nodes[j].lock()->id
                       << "\t+\t0M" << std::endl;
            }
        }
        handle.close();
    } else {
        BOOST_LOG_TRIVIAL(error) << "Unable to open kmergraph file " << filepath;
        std::exit(EXIT_FAILURE);
    }
}

void KmerGraph::load(const fs::path& filepath)
{
    clear();

    std::string line;
    std::vector<std::string> split_line;
    std::stringstream ss;
    uint32_t id = 0, from, to;
    prg::Path p;
    uint32_t num_nodes = 0;

    fs::ifstream myfile(filepath);
    if (myfile.is_open()) {

        while (getline(myfile, line).good()) {
            if (line[0] == 'S') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 4);
                id = std::stoi(split_line[1]);
                num_nodes = std::max(num_nodes, id);
            }
        }
        myfile.clear();
        myfile.seekg(0, myfile.beg);
        nodes.reserve(num_nodes);
        std::vector<uint16_t> outnode_counts(num_nodes + 1, 0),
            innode_counts(num_nodes + 1, 0);

        while (getline(myfile, line).good()) {
            if (line[0] == 'S') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 4);
                id = stoi(split_line[1]);
                ss << split_line[2];
                char c = ss.peek();
                assert(isdigit(c)
                    or assert_msg("Cannot read in this sort of kmergraph GFA as it "
                                  "does not label nodes "
                                  "with their PRG path"));
                ss >> p;
                ss.clear();
                // add_node(p);
                KmerNodePtr kmer_node = std::make_shared<KmerNode>(id, p);
                assert(kmer_node != nullptr);
                assert(id == nodes.size() or num_nodes - id == nodes.size()
                    or assert_msg("id " << id << " != " << nodes.size()
                                        << " nodes.size() for kmergraph "));
                nodes.push_back(kmer_node);
                sorted_nodes.insert(kmer_node);
                if (k == 0 and p.length() > 0) {
                    k = p.length();
                }
                if (split_line.size() >= 6) {
                    kmer_node->num_AT = std::stoi(split_line[5]);
                }
            } else if (line[0] == 'L') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 5);
                assert(stoi(split_line[1]) < (int)outnode_counts.size()
                    or assert_msg(
                        stoi(split_line[1]) << ">=" << outnode_counts.size()));
                assert(stoi(split_line[3]) < (int)innode_counts.size()
                    or assert_msg(stoi(split_line[3]) << ">=" << innode_counts.size()));
                outnode_counts[stoi(split_line[1])] += 1;
                innode_counts[stoi(split_line[3])] += 1;
            }
        }

        if (id == 0) {
            reverse(nodes.begin(), nodes.end());
        }

        id = 0;
        for (const auto& n : nodes) {
            assert(nodes[id]->id == id);
            id++;
            assert(n->id < outnode_counts.size()
                or assert_msg(n->id << ">=" << outnode_counts.size()));
            assert(n->id < innode_counts.size()
                or assert_msg(n->id << ">=" << innode_counts.size()));
            n->out_nodes.reserve(outnode_counts[n->id]);
            n->in_nodes.reserve(innode_counts[n->id]);
        }

        myfile.clear();
        myfile.seekg(0, myfile.beg);

        while (getline(myfile, line).good()) {
            if (line[0] == 'L') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 5);
                if (split_line[2] == split_line[4]) {
                    from = std::stoi(split_line[1]);
                    to = std::stoi(split_line[3]);
                } else {
                    // never happens
                    from = std::stoi(split_line[3]);
                    to = std::stoi(split_line[1]);
                }
                add_edge(nodes[from], nodes[to]);
            }
        }
    } else {
        BOOST_LOG_TRIVIAL(error) << "Unable to open kmergraph file " << filepath;
        exit(1);
    }
}

uint32_t KmerGraph::min_path_length()
{
    // TODO: FIX THIS INNEFICIENCY I INTRODUCED
    std::vector<KmerNodePtr> sorted_nodes(
        this->sorted_nodes.begin(), this->sorted_nodes.end());

    if (shortest_path_length > 0) {
        return shortest_path_length;
    }

#ifndef NDEBUG
    // TODO: check if tests must be updated or not due to this (I think not -
    // sorted_nodes is always sorted) if this is added, some tests bug, since it was not
    // executed before... check();
#endif

    std::vector<uint32_t> len(
        sorted_nodes.size(), 0); // length of shortest path from node i to end of graph
    for (uint32_t j = sorted_nodes.size() - 1; j != 0; --j) {
        for (uint32_t i = 0; i != sorted_nodes[j - 1]->out_nodes.size(); ++i) {
            if (len[sorted_nodes[j - 1]->out_nodes[i].lock()->id] + 1 > len[j - 1]) {
                len[j - 1] = len[sorted_nodes[j - 1]->out_nodes[i].lock()->id] + 1;
            }
        }
    }
    shortest_path_length = len[0];
    return len[0];
}

bool KmerGraph::operator==(const KmerGraph& other_graph) const
{
    // false if have different numbers of nodes
    if (other_graph.nodes.size() != nodes.size()) {
        return false;
    }

    // false if have different nodes
    for (const auto& kmer_node_ptr : nodes) {
        const auto& kmer_node = *kmer_node_ptr;
        // if node not equal to a node in other_graph, then false
        auto found = find_if(other_graph.nodes.begin(), other_graph.nodes.end(),
            condition(kmer_node.path));
        if (found == other_graph.nodes.end()) {
            return false;
        }

        // if the node is found but has different edges, then false
        if (kmer_node.out_nodes.size() != (*found)->out_nodes.size()) {
            return false;
        }
        if (kmer_node.in_nodes.size() != (*found)->in_nodes.size()) {
            return false;
        }
        for (uint32_t j = 0; j != kmer_node.out_nodes.size(); ++j) {
            if ((*found)->find_node_in_out_nodes(*(kmer_node.out_nodes[j].lock()))
                == (*found)->out_nodes.end()) {
                return false;
            }
        }
    }
    return true;
}

bool pCompKmerNode::operator()(KmerNodePtr lhs, KmerNodePtr rhs)
{
    return (lhs->path) < (rhs->path);
}