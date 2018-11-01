#include <sstream>
#include <fstream>
#include <cassert>
#include <limits>
#include <cstdio>      /* NULL */
#include <cstdlib>     /* srand, rand */
#include <cmath>

#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/log/trivial.hpp>

#include "utils.h"
#include "kmernode.h"
#include "kmergraph.h"
#include "localPRG.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


using namespace prg;

KmerGraph::KmerGraph() {
    nodes.reserve(60000);
    reserved_size = 60000;
    num_reads = 0;
    shortest_path_length = 0;
    k = 0; // nb the kmer size is determined by the first non-null node added
    p = 1;
    nb_p = 0.015;
    nb_r = 2;
    thresh = -25;
    exp_depth_covg = 0;
}

// copy constructor
KmerGraph::KmerGraph(const KmerGraph &other) {
    num_reads = other.num_reads;
    shortest_path_length = other.shortest_path_length;
    k = other.k;
    p = other.p;
    nb_p = other.nb_p;
    nb_r = other.nb_r;
    thresh = other.thresh;
    exp_depth_covg = other.exp_depth_covg;
    KmerNodePtr n;

    // first we need to deallocate for any nodes already got!
    clear();
    nodes.reserve(other.nodes.size());

    // create deep copies of the nodes, minus the edges
    for (const auto &node : other.nodes) {
        n = std::make_shared<KmerNode>(*node);
        assert(nodes.size() == n->id);
        nodes.push_back(n);
    }

    // now need to copy the edges
    for (const auto &node : other.nodes) {
        for (uint32_t j = 0; j < node->outNodes.size(); ++j) {
            add_edge(nodes.at(node->id), nodes.at(node->outNodes[j]->id));
        }
    }

    // finally create sorted_nodes
    sort_topologically();

}

// Assignment operator
KmerGraph &KmerGraph::operator=(const KmerGraph &other) {
    // check for self-assignment
    if (this == &other)
        return *this;

    // first we need to deallocate for any nodes already got!
    clear();
    nodes.reserve(other.nodes.size());

    // shallow copy no pointers
    num_reads = other.num_reads;
    shortest_path_length = other.shortest_path_length;
    k = other.k;
    p = other.p;
    nb_p = other.nb_p;
    nb_r = other.nb_r;
    thresh = other.thresh;
    exp_depth_covg = other.exp_depth_covg;
    KmerNodePtr n;

    // create deep copies of the nodes, minus the edges
    for (const auto &node : other.nodes) {
        n = std::make_shared<KmerNode>(*node);
        assert(nodes.size() == n->id);
        nodes.push_back(n);
    }

    // now need to copy the edges
    for (const auto &node : other.nodes) {
        for (uint32_t j = 0; j < node->outNodes.size(); ++j) {
            add_edge(nodes.at(node->id), nodes.at(node->outNodes[j]->id));
        }
    }

    // finally create sorted_nodes
    sort_topologically();

    return *this;
}

KmerGraph::~KmerGraph() {
    clear();
}

void KmerGraph::clear() {
    nodes.clear();
    assert(nodes.empty());

    sorted_nodes.clear();
    assert(sorted_nodes.empty());

    num_reads = 0;
    shortest_path_length = 0;
    k = 0;
    p = 1;
    nb_p = 0.015;
    nb_r = 2;
    thresh = -25;
    exp_depth_covg = 0;
}

KmerNodePtr KmerGraph::add_node(const Path &p) {
    for (const auto &c : nodes) {
        if (c->path == p) {
            return c;
        }
    }

    // if we didn't find an existing node
    KmerNodePtr n(std::make_shared<KmerNode>(nodes.size(), p));
    nodes.push_back(n);
    //nodes[next_id] = make_shared<KmerNode>(next_id, p);
    assert(k == 0 or p.length() == 0 or p.length() == k);
    if (k == 0 and p.length() > 0) {
        k = p.length();
    }
    assert(nodes.size() < std::numeric_limits<uint32_t>::max() ||
           assert_msg("WARNING, reached max kmer graph node size"));
    if (nodes.size() == reserved_size) {
        reserved_size *= 2;
        nodes.reserve(reserved_size);
    }
    return n;
}

KmerNodePtr KmerGraph::add_node_with_kh(const Path &p, const uint64_t &kh, const uint8_t &num) {
    KmerNodePtr n = add_node(p);
    n->khash = kh;
    n->num_AT = num;
    assert(n->khash < std::numeric_limits<uint64_t>::max());
    return n;
}

condition::condition(const Path &p) : q(p) {};

bool condition::operator()(const KmerNodePtr kn) const { return kn->path == q; }

void KmerGraph::add_edge(KmerNodePtr from, KmerNodePtr to) {
    assert(from->id < nodes.size() and nodes[from->id] == from);
    assert(to->id < nodes.size() and nodes[to->id] == to);
    assert(from->path < to->path
           or assert_msg(
            "Cannot add edge from " << from->id << " to " << to->id << " because " << from->path << " is not less than "
                                    << to->path));

    if (find(from->outNodes.begin(), from->outNodes.end(), to) == from->outNodes.end()) {
        from->outNodes.emplace_back(to);
        to->inNodes.emplace_back(from);
        //cout << "added edge from " << from->id << " to " << to->id << endl;
    }
}

void KmerGraph::remove_shortcut_edges() {
    BOOST_LOG_TRIVIAL(debug) << "Remove 'bad' edges from kmergraph";
    Path temp_path;
    uint32_t num_removed_edges = 0;
    std::vector<KmerNodePtr> v = {};
    std::deque<std::vector<KmerNodePtr>> d;

    for (const auto &n : nodes) {
        //cout << n.first << endl;
        for (const auto &out : n->outNodes) {
            for (auto next_out = out->outNodes.begin(); next_out != out->outNodes.end();) {
                // if the outnode of an outnode of A is another outnode of A
                if (find(n->outNodes.begin(), n->outNodes.end(), *next_out) != n->outNodes.end()) {
                    temp_path = get_union(n->path, (*next_out)->path);
                    //cout << "found the union of " << n->path << " and " << (*next_out)->path << endl;
                    //cout << "check if " << temp_path << " contains " << out->path << endl;
                    if (out->path.is_subpath(temp_path)) {
                        //remove it from the outnodes
                        (*next_out)->inNodes.erase(find((*next_out)->inNodes.begin(), (*next_out)->inNodes.end(), out));
                        next_out = out->outNodes.erase(next_out);
                        num_removed_edges += 1;
                        break;
                    } else {
                        next_out++;
                    }
                } else {
                    next_out++;
                }
            }
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "Found and removed " << num_removed_edges << " edges from the kmergraph";
}

void KmerGraph::sort_topologically() {
    sorted_nodes.clear();
    sorted_nodes.reserve(nodes.size());
    sorted_nodes = nodes;
    sort(sorted_nodes.begin(), sorted_nodes.end(), pCompKmerNode());
}

void KmerGraph::check() {
    if (sorted_nodes.empty()) {
        sort_topologically();
    }

    // should not have any leaves, only nodes with degree 0 are start and end
    for (auto c = sorted_nodes.begin(); c != sorted_nodes.end(); ++c) {
        assert(!(*c)->inNodes.empty() or (*c) == sorted_nodes[0] ||
               assert_msg("node" << **c << " has inNodes size " << (*c)->inNodes.size()));
        assert(!(*c)->outNodes.empty() or (*c) == sorted_nodes.back() || assert_msg(
                "node" << **c << " has outNodes size " << (*c)->outNodes.size() << " and isn't equal to back node "
                       << *sorted_nodes.back()));
        for (const auto &d: (*c)->outNodes) {
            assert((*c)->path < d->path || assert_msg((*c)->path << " is not less than " << d->path));
            assert(find(c, sorted_nodes.end(), d) != sorted_nodes.end() ||
                   assert_msg(d->id << " does not occur later in sorted list than " << (*c)->id));
        }
    }
}

void KmerGraph::set_exp_depth_covg(const uint32_t edp) {
    assert(edp > 0);
    exp_depth_covg = edp;
}

void KmerGraph::set_p(const float e_rate) {
    assert(k != 0);
    assert(0 < e_rate and e_rate < 1);
    p = 1 / exp(e_rate * k);
    //cout << "using p: " << p << endl;
}

void KmerGraph::set_nb(const float &nb_prob, const float &nb_fail) {
    if (nb_prob == 0 and nb_fail == 0)
        return;
    //qcout << "set nb" << endl;
    assert((nb_p > 0 and nb_p < 1) || assert_msg("nb_p " << nb_p << " was not set in kmergraph"));
    assert(nb_r > 0 || assert_msg("nb_r was not set in kmergraph"));
    nb_p += nb_prob;
    nb_r += nb_fail;
}

float KmerGraph::nb_prob(uint32_t j) {
    auto k = nodes[j]->covg[0] + nodes[j]->covg[1];
    //cout << "j: " << j << " " << nodes[j]->covg[0] << " + " << nodes[j]->covg[1] << " = "
    //     << nodes[j]->covg[0] + nodes[j]->covg[1] << " = " << k << endl;
    //cout << "nb_prob(" << nb_r << "," << nb_p << "," << k << ") = ";
    float ret = log(pdf(boost::math::negative_binomial(nb_r, nb_p), k));
    ret = std::max(ret, std::numeric_limits<float>::lowest() / 1000);
    //cout << ret << endl;
    return ret;
}

float KmerGraph::prob(uint32_t j) {
    assert(num_reads != 0);
    return prob(j, num_reads);
}

float KmerGraph::prob(uint32_t j, uint32_t num) {
    //prob of node j where j is node id (hence pos in nodes)
    assert(p != 1);
    assert(j < nodes.size());
    if (sorted_nodes.empty() and !nodes.empty()) {
        sort_topologically();
        check();
    }

    //cout << "prob of node " << j << " given " << num << " reads covering and covg " << nodes[j]->covg[0] << " , " << nodes[j]->covg[1];
    float ret;
    if (j == sorted_nodes[0]->id or j == sorted_nodes.back()->id) {
        ret = 0; // is really undefined
    } else if (nodes[j]->covg[0] + nodes[j]->covg[1] > num) {
        // under model assumptions this can't happen, but it inevitably will, so bodge
        ret = lognchoosek2(nodes[j]->covg[0] + nodes[j]->covg[1], nodes[j]->covg[0], nodes[j]->covg[1]) +
              (nodes[j]->covg[0] + nodes[j]->covg[1]) * log(p / 2);
        // note this may give disadvantage to repeat kmers
    } else {
        ret = lognchoosek2(num, nodes[j]->covg[0], nodes[j]->covg[1]) +
              (nodes[j]->covg[0] + nodes[j]->covg[1]) * log(p / 2) +
              (num - (nodes[j]->covg[0] + nodes[j]->covg[1])) * log(1 - p);
    }
    //cout << " is " << ret << endl;
    return ret;
}

float KmerGraph::find_max_path(std::vector<KmerNodePtr> &maxpath) {
    // finds a max likelihood path

    // check if p set
    assert(p < 1 || assert_msg("p was not set in kmergraph"));
    assert(num_reads > 0 || assert_msg("num_reads was not set in kmergraph"));
    //cout << "Kmer graph has " << nodes.size() << " nodes" << endl;

    // need to catch if thesh not set too...

    check();

    // also check not all 0 covgs
    bool not_all_zero = false;
    for (const auto &n : nodes) {
        if (n->covg[0] + n->covg[1] > 0) {
            not_all_zero = true;
            break;
        }
    }
    if (!not_all_zero) {
        BOOST_LOG_TRIVIAL(debug) << "ALL ZEROES";
    }

    // create vectors to hold the intermediate values
    std::vector<float> M(nodes.size(), 0); // max log prob pf paths from pos i to end of graph
    std::vector<int> len(nodes.size(), 0); // length of max log path from pos i to end of graph
    std::vector<uint32_t> prev(nodes.size(), nodes.size() - 1); // prev node along path
    float max_mean;
    int max_len;

    for (uint32_t j = nodes.size() - 1; j != 0; --j) {
        max_mean = std::numeric_limits<float>::lowest();
        max_len = 0; // tie break with longest kmer path
        //cout << "node " << sorted_nodes[j-1]->id << " has " << sorted_nodes[j-1]->outNodes.size() << " outnodes" << endl;
        for (uint32_t i = 0; i != sorted_nodes[j - 1]->outNodes.size(); ++i) {
            if ((sorted_nodes[j - 1]->outNodes[i]->id == sorted_nodes.back()->id and thresh > max_mean + 0.000001) or
                (M[sorted_nodes[j - 1]->outNodes[i]->id] / len[sorted_nodes[j - 1]->outNodes[i]->id] >
                 max_mean + 0.000001) or
                (max_mean - M[sorted_nodes[j - 1]->outNodes[i]->id] / len[sorted_nodes[j - 1]->outNodes[i]->id] <=
                 0.000001 and
                 len[sorted_nodes[j - 1]->outNodes[i]->id] > max_len)) {
                M[sorted_nodes[j - 1]->id] = prob(sorted_nodes[j - 1]->id) + M[sorted_nodes[j - 1]->outNodes[i]->id];
                len[sorted_nodes[j - 1]->id] = 1 + len[sorted_nodes[j - 1]->outNodes[i]->id];
                prev[sorted_nodes[j - 1]->id] = sorted_nodes[j - 1]->outNodes[i]->id;
                //cout << sorted_nodes[j-1]->id << " has prob: " << prob(sorted_nodes[j-1]->id) << "  M: " << M[sorted_nodes[j - 1]->id];
                //cout << " len: " << len[sorted_nodes[j - 1]->id] << " prev: " << prev[sorted_nodes[j - 1]->id] << endl;
                if (sorted_nodes[j - 1]->outNodes[i]->id != sorted_nodes.back()->id) {
                    max_mean = M[sorted_nodes[j - 1]->outNodes[i]->id] / len[sorted_nodes[j - 1]->outNodes[i]->id];
                    max_len = len[sorted_nodes[j - 1]->outNodes[i]->id];
                    //cout << " and new max_mean: " << max_mean;
                } else {
                    max_mean = thresh;
                }
                //cout << endl;
            }
        }
        //cout << sorted_nodes[j-1]->id << "  M: " << M[sorted_nodes[j-1]->id] << " len: " << len[sorted_nodes[j-1]->id] << " prev: " << prev[sorted_nodes[j-1]->id] << endl;
    }
    // remove the final length added for the null start node
    len[0] -= 1;

    // extract path
    uint32_t prev_node = prev[sorted_nodes[0]->id];
    while (prev_node < sorted_nodes.size() - 1) {
        //cout << prev_node << "->";
        maxpath.push_back(nodes[prev_node]);
        prev_node = prev[prev_node];
    }
    //cout << endl;
    //cout << "len[0]: " << len[0] << " maxpath.size(): " << maxpath.size() << " maxpath.back()->id: " << maxpath.back()->id << endl;

    assert(len[0] > 0 || assert_msg("found no path through kmer prg"));
    return M[0] / len[0];
}

float KmerGraph::find_nb_max_path(std::vector<KmerNodePtr> &maxpath) {
    // finds a max likelihood path
    check();

    // also check not all 0 covgs
    bool not_all_zero = false;
    for (const auto &n : nodes) {
        if (n->covg[0] + n->covg[1] > 0) {
            not_all_zero = true;
            break;
        }
    }
    if (!not_all_zero) {
        BOOST_LOG_TRIVIAL(debug) << "ALL ZEROES";
    }

    // create vectors to hold the intermediate values
    std::vector<float> M(nodes.size(), 0); // max log prob pf paths from pos i to end of graph
    std::vector<int> len(nodes.size(), 0); // length of max log path from pos i to end of graph
    std::vector<uint32_t> prev(nodes.size(), nodes.size() - 1); // prev node along path
    float max_mean;
    int max_len;

    for (uint32_t j = nodes.size() - 1; j != 0; --j) {
        max_mean = std::numeric_limits<float>::lowest();
        max_len = 0; // tie break with longest kmer path
        //cout << "node " << sorted_nodes[j-1]->id << " has " << sorted_nodes[j-1]->outNodes.size() << " outnodes" << endl;
        for (uint32_t i = 0; i != sorted_nodes[j - 1]->outNodes.size(); ++i) {
            if ((sorted_nodes[j - 1]->outNodes[i]->id == sorted_nodes.back()->id and thresh > max_mean + 0.000001) or
                (M[sorted_nodes[j - 1]->outNodes[i]->id] / len[sorted_nodes[j - 1]->outNodes[i]->id] >
                 max_mean + 0.000001) or
                (max_mean - M[sorted_nodes[j - 1]->outNodes[i]->id] / len[sorted_nodes[j - 1]->outNodes[i]->id] <=
                 0.000001 and
                 len[sorted_nodes[j - 1]->outNodes[i]->id] > max_len)) {
                M[sorted_nodes[j - 1]->id] = nb_prob(sorted_nodes[j - 1]->id) + M[sorted_nodes[j - 1]->outNodes[i]->id];
                len[sorted_nodes[j - 1]->id] = 1 + len[sorted_nodes[j - 1]->outNodes[i]->id];
                prev[sorted_nodes[j - 1]->id] = sorted_nodes[j - 1]->outNodes[i]->id;
                //cout << sorted_nodes[j-1]->id << " has nb_prob: " << nb_prob(sorted_nodes[j-1]->id) << "  M: " << M[sorted_nodes[j - 1]->id];
                //cout << " len: " << len[sorted_nodes[j - 1]->id] << " prev: " << prev[sorted_nodes[j - 1]->id] << endl;
                if (sorted_nodes[j - 1]->outNodes[i]->id != sorted_nodes.back()->id) {
                    max_mean = M[sorted_nodes[j - 1]->outNodes[i]->id] / len[sorted_nodes[j - 1]->outNodes[i]->id];
                    max_len = len[sorted_nodes[j - 1]->outNodes[i]->id];
                    //cout << " and new max_mean: " << max_mean;
                } else {
                    max_mean = thresh;
                }
                //cout << endl;
            }
        }
        //cout << sorted_nodes[j-1]->id << "  M: " << M[sorted_nodes[j-1]->id] << " len: " << len[sorted_nodes[j-1]->id] << " prev: " << prev[sorted_nodes[j-1]->id] << endl;
    }
    // remove the final length added for the null start node
    len[0] -= 1;

    // extract path
    uint32_t prev_node = prev[sorted_nodes[0]->id];
    while (prev_node < sorted_nodes.size() - 1) {
        //cout << prev_node << "->";
        maxpath.push_back(nodes[prev_node]);
        prev_node = prev[prev_node];
    }
    //cout << endl;
    //cout << "len[0]: " << len[0] << " maxpath.size(): " << maxpath.size() << " maxpath.back()->id: " << maxpath.back()->id << endl;

    assert(len[0] > 0 || assert_msg("found no path through kmer prg"));
    return M[0] / len[0];
}

std::vector<std::vector<KmerNodePtr>> KmerGraph::find_max_paths(uint32_t num) {

    // save original coverges so can put back at the end
    std::vector<uint32_t> original_covgs0, original_covgs1;
    for (uint32_t i = 0; i != nodes.size(); ++i) {
        original_covgs0.push_back(nodes[i]->covg[0]);
        original_covgs1.push_back(nodes[i]->covg[1]);
    }

    // find num max paths
    //cout << "expected covg " << (uint)(p*num_reads/num) << endl;
    std::vector<std::vector<KmerNodePtr>> paths;
    std::vector<KmerNodePtr> maxpath;
    find_max_path(maxpath);
    //uint min_covg;
    paths.push_back(maxpath);

    while (paths.size() < num) {
        for (uint32_t i = 0; i != maxpath.size(); ++i) {
            maxpath[i]->covg[0] -= std::min(maxpath[i]->covg[0], (uint32_t) (p * num_reads / num));
            maxpath[i]->covg[1] -= std::min(maxpath[i]->covg[1], (uint32_t) (p * num_reads / num));
        }
        maxpath.clear();
        find_max_path(maxpath);
        paths.push_back(maxpath);
    }

    // put covgs back
    for (uint32_t i = 0; i != nodes.size(); ++i) {
        nodes[i]->covg[0] = original_covgs0[i];
        nodes[i]->covg[1] = original_covgs1[i];
    }

    return paths;
}

std::vector<std::vector<KmerNodePtr>> KmerGraph::get_random_paths(uint32_t num_paths) {
    // find a random path through kmergraph picking ~uniformly from the outnodes at each point
    std::vector<std::vector<KmerNodePtr>> rpaths;
    std::vector<KmerNodePtr> rpath;
    uint32_t i;

    time_t now;
    now = time(nullptr);
    srand((unsigned int) now);

    if (!nodes.empty()) {
        for (uint32_t j = 0; j != num_paths; ++j) {
            i = rand() % nodes[0]->outNodes.size();
            rpath.push_back(nodes[0]->outNodes[i]);
            while (rpath.back() != nodes[nodes.size() - 1]) {
                if (rpath.back()->outNodes.size() == 1) {
                    rpath.push_back(rpath.back()->outNodes[0]);
                } else {
                    i = rand() % rpath.back()->outNodes.size();
                    rpath.push_back(rpath.back()->outNodes[i]);
                }
            }
            rpath.pop_back();
            rpaths.push_back(rpath);
            rpath.clear();
        }
    }
    return rpaths;
}

float KmerGraph::prob_path(const std::vector<KmerNodePtr> &kpath) {
    float ret_p = 0;
    for (uint32_t i = 0; i != kpath.size(); ++i) {
        ret_p += prob(kpath[i]->id);
    }
    uint32_t len = kpath.size();
    if (kpath[0]->path.length() == 0) {
        len -= 1;
    }
    if (kpath.back()->path.length() == 0) {
        len -= 1;
    }
    if (len == 0) {
        len = 1;
    }
    return ret_p / len;
}

float KmerGraph::prob_paths(const std::vector<std::vector<KmerNodePtr>> &kpaths) {
    if (kpaths.empty()) {
        return 0; // is this the correct default?
    }

    // collect how much coverage we expect on each node from these paths
    std::vector<uint32_t> path_node_covg(nodes.size(), 0);
    for (uint32_t i = 0; i != kpaths.size(); ++i) {
        for (uint32_t j = 0; j != kpaths[i].size(); ++j) {
            path_node_covg[kpaths[i][j]->id] += 1;
        }
    }

    // now calculate max likelihood assuming independent paths
    float ret_p = 0;
    uint32_t len = 0;
    for (uint32_t i = 0; i != path_node_covg.size(); ++i) {
        if (path_node_covg[i] > 0) {
            //cout << "prob of node " << nodes[i]->id << " which has path covg " << path_node_covg[i] << " and so we expect to see " << num_reads*path_node_covg[i]/kpaths.size() << " times IS " << prob(nodes[i]->id, num_reads*path_node_covg[i]/kpaths.size()) << endl;
            ret_p += prob(nodes[i]->id, num_reads * path_node_covg[i] / kpaths.size());
            if (nodes[i]->path.length() > 0) {
                len += 1;
            }
        }
    }

    if (len == 0) {
        len = 1;
    }

    //cout << "len " << len << endl;
    return ret_p / len;
}

void KmerGraph::save_covg_dist(const std::string &filepath) {

    std::ofstream handle;
    handle.open(filepath);

    for (const auto &c : nodes) {
        handle << c->covg[0] << "," << c->covg[1] << "," << (unsigned) c->num_AT << " ";
    }
    handle.close();
    return;
}

uint32_t KmerGraph::min_path_length() {
    if (shortest_path_length > 0) {
        return shortest_path_length;
    }

    if (sorted_nodes.empty()) {
        sort_topologically();
        check();
    }

    std::vector<uint32_t> len(sorted_nodes.size(), 0); // length of shortest path from node i to end of graph
    for (uint32_t j = sorted_nodes.size() - 1; j != 0; --j) {
        for (uint32_t i = 0; i != sorted_nodes[j - 1]->outNodes.size(); ++i) {
            if (len[sorted_nodes[j - 1]->outNodes[i]->id] + 1 > len[j - 1]) {
                len[j - 1] = len[sorted_nodes[j - 1]->outNodes[i]->id] + 1;
            }
        }
    }
    shortest_path_length = len[0];
    return len[0];
}

void KmerGraph::save(const std::string &filepath, const std::shared_ptr<LocalPRG> localprg) {
    std::ofstream handle;
    handle.open(filepath);
    if (handle.is_open()) {
        handle << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << std::endl;
        for (const auto &c : nodes) {
            handle << "S\t" << c->id << "\t";

            if (localprg != nullptr) {
                handle << localprg->string_along_path(c->path);
            } else {
                handle << c->path;
            }

            handle << "\tFC:i:" << c->covg[0] << "\t" << "\tRC:i:"
                   << c->covg[1] << std::endl;//"\t" << (unsigned)nodes[i].second->num_AT << endl;

            for (uint32_t j = 0; j < c->outNodes.size(); ++j) {
                handle << "L\t" << c->id << "\t+\t" << c->outNodes[j]->id << "\t+\t0M" << std::endl;
            }
        }
        handle.close();
    } else {
        std::cerr << "Unable to open kmergraph file " << filepath << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void KmerGraph::load(const std::string &filepath) {

    clear();

    std::string line;
    std::vector<std::string> split_line;
    std::stringstream ss;
    uint32_t id = 0, covg, from, to;
    Path p;
    uint32_t num_nodes = 0;

    std::ifstream myfile(filepath);
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
        std::vector<uint16_t> outnode_counts(num_nodes + 1, 0), innode_counts(num_nodes + 1, 0);

        while (getline(myfile, line).good()) {
            if (line[0] == 'S') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 4);
                id = stoi(split_line[1]);
                ss << split_line[2];
                char c = ss.peek();
                assert (isdigit(c) or assert_msg("Cannot read in this sort of kmergraph GFA as it does not label nodes "
                                                 "with their PRG path"));
                ss >> p;
                ss.clear();
                //add_node(p);
                KmerNodePtr n = std::make_shared<KmerNode>(id, p);
                assert(n != nullptr);
                assert(id == nodes.size() or num_nodes - id == nodes.size() or
                       assert_msg("id " << id << " != " << nodes.size() << " nodes.size() for kmergraph "));
                nodes.push_back(n);
                if (k == 0 and p.length() > 0) {
                    k = p.length();
                }
                covg = stoi(split(split_line[3], "FC:i:")[0]);
                n->covg[0] = covg;
                covg = stoi(split(split_line[4], "RC:i:")[0]);
                n->covg[1] = covg;
                if (split_line.size() >= 6) {
                    n->num_AT = std::stoi(split_line[5]);
                }
            } else if (line[0] == 'L') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 5);
                assert(stoi(split_line[1]) < (int) outnode_counts.size() or
                       assert_msg(stoi(split_line[1]) << ">=" << outnode_counts.size()));
                assert(stoi(split_line[3]) < (int) innode_counts.size() or
                       assert_msg(stoi(split_line[3]) << ">=" << innode_counts.size()));
                outnode_counts[stoi(split_line[1])] += 1;
                innode_counts[stoi(split_line[3])] += 1;
            }
        }

        if (id == 0) {
            reverse(nodes.begin(), nodes.end());
        }

        id = 0;
        for (const auto &n : nodes) {
            assert(nodes[id]->id == id);
            id++;
            assert(n->id < outnode_counts.size() or assert_msg(n->id << ">=" << outnode_counts.size()));
            assert(n->id < innode_counts.size() or assert_msg(n->id << ">=" << innode_counts.size()));
            n->outNodes.reserve(outnode_counts[n->id]);
            n->inNodes.reserve(innode_counts[n->id]);
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
                //nodes[from]->outNodes.push_back(nodes.at(to));
                //nodes[to]->inNodes.push_back(nodes.at(from));
            }
        }
    } else {
        std::cerr << "Unable to open kmergraph file " << filepath << std::endl;
        exit(1);
    }
}

void KmerGraph::append_coverages_to_file(const std::string &filepath, const std::string &sample_name) {
    std::ofstream handle;
    handle.open(filepath, std::ios::app);
    assert (!handle.fail() or assert_msg("Could not open file " << filepath));

    handle << sample_name;
    for (const auto &c : nodes) {
        handle << "\t" << c->covg[0] << "," << c->covg[1];
    }
    handle << std::endl;
    handle.close();
}

bool KmerGraph::operator==(const KmerGraph &y) const {
    // false if have different numbers of nodes
    if (y.nodes.size() != nodes.size()) {//cout << "different numbers of nodes" << endl; 
        return false;
    }

    // false if have different nodes
    for (const auto &c : nodes) {
        // if node not equal to a node in y, then false
        auto found = find_if(y.nodes.begin(), y.nodes.end(), condition(c->path));
        if (found == y.nodes.end()) {
            return false;
        }

        // if the node is found but has different edges, then false
        if (c->outNodes.size() != (*found)->outNodes.size()) { return false; }
        if (c->inNodes.size() != (*found)->inNodes.size()) { return false; }
        for (uint32_t j = 0; j != c->outNodes.size(); ++j) {
            spointer_values_equal<KmerNode> eq = {c->outNodes[j]};
            if (find_if((*found)->outNodes.begin(), (*found)->outNodes.end(), eq) ==
                (*found)->outNodes.end()) { return false; }
        }

    }
    // otherwise is true
    return true;
}

bool pCompKmerNode::operator()(KmerNodePtr lhs, KmerNodePtr rhs) {
    return (lhs->path) < (rhs->path);
}

std::ostream &operator<<(std::ostream &out, KmerGraph const &data) {
    for (const auto &c: data.nodes) {
        out << *(c);
    }
    return out;
}
