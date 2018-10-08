#include <iostream>
#include <fstream>
#include <cassert>
#include <climits>
#include <unordered_set>
#include <set>
#include <memory>
#include <utility>
#include <algorithm>
#include "pangenome/panread.h"
#include "pangenome/pannode.h"
#include "minihits.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace pangenome;

Read::Read(const uint32_t i) : id(i) {}

void Read::add_hits(const uint32_t prg_id, set<MinimizerHitPtr, pComp> &cluster) {
    if (hits.find(prg_id) == hits.end())
        hits[prg_id] = {};

    auto before_size = hits[prg_id].size();
    hits[prg_id].insert(cluster.begin(), cluster.end());
    assert(hits[prg_id].size() == before_size + cluster.size());
}

// find the index i in the nodes and node_orientations vectors such that [i,i+v.size()]
// corresponds to these vectors of nodes or some vector overlapping end of read
// NB will find the first such instance if there is more than one

// find the position range where overlaps node_ids and node_orients in read
pair<uint32_t, uint32_t>
Read::find_position(const vector<uint_least32_t> &node_ids, const vector<bool> &node_orients,
                    const uint16_t min_overlap) {
    /*cout << "searching for ";
    for (const auto &n : node_ids)
    {
        cout << n << " ";
    }
    cout << " in ";
    for (const auto &n : nodes)
    {
        cout << n->node_id << " ";
    }
    cout << endl;*/

    assert(node_ids.size() == node_orients.size());
    assert(not node_ids.empty());
    uint32_t search_pos = 0;
    uint32_t found_pos = 0;

    //cout << "searching forwards for " << node_ids[0] << " " << node_orients[0];
    //cout << " and backwards for " << node_ids.back() << " " << !node_orients.back() << endl;

    for (uint32_t i = 0; i < nodes.size(); ++i) {
        //cout << "compare nodes pos " << i << " with node_ids 0" << endl;
        // if first node matches at position i going forwards...
        if (nodes[i]->node_id == node_ids[0] and node_orientations[i] == node_orients[0]) {
            //cout << "start node " << i << " fwd " << nodes[i]->node_id << " " << node_orientations[i];
            //cout << " matches " << node_ids[0] << " and " << node_orients[0] << endl;

            search_pos = 0;
            found_pos = 0;
            while (i + found_pos < nodes.size()
                   and nodes[i + found_pos]->node_id == node_ids[search_pos]
                   and node_orientations[i + found_pos] == node_orients[search_pos]) {
                //cout << "fwd " << search_pos << " " << i + found_pos << endl;
                if (search_pos == node_ids.size() - 1 or i + found_pos == nodes.size() - 1) {
                    if (found_pos + 1 >= min_overlap) {
                        return make_pair(i, i + found_pos);
                    } else {
                        break;
                    }
                }
                search_pos++;
                found_pos++;
            }
            //cout << "end fwd" << endl;
        }

        // if i+node_ids.size() is over the end of nodes, consider partial matches which skip the
        // first <size_overhang> nodes
        //cout << "compare nodes pos " << 0 << " with node_ids " << i + node_ids.size() - nodes.size() << endl;
        if (i + node_ids.size() > nodes.size()
            and nodes[0]->node_id == node_ids[i + node_ids.size() - nodes.size()]
            and node_orientations[0] == node_orients[i + node_orients.size() - nodes.size()]) {

            //cout << "start node " << i << " truncated fwd " << nodes[0]->node_id;
            //cout << " " << node_orientations[0];
            //cout << " matches " << node_ids[i + node_ids.size() - nodes.size()];
            //cout << " and " << node_orients[i + node_ids.size() - nodes.size()] << endl;

            search_pos = i + node_ids.size() - nodes.size();
            found_pos = 0;
            while (found_pos < nodes.size()
                   and nodes[found_pos]->node_id == node_ids[search_pos]
                   and node_orientations[found_pos] == node_orients[search_pos]) {
                //cout << "fwd " << search_pos << " " << found_pos << endl;
                if (search_pos == node_ids.size() - 1 or found_pos == nodes.size() - 1) {
                    if (found_pos + 1 >= min_overlap) {
                        return make_pair(0, found_pos);
                    } else {
                        break;
                    }
                }
                search_pos++;
                found_pos++;
            }
        }

        //cout << "compare nodes pos " << nodes.size() -1 -i << " with node_ids " << 0 << endl;
        if (nodes[nodes.size() - 1 - i]->node_id == node_ids[0]
            and node_orientations[node_orientations.size() - 1 - i] == !node_orients[0]) {
            //cout << "start node " << i << " bwd " << nodes[nodes.size() -1 -i]->node_id;
            //cout << " " << node_orientations[nodes.size() -1 -i];
            //cout << " matches " << node_ids[0] << " and " << !node_orients[0] << endl;

            search_pos = 0;
            found_pos = 0;
            while (i + found_pos < nodes.size()
                   and nodes[nodes.size() - 1 - i - found_pos]->node_id == node_ids[search_pos]
                   and node_orientations[nodes.size() - 1 - i - found_pos] == !node_orients[search_pos]) {
                //cout << "bwd " << search_pos << " " << nodes.size() -1 -i -found_pos << endl;
                if (search_pos == node_ids.size() - 1 or i + 1 + found_pos == nodes.size()) {
                    if (found_pos + 1 >= min_overlap) {
                        return make_pair(nodes.size() - 1 - i - found_pos, nodes.size() - 1 - i);
                    } else {
                        break;
                    }
                }
                search_pos++;
                found_pos++;
            }
            //cout << "end bwd" << endl;
        }

        // if we are considering matches which overlap the start backwards, also consider ones which overlap
        // the end backwards by the same amount
        //cout << "compare nodes pos " << nodes.size() -1 << " with node_ids " << nodes.size() -1 - i << endl;
        if (i + node_ids.size() > nodes.size()
            //and nodes[nodes.size() -1 -i]->node_id == node_ids[i + node_ids.size() - nodes.size()]
            //and node_orientations[node_orientations.size() -1  -i] == !node_orients[i + node_orients.size() - nodes.size()])
            and nodes.back()->node_id == node_ids[i + node_ids.size() - nodes.size()]
            and node_orientations.back() == !node_orients[i + node_orients.size() - nodes.size()]) {
            //cout << "start node " << i << " truncated bwd " << nodes.back()->node_id;
            //cout << " " << node_orientations.back();
            // //cout << " matches " << node_ids[node_ids.size() - nodes.size()-1 - i];
            // //cout << " and " << !node_orients[node_ids.size() - nodes.size()-1 - i] << endl;
            //cout << " matches " << node_ids[i + node_ids.size() - nodes.size()];
            //cout << " and " << !node_orients[i + node_ids.size() - nodes.size()] << endl;

            search_pos = i + node_ids.size() - nodes.size();
            found_pos = 0;
            while (found_pos < nodes.size()
                   and nodes[nodes.size() - 1 - found_pos]->node_id == node_ids[search_pos]
                   and node_orientations[nodes.size() - 1 - found_pos] == !node_orients[search_pos]) {
                //cout << "bwd " << search_pos << " " << found_pos << endl;
                if (search_pos == node_ids.size() - 1 or i + 1 + found_pos == nodes.size()) {
                    if (found_pos + 1 >= min_overlap) {
                        return make_pair(nodes.size() - 1 - found_pos, nodes.size() - 1);
                    } else {
                        break;
                    }
                }
                search_pos++;
                found_pos++;
            }
        }
    }
    return make_pair(std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max());
}

void Read::remove_node(NodePtr n_original) {
    //removes all copies of node
    auto it = find(nodes.begin(), nodes.end(), n_original);
    while (it != nodes.end()) {
        uint32_t d = distance(nodes.begin(), it);
        nodes.erase(it);
        node_orientations.erase(node_orientations.begin() + d);
        it = find(nodes.begin(), nodes.end(), n_original);
    }
}

vector<NodePtr>::iterator Read::remove_node(vector<NodePtr>::iterator nit) {
    //(*nit)->covg -= 1;
    uint32_t d = distance(nodes.begin(), nit);
    node_orientations.erase(node_orientations.begin() + d);
    nit = nodes.erase(nit);
    return nit;
}

void Read::replace_node(vector<NodePtr>::iterator n_original, NodePtr n) {
    //hits[n->node_id].insert(hits[(*n_original)->node_id].begin(),hits[(*n_original)->node_id].end() );
    auto it = nodes.erase(n_original);
    nodes.insert(it, n);
    //hits.erase((*n_original)->node_id);

}

bool Read::operator==(const Read &y) const {
    if (id != y.id) { return false; }

    return true;
}

bool Read::operator!=(const Read &y) const {
    return !(*this == y);
}

bool Read::operator<(const Read &y) const {
    return (id < y.id);
}


std::ostream &pangenome::operator<<(std::ostream &out, const pangenome::Read &r) {
    out << r.id << "\t";

    for (const auto i : r.nodes) {
        out << *i << " ";
    }
    out << endl;
    return out;
}

